/*
 * RXNFileManipulator.java
 *
 * Created on September 12, 2007, 3:42 AM
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.File;
import java.io.IOException;
import static java.lang.Integer.parseInt;
import static java.lang.System.out;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import static java.util.logging.Level.SEVERE;

import java.util.regex.Pattern;
import static java.util.regex.Pattern.compile;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.BACKWARD;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import static org.openscience.cdk.interfaces.IReaction.Direction.FORWARD;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.containers.FingerPrintContainer;
import uk.ac.ebi.reactionblast.containers.MolContainer;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class RXNFileManipulator extends BasicDebugger {

    private static FingerPrintContainer FPC = FingerPrintContainer.getInstance();
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(RXNFileManipulator.class);
    private MolContainer atomContainer = MolContainer.getInstance();
    private Integer moleculeCounter = 0;

    /**
     *
     *
     */
    public RXNFileManipulator() {
        FPC = FingerPrintContainer.getInstance();
        atomContainer = MolContainer.getInstance();
        moleculeCounter = 0;
        try {
            FPC.Clear();
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        try {
            atomContainer.Clear();
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction process(File rxnFile) throws Exception {

        IReaction cleanedReaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                out.println("Reading the input RXN File");
                String RXNFileName = rxnFile.getCanonicalPath();

                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);

                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);

                ReactionID += fields[0];

//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");
                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();

                out.println("Number of Reactant " + _imoledu.getAtomContainerCount());
                out.println("Number of Product " + _imolpro.getAtomContainerCount());
                cleanedReaction.setID(ReactionID);
                _Stoichiometry.clear();

                Double _TempStoic = 0.0;

                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {

                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }

                }

                setReactantMolecule(cleanedReaction, _Stoichiometry, _Metabolites);

                //System.out.println();
                //System.out.println("****************************");
                //System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;

                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }

                }

                setProductMolecule(cleanedReaction, _Stoichiometry, _Metabolites);
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                cleanedReaction.setDirection(BIDIRECTIONAL);
            }//end of if

        } catch (IOException ex) {
            LOGGER.debug("Error: Generation of FingerPrint Failed the AtomContainer: ");
            ex.printStackTrace();
        }

        FPC.Clear();
        atomContainer.Clear();

        out.println("Reaction is read and initialized");
        return cleanedReaction;

    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processMACiE(File rxnFile) throws Exception {

        IReaction IR = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        if (FileName.endsWith(".rxn")) {
            try {

                String RXNFileName = rxnFile.getCanonicalPath();
//                    molID = null;
                //System.out.println();
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "R";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
//                System.out.println(fields);
                if (!FileName.endsWith(".ov.rxn")) {
                    out.println("Reading the input MACiE rxn (stage reaction) File");
                    ReactionID = ReactionID + fields[0] + fields[1].replace("stg", ".");
                } else {
                    out.println("Reading the input MACiE rxn (overall reaction) File");
                    ReactionID = fields[0] + fields[1].replace("ov", "");
                }
//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");

                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();
//
//                System.out.println(_imoledu.getAtomContainerCount());
//                System.out.println(_imolpro.getAtomContainerCount());

                IR.setID(ReactionID);

                Double _TempStoic = 0.0;
                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
//                System.out.println(_Metabolites.size());
                setReactantMolecule(IR, _Stoichiometry, _Metabolites);
//                System.out.println();
//                System.out.println("****************************");
//                System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setProductMolecule(IR, _Stoichiometry, _Metabolites);
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                IR.setDirection(BIDIRECTIONAL);
                //Print the reaction
//                printReaction(cleanedReaction);
            } catch (IOException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }//end of if

        FPC.Clear();
        atomContainer.Clear();
        out.println("Reaction is read and initialized");
        return IR;
    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processIntEnz(File rxnFile) throws Exception {

        IReaction IR = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                out.println("Reading the input RehA rxn File");
                String RXNFileName = rxnFile.getCanonicalPath();
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
                ReactionID += fields[0];
                int direction = parseInt(fields[1]);
//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");

                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();
                IR.setID(ReactionID);
                _Stoichiometry.clear();
                Double _TempStoic = 0.0;
                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setReactantMolecule(IR, _Stoichiometry, _Metabolites);
                //System.out.println();
                //System.out.println("****************************");
                //System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setProductMolecule(IR, _Stoichiometry, _Metabolites);
                //As per IntEnz 0 for undefined direction, 1 for LR, 2 for RL and 3 for bidirectional
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                switch (direction) {
                    case 1:
                        IR.setDirection(FORWARD);
                        break;
                    case 2:
                        IR.setDirection(BACKWARD);
                        break;
                    case 3:
                    case 0:
                        IR.setDirection(BIDIRECTIONAL);
                        break;
                    default:
                        break;
                }

            } //end of if
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        FPC.Clear();
        atomContainer.Clear();
        out.println("Reaction is read and initialized");
        return IR;
//
    }

    private void setReactantMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {
////
//        System.out.println("-------------------");
//        System.out.println("Reactants");
//        System.out.println("-------------------");

        Iterator<IAtomContainer> it = _Metabolites.iterator();
//        System.out.println("Stoic Map Size: " + _Stoichiometry.size());
        while (it.hasNext()) {
            IAtomContainer mol = it.next();
//            System.out.println("Insertion: " + mol.getID() + ", " + _Stoichiometry.get(mol.getID())
//                    + ", ac:" + mol.getAtomCount());
//            mol.setProperty("STOICHIOMETRY", _Stoichiometry.get(mol.getID()));
            IR.addReactant(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private void setProductMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {

//        System.out.println("-------------------");
//        System.out.println("Products");
//        System.out.println("-------------------");
        Iterator<IAtomContainer> it = _Metabolites.iterator();

        while (it.hasNext()) {

            IAtomContainer mol = it.next();
//            System.out.println("Insertion: " + mol.getID() + ", " + _Stoichiometry.get(mol.getID())
//                    + ", ac:" + mol.getAtomCount());
//            mol.setProperty("STOICHIOMETRY", _Stoichiometry.get(mol.getID()));
            IR.addProduct(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private IAtomContainer setProperty(IAtomContainer mol) throws Exception {
        String molID = mol.getID();
//         System.out.println("After Cleanup");
//        System.out.println(CDKChemaxonIOConveter.getChemAxonMolecule(mol).exportToFormat("mol:V2"));

        IAtomContainer molecule = new AtomContainer(mol);
        removeHydrogensExceptSingleAndPreserveAtomID(mol);
        percieveAtomTypesAndConfigureAtoms(molecule);
        aromatizeMolecule(molecule);
        molecule.setID(molID);

        try {
            try {
                if (molecule.getAtomCount() > 0) {
                    IFingerprintGenerator fpr = new FingerprintGenerator();
                    BitSet fingerprint_Present_Mol = fpr.getFingerprint(molecule);
                    /*
                    Single Atom fingerprints
                     */
                    if (fingerprint_Present_Mol.isEmpty()) {
                        Fingerprinter fingerprinter = new Fingerprinter();
                        fingerprint_Present_Mol = fingerprinter.getBitFingerprint(molecule).asBitSet();
                    }
                    //Loop for Unique Mol ID Creation
                    if (!fingerprint_Present_Mol.isEmpty()) {
//                        System.out.println("FP1 " + fingerprint_Present_Mol.cardinality());
                        if (!FPC.isValuePresent(fingerprint_Present_Mol)) {
//                            System.out.println("FP1 " + fingerprint_Present_Mol.cardinality());
                            if (molID == null) {
                                int val = ++moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
//                            System.out.println("FP " + fingerprint_Present_Mol);
//                            System.out.println("ID " + molID);
                            FPC.setValue(molID, fingerprint_Present_Mol);
                            atomContainer.add(molID, molecule);

//                            System.out.println("\n\n 1: AFTER SET ID " + molID);
//                            System.out.println("F ID " + FPC.getFingerPrint(molID));
//                            System.out.println("Mol ID FROM CONT: " + atomContainer.getMoleculeID(molecule));
//                            System.out.println("Mol  " + molecule.getID());
                        } else if (FPC.isValuePresent(fingerprint_Present_Mol)
                                && atomContainer.compareAtomContainer(FPC.getMoleculeID(fingerprint_Present_Mol), molecule)) {
                            if (molID == null) {
                                molID = FPC.getMoleculeID(fingerprint_Present_Mol);
                                molecule.setID(molID);
                            }
//                            System.out.println("\n\n 2: AFTER SET ID " + molID);
//                            System.out.println("F ID " + FPC.getFingerPrint(molID));
//                            System.out.println("Mol ID FROM CONT: " + atomContainer.getMoleculeID(molecule));
//                            System.out.println("Mol  " + molecule.getID());
                        } else {
                            if (molID == null) {
                                int val = ++moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
                            FPC.setValue(molID, fingerprint_Present_Mol);
                            atomContainer.add(molID, molecule);

//                            System.out.println("\n\n 3: AFTER SET ID " + molID);
//                            System.out.println("F ID " + FPC.getFingerPrint(molID));
//                            System.out.println("Mol ID FROM CONT: " + atomContainer.getMoleculeID(molecule));
//                            System.out.println("Mol  " + molecule.getID());
                        }
                    } else {
                        LOGGER.debug("error: Fingerprint can't be generated for this molecule " + SmilesGenerator.generic().create(molecule));
                    }
                } else {
                    LOGGER.debug("error: Mol file should contain atleast one atom! " + SmilesGenerator.generic().create(molecule));
                }
            } catch (CDKException ex) {
                LOGGER.error(ex);
            } catch (Exception ex) {
                LOGGER.error(ex);
            }
            if (molecule.getID() == null) {
                try {
                    throw new CDKException("Mol ID is NULL");
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
//        printAtoms(molecule);
//        System.out.println(molecule.getID());
        return molecule;
    }
}
