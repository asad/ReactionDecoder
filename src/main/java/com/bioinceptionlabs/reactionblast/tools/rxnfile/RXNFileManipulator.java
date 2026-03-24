/*
 * RXNFileManipulator.java
 *
 * Created on September 12, 2007, 3:42 AM
 *
 *
 * @author Syed Asad Rahman, BioInception
 * @contact asad.rahman@bioinceptionlabs.com
 *
 */
package com.bioinceptionlabs.reactionblast.tools.rxnfile;

import java.io.File;
import java.io.IOException;
import static java.lang.Integer.parseInt;
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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.BACKWARD;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import static org.openscience.cdk.interfaces.IReaction.Direction.FORWARD;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import com.bioinceptionlabs.reactionblast.containers.FingerPrintContainer;
import com.bioinceptionlabs.reactionblast.containers.MolContainer;
import com.bioinceptionlabs.reactionblast.fingerprints.FingerprintGenerator;
import com.bioinceptionlabs.reactionblast.fingerprints.interfaces.IFingerprintGenerator;
import com.bioinceptionlabs.reactionblast.tools.BasicDebugger;

/**
 *
 * @author Syed Asad Rahman, BioInception
 * @contact asad.rahman@bioinceptionlabs.com
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

        IReaction cleanedReaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                LOGGER.debug("Reading the input RXN File");
                String RXNFileName = rxnFile.getCanonicalPath();

                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);

                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);

                ReactionID += fields[0];

                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();

                LOGGER.debug("Number of Reactant " + _imoledu.getAtomContainerCount());
                LOGGER.debug("Number of Product " + _imolpro.getAtomContainerCount());
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
            LOGGER.error(ex);
        }

        FPC.Clear();
        atomContainer.Clear();

        LOGGER.debug("Reaction is read and initialized");
        return cleanedReaction;

    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processMACiE(File rxnFile) throws Exception {

        IReaction IR = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        if (FileName.endsWith(".rxn")) {
            try {

                String RXNFileName = rxnFile.getCanonicalPath();
//                    molID = null;
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "R";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
                if (!FileName.endsWith(".ov.rxn")) {
                    LOGGER.debug("Reading the input MACiE rxn (stage reaction) File");
                    ReactionID = ReactionID + fields[0] + fields[1].replace("stg", ".");
                } else {
                    LOGGER.debug("Reading the input MACiE rxn (overall reaction) File");
                    ReactionID = fields[0] + fields[1].replace("ov", "");
                }
                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();

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
                setReactantMolecule(IR, _Stoichiometry, _Metabolites);
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
            } catch (IOException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }//end of if

        FPC.Clear();
        atomContainer.Clear();
        LOGGER.debug("Reaction is read and initialized");
        return IR;
    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processIntEnz(File rxnFile) throws Exception {

        IReaction IR = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                LOGGER.debug("Reading the input RehA rxn File");
                String RXNFileName = rxnFile.getCanonicalPath();
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
                ReactionID += fields[0];
                int direction = parseInt(fields[1]);

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
        LOGGER.debug("Reaction is read and initialized");
        return IR;
    }

    private void setReactantMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {

        Iterator<IAtomContainer> it = _Metabolites.iterator();
        while (it.hasNext()) {
            IAtomContainer mol = it.next();
            IR.addReactant(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private void setProductMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {

        Iterator<IAtomContainer> it = _Metabolites.iterator();

        while (it.hasNext()) {

            IAtomContainer mol = it.next();
            IR.addProduct(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private IAtomContainer setProperty(IAtomContainer mol) throws Exception {
        String molID = mol.getID();

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
                        if (!FPC.isValuePresent(fingerprint_Present_Mol)) {
                            if (molID == null) {
                                int val = ++moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
                            FPC.setValue(molID, fingerprint_Present_Mol);
                            atomContainer.add(molID, molecule);

                        } else if (FPC.isValuePresent(fingerprint_Present_Mol)
                                && atomContainer.compareAtomContainer(FPC.getMoleculeID(fingerprint_Present_Mol), molecule)) {
                            if (molID == null) {
                                molID = FPC.getMoleculeID(fingerprint_Present_Mol);
                                molecule.setID(molID);
                            }
                        } else {
                            if (molID == null) {
                                int val = ++moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
                            FPC.setValue(molID, fingerprint_Present_Mol);
                            atomContainer.add(molID, molecule);

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
        return molecule;
    }
}
