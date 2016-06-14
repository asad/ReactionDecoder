/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import static java.lang.System.out;
import java.net.URISyntaxException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.BOND_CHANGE_INFORMATION;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool;
import static uk.ac.ebi.reactionblast.tools.ImageGenerator.TopToBottomReactionLayoutImage;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import org.openscience.cdk.exception.CDKException;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCanonicalisedBondChangePattern;

/**
 * Test top 10 USPTO Reactions
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class USPTOTest extends MappingUtility {

    private static final Logger LOG = getLogger(USPTOTest.class.getName());

    public USPTOTest() {
        super();
    }

    /**
     * @param args the command line arguments
     * @throws java.net.URISyntaxException
     */
    public static void main(String[] args) throws URISyntaxException {
        // TODO code application logic here
        String path = new SMILES2AAMTest().resource.toURI().getPath();
        File dir = new File(path, USPTO_RXN);
        try {
            out.println("Root Directory: " + path);
            out.println("USPTO RXN File Directory: " + dir.getCanonicalPath());
        } catch (IOException ex) {
            Logger.getLogger(SMILES2AAMTest.class.getName()).log(Level.SEVERE, null, ex);
        }

        int counter = 1;
        /*
         * Instance of SMILES with AAM
         */
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        if (dir.isDirectory()) {
            String readLine = null;
            File reactionFile = new File(dir, "2001-2013_USPTOapplications_reactionSmiles_feb2014filters.rsmi.zip");

            try (ZipInputStream zipIn = new ZipInputStream(new FileInputStream(reactionFile))) {
                ZipEntry entry = zipIn.getNextEntry();
                Scanner sc = new Scanner(zipIn);
                while (sc.hasNextLine()) {
                    readLine = sc.nextLine();
//                    System.out.println(readLine);
                    String[] split = readLine.split("\\s+");
                    String reactionSMILES = split[0].trim();
                    String reactionID = split[split.length - 1].trim();
                    //System.out.println(reactionID+ ", Parsing Input Reaction SMILES " + reactionSMILES);
                    try {
                        IReaction inputReaction = sp.parseReactionSmiles(reactionSMILES);
                        inputReaction.setID(reactionID);
                        IPatternFingerprinter formedCleavedWFingerprint = new PatternFingerprinter();
                        Set<IBond> bondChanges = getBondChanges(inputReaction);
                        for (IBond bond : bondChanges) {
                            try {
                                formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
                                formedCleavedWFingerprint.setFingerprintID(reactionID);
                            } catch (CDKException ex) {
                                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                        System.out.println("F/C " + formedCleavedWFingerprint.toString());
                        
                        try {
                            ReactionMechanismTool rmt = new ReactionMechanismTool(inputReaction, true, true, false, new StandardizeReaction());
                            MappingSolution s = rmt.getSelectedSolution();
                            System.out.println("Reaction ID: " + s.getReaction().getID() + ", Selected Algorithm: " + s.getAlgorithmID());
                            IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
                            System.out.println("F/C RDT " + s.getBondChangeCalculator().getFormedCleavedWFingerprint().toString());
                            /*
                             * Code for decipt Image generation
                             */
                            TopToBottomReactionLayoutImage(reactionWithCompressUnChangedHydrogens, (s.getReaction().getID() + s.getAlgorithmID()), "Output");

//                            /*
//                             * Depict all 4 mappings
//                             */
//                            for (MappingSolution m : rmt.getAllSolutions()) {
//                                out.println("--------------------------------------");
//                                BondChangeCalculator bcc = m.getBondChangeCalculator();
//                                out.println(m.getAlgorithmID() + ", fp " + bcc.getFormedCleavedWFingerprint().toString());
//                                out.println(m.getAlgorithmID() + ", fp " + bcc.getOrderChangesWFingerprint().toString());
//                                
//                                out.println("BE " + m.getBondEnergySum() + ", Fragment " + m.getTotalFragmentChanges());
//                                new ImageGenerator().drawLeftToRightReactionLayout("Output", bcc.getReactionWithCompressUnChangedHydrogens(), ("Map_" + s.getReaction().getID() + m.getAlgorithmID()));
//                                out.println();
//                                out.println("--------------------------------------");
//                            }
                        } catch (AssertionError | Exception ex) {
                            Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    } catch (InvalidSmilesException ex) {
                        Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    /*
                     * Only top 10 reactions tested
                     */
                    counter++;
                    if (counter > 1) {
                        break;
                    }
                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    static Set<IBond> getBondChanges(IReaction mappedReaction) {

        HashMap<IAtom, IAtom> mappings = new HashMap<>();

        System.out.println("Mapping count: " + mappedReaction.getMappingCount());

        for (IAtomContainer container : ExtReactionManipulatorTool.getAllReactants(mappedReaction).atomContainers()) {
            for (IAtom a : container.atoms()) {
                Integer mappingNumber = a.getProperty(ATOM_ATOM_MAPPING);
                if (mappingNumber != null) {
                    mappings.put(a, null);
                }
            }
        }

        for (IAtomContainer container : ExtReactionManipulatorTool.getAllProducts(mappedReaction).atomContainers()) {
            for (IAtom a : container.atoms()) {
                Integer mappingNumber = a.getProperty(ATOM_ATOM_MAPPING);
                if (mappingNumber != null) {
                    for (IAtom mappedAtom : mappings.keySet()) {
                        Integer storedAtomMappingNumber = mappedAtom.getProperty(ATOM_ATOM_MAPPING);
                        if (mappingNumber.intValue() == storedAtomMappingNumber) {
                            mappings.put(mappedAtom, a);
                            break;
                        }
                    }
                }
            }
        }

        Set<IAtom> unMappedAtomsToBeRemoved = new HashSet<>();
        for (Map.Entry<IAtom, IAtom> map : mappings.entrySet()) {
            if (map.getValue() == null) {
                unMappedAtomsToBeRemoved.add(map.getKey());
            }
        }

        for (IAtom removeAtom : unMappedAtomsToBeRemoved) {
            mappings.remove(removeAtom);
        }

//        for (int i = 0; i < mappedReaction.getMappingCount(); i++) {
//            IMapping mapping = mappedReaction.getMapping(i);
//            IAtom a1 = (IAtom) mapping.getChemObject(0);
//            IAtom a2 = (IAtom) mapping.getChemObject(1);
//            mappings.put(a1, a2);
//        }
        System.out.println("Mapping Size: " + mappings.size());

        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        Set<IBond> reactantsbonds = getBonds(allReactants);
        Set<IBond> productsbonds = getBonds(allProducts);

        Set<IBond> bondChange = detectBondsCleavedAndFormed(reactantsbonds, productsbonds, mappings);

        return bondChange;
    }

    static Set<IBond> getBonds(IAtomContainerSet containers) {
        Set<IBond> bonds = new LinkedHashSet<>();
        for (IAtomContainer container : containers.atomContainers()) {
            IBond[] bondArray = ExtAtomContainerManipulator.getBondArray(container);
            bonds.addAll(Arrays.asList(bondArray));
        }
        return bonds;
    }

    static Set<IBond> detectBondsCleavedAndFormed(Set<IBond> reactantsbonds, Set<IBond> productsbonds, HashMap<IAtom, IAtom> mappings) {
        Set<IBond> bondChange = new LinkedHashSet<>();

        for (IBond rb : reactantsbonds) {
            if (mappings.containsKey(rb.getAtom(0)) && mappings.containsKey(rb.getAtom(1))) {
                boolean bondBroken = true;
                for (IBond pb : productsbonds) {
                    if (pb.contains(mappings.get(rb.getAtom(0))) && pb.contains(mappings.get(rb.getAtom(1)))) {
                        bondBroken = false;
                        break;
                    }
                }
                if (bondBroken) {
                    rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                    rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                    bondChange.add(rb);
                }
            } else if (mappings.containsKey(rb.getAtom(0)) && !mappings.containsKey(rb.getAtom(1))) {
                rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                bondChange.add(rb);
            } else if (!mappings.containsKey(rb.getAtom(0)) && mappings.containsKey(rb.getAtom(1))) {
                rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                bondChange.add(rb);
            }
        }

        for (IBond pb : productsbonds) {
            if (mappings.containsKey(pb.getAtom(0)) && mappings.containsKey(pb.getAtom(1))) {
                boolean bondBroken = true;
                for (IBond rb : reactantsbonds) {
                    if (rb.contains(mappings.get(pb.getAtom(0))) && rb.contains(mappings.get(pb.getAtom(1)))) {
                        bondBroken = false;
                        break;
                    }
                }
                if (bondBroken) {
                    pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    bondChange.add(pb);
                }
            } else if (mappings.containsKey(pb.getAtom(0)) && !mappings.containsKey(pb.getAtom(1))) {
                pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                bondChange.add(pb);
            } else if (!mappings.containsKey(pb.getAtom(0)) && mappings.containsKey(pb.getAtom(1))) {
                pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                bondChange.add(pb);
            }
        }

        return bondChange;
    }

}
