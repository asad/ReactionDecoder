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
import java.util.Scanner;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import org.openscience.cdk.exception.CDKException;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;
import static java.util.logging.Logger.getLogger;

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
                        Set<IBond> bondChanges = ReactionMappingUtility.getBondChanges(inputReaction);
                        for (IBond bond : bondChanges) {
                            try {
                                formedCleavedWFingerprint.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond), 1.0));
                                formedCleavedWFingerprint.setFingerprintID(reactionID);
                            } catch (CDKException ex) {
                                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                        System.out.println("F/C " + formedCleavedWFingerprint.toString());

                        try {
                            ReactionMechanismTool rmt = new ReactionMechanismTool(inputReaction, true, new StandardizeReaction());
                            MappingSolution s = rmt.getSelectedSolution();
                            System.out.println("Reaction ID: " + s.getReaction().getID() + ", Selected Algorithm: " + s.getAlgorithmID());
                            IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
                            System.out.println("F/C RDT " + s.getBondChangeCalculator().getFormedCleavedWFingerprint().toString());
                            /*
                             * Code for decipt Image generation
                             */
//                            TopToBottomReactionLayoutImage(reactionWithCompressUnChangedHydrogens, (s.getReaction().getID() + s.getAlgorithmID()), "Output");

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

}
