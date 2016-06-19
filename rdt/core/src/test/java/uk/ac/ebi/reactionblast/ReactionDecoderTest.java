/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad@ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast;

import static uk.ac.ebi.reactionblast.TestUtility.BRENDA_RXN_DIR;
import java.io.FileNotFoundException;
import static java.lang.System.out;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;
import static uk.ac.ebi.reactionblast.TestUtility.BUG_RXN_DIR;
import static uk.ac.ebi.reactionblast.TestUtility.KEGG_RXN_DIR;
import static uk.ac.ebi.reactionblast.TestUtility.OTHER_RXN;
import static uk.ac.ebi.reactionblast.TestUtility.RHEA_RXN_DIR;
import org.openscience.cdk.interfaces.IAtom;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.smiles.SmilesGenerator;
import uk.ac.ebi.reactionblast.mechanism.interfaces.IBondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;
import static uk.ac.ebi.reactionblast.tools.ReactionSimilarityTool.getSimilarity;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionDecoderTest extends MappingUtility {

    private static final Logger LOG = getLogger(ReactionDecoderTest.class.getName());

    /*
     * Test case for Reaction SMILES
     * Expected Solution
     * MIN, fp 
     * ID=TestReaction:Bond Cleaved and Formed (1)
     * [C%C:2.0]
     *
     * BE 682.0, Fragment 0
     * @throws Exception
     */
    @Test
    public void TestStereo() throws Exception {

        String reactionSM = "C[C@H](N)C(O)=O>>C[C@@H](N)C(O)=O";

        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        parseReactionSmiles.setID("TestReaction");

        ReactionMechanismTool testReactions = getAnnotation(parseReactionSmiles);

        IReaction mappedReaction = smilesParser.parseReactionSmiles(testReactions.getMappedReactionSMILES());
        mappedReaction.setID("TestNewReaction");

        System.out.println("Reaction " + SmilesGenerator.isomeric().aromatic().createReactionSMILES(mappedReaction));

        Map<IAtom, IAtom> mappings = ReactionMappingUtility.getMappings(mappedReaction);
        IPatternFingerprinter stereoWFingerprintNew = new PatternFingerprinter();
        stereoWFingerprintNew.setFingerprintID("TestNewReaction ST");
        Set<IAtom> atomStereoChanges = ReactionMappingUtility.getAtomStereoChanges(mappedReaction, mappings);

        for (IAtom atom : atomStereoChanges) {
            stereoWFingerprintNew.add(new Feature(ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom), 1.0));
        }

        System.out.println("ST " + stereoWFingerprintNew.toString());

        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator().getStereoChangesWFingerprint();
        System.out.println(
                "OLD ST " + stereoChangesWFingerprint);

        assertEquals(1, stereoChangesWFingerprint.getFeatureCount());
    }


    /*
     * Test case for Reaction SMILES
     * Expected Solution
     * MIN, fp 
     * ID=TestReaction:Bond Cleaved and Formed (1)
     * [C%C:2.0]
     *
     * BE 682.0, Fragment 0
     * @throws Exception
     */
    @Test
    public void Test() throws Exception {

        String reactionSM = "CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O";//C[C@H](N)C(O)=O>>C[C@@H](N)C(O)=O

        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        parseReactionSmiles.setID("TestReaction");

        ReactionMechanismTool testReactions = getAnnotation(parseReactionSmiles);
        System.out.println("SM " + testReactions.getMappedReactionSMILES());

        IReaction mappedReaction = smilesParser.parseReactionSmiles(testReactions.getMappedReactionSMILES());
        mappedReaction.setID("TestNewReaction");

        IPatternFingerprinter formedCleavedWFingerprintNew = new PatternFingerprinter();
        formedCleavedWFingerprintNew.setFingerprintID("TestNewReaction FC");
        Map<IAtom, IAtom> mappings = ReactionMappingUtility.getMappings(mappedReaction);
        Set<IBond> bondCleavedFormedChanges = ReactionMappingUtility.getBondCleavedFormedChanges(mappedReaction, mappings);

        System.out.println("Changes NEW COUNT " + bondCleavedFormedChanges.size());

        for (IBond bond : bondCleavedFormedChanges) {
            try {
                formedCleavedWFingerprintNew.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond), 1.0));
            } catch (CDKException ex) {
                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        System.out.println("F/C " + formedCleavedWFingerprintNew.toString());

        IPatternFingerprinter orderWFingerprintNew = new PatternFingerprinter();
        orderWFingerprintNew.setFingerprintID("TestNewReaction FC");

        Map<IBond, IBond> bondOrderChanges = ReactionMappingUtility.getBondOrderChanges(mappedReaction, mappings);

        for (Map.Entry<IBond, IBond> bond : bondOrderChanges.entrySet()) {
            try {
                orderWFingerprintNew.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond.getKey(), bond.getValue()), 1.0));
            } catch (CDKException ex) {
                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        System.out.println("OC " + orderWFingerprintNew.toString());

        IPatternFingerprinter stereoWFingerprintNew = new PatternFingerprinter();
        orderWFingerprintNew.setFingerprintID("TestNewReaction ST");

        Set<IAtom> atomStereoChanges = ReactionMappingUtility.getAtomStereoChanges(mappedReaction, mappings);

        for (IAtom atom : atomStereoChanges) {

            stereoWFingerprintNew.add(new Feature(ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom), 1.0));
        }

        System.out.println("ST " + stereoWFingerprintNew.toString());

        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        System.out.println(
                "OLD F/C " + formedCleavedWFingerprint);
        System.out.println(
                "OLD OC " + testReactions
                .getSelectedSolution()
                .getBondChangeCalculator().getOrderChangesWFingerprint());
        assertEquals(
                1, formedCleavedWFingerprint.getFeatureCount());
        String inputRankLabelledAtomsReactant = testReactions.getSelectedSolution().getReactor().getInputRankLabelledAtomsReactant();

        assertEquals(
                "{\"M00001\": O-3(1), C-2(2), C-4(3), C-5(4), C-1(5)} {\"M00002\": C-5(6), C-4(7), C-3(8), C-2(9), C-1(10)}", inputRankLabelledAtomsReactant);
    }


    /*
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     * Ring rearrangment case
     *
     *
     * @throws Exception
     */
    @Test
    public void R01081() throws Exception {

        String reactionID = "R01081";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R01081:Bond Cleaved and Formed (1)  C%O:2; 
         * BE 706.0, Fragment 0
         */
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     *
     * It should cleave the bonds between the ring systems (bridge)
     *
     * MIN, fp ID=R01557:Bond Cleaved and Formed (2) [C-O:2.0, H-O:2.0]
     *
     * BE 716.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R01557() throws Exception {
        System.gc();
        String reactionID = "R01557";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     *
     * MIN, fp ID=R02007:Bond Cleaved and Formed (3) [C%C:2.0, C-C:2.0, C-H:4.0]
     *
     * BE 1374.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R02007() throws Exception {

        String reactionID = "R02007";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * *******************************************************************
     *
     * SUCESSFUL TEST CASES
     *
     * *******************************************************************
     *
     * @throws java.lang.Exception
     */
    /*
     *@Ring Breaking
     *
     * MIN, fp 
     ID=R01819:Bond Cleaved and Formed (2)
     [C%O:2.0, C-H:2.0]

     BE 706.0, Fragment 0
     *@throws java.lang.Exception
     */
    @Test
    public void R01819() throws Exception {

        String reactionID = "R01819";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *
     *
     * MAX, fp 
     * ID=R09087:Bond Cleaved and Formed (2)
     * [H-O:2.0, O-P:2.0]
     * 
     * BE 670.0, Fragment 0
     *
     * @throws java.lang.Exception
     */
    /**
     *
     * @throws Exception
     */
    @Test
    public void R09087() throws Exception {

        String reactionID = "R09087";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @COMPLEX CASE.
     * RINGS, fp 
     * ID=R05219:Bond Cleaved and Formed (4)
     * [C-C:2.0, C-O:1.0, C-S:1.0, H-O:2.0]
     * 
     * BE 1322.0, Fragment 0
     * 
     * @throws java.lang.Exception
     *
     *
     * @throws Exception
     */
    @Test
    public void R05219() throws Exception {

        String reactionID = "R05219";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Deadlock Test
     *
     *  MIN
     *  ID=R00114:Bond Cleaved and Formed (3)
     *  [C-N:2.0, C-O:1.0, C=O:1.0]
     * 
     *  BE 1767.0, Fragment 0
     */
    /**
     *
     * @throws Exception
     */
    @Test
    public void R00114() throws Exception {

        String reactionID = "R00114";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R05217:Bond Cleaved and Formed (5)
     *
     * [C%O:1.0, C-H:1.0, C-O:1.0, H-O:4.0, O=O:1.0]
     *
     * BE 1205.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R05217() throws Exception {

        String reactionID = "R05217";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @Long Job
     *
     * MIN, fp ID=R00097:Bond Cleaved and Formed (1) [Co-O:2.0]
     *
     * BE 0.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00097() throws Exception {

        String reactionID = "R00097";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @Flips
     * NAP and NAPH mapping is important, check if P flipping takes place
     *
     * @throws Exception
     */
    @Test
    public void R00023() throws Exception {

        String reactionID = "R00023";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @ Carbon
     * NAP and NAPH mapping is important, check if P flipping takes place
     * MIN, fp 
     *ID=R00014:Bond Cleaved and Formed (3)
     *[C-C:2.0, C-H:2.0, H-O:2.0]
     *
     *BE 692.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R00014() throws Exception {

        String reactionID = "R00014";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }


    /*
     * @
     * Phosphate and Water mapping is important, check if P flipping takes place
     * MIN, fp 
     *  ID=R00004:Bond Cleaved and Formed (2)
     *  [H-O:2.0, O-P:2.0]
     * 
     * BE 670.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void TestR00004_time() throws Exception {

        ReactionMechanismTool testReactions = testReactions("R00004", KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *  MAX, fp 
     *  ID=R00025:Bond Cleaved and Formed (5)
     *  [C=N:1.0, C=O:1.0, H-N:2.0, H-O:3.0, O=O:1.0]
     * 
     *  BE 1908.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R00025() throws Exception {

        String reactionID = "R00025";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R03627:Bond Order Changed (1)
     * 
     *  BE 0.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R03627() throws Exception {
        String reactionID = "R03627";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();
        assertEquals(1, orderChangesWFingerprint.getFeatureCount());
    }

    /**
     * @Ring Breaking and Assimilation
     *
     * MIN, fp ID=R08639:Bond Cleaved and Formed (2) [H-O:2.0, O-P:2.0]
     *
     * BE 670.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R08639() throws Exception {

        ReactionMechanismTool testReactions = testReactions("R08639", KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @Ring Breaking and Assimilation
     *
     * MIN, fp ID=R01715:Bond Cleaved and Formed (2) [C-C:1.0, C-O:1.0]
     *
     * BE 704.0, Fragment 0
     *
     * 2 stereo changes, instead of 3
     *
     * @throws Exception
     */
    @Test
    public void R01715() throws Exception {

        String reactionID = "R01715";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        IPatternFingerprinter stereoFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        assertEquals(1, stereoFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R05137:Bond Cleaved and Formed (2)
     *  [C%O:1.0, H-O:2.0]
     * 
     *  BE 353.0, Fragment 0
     * Ring match change missing in the reactant side
     *
     *
     * @throws Exception
     */
    @Test
    public void R05137() throws Exception {

        String reactionID = "R05137";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Test for detecting Stereo Changes
     * @throws Exception
     */
    @Test
    public void R03673() throws Exception {

        String reactionID = "R03673";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        assertEquals(1, stereoChangesWFingerprint.getFeatureCount());
    }

    /*
     *  MIN, fp 
     *  ID=R01068:Bond Cleaved and Formed (2)
     *  [C%C:1.0, C%O:1.0]
     * 
     *  BE 694.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R01068() throws Exception {

        String reactionID = "R01068";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * 
     *  MIN, fp 
     *  ID=R03959:Bond Cleaved and Formed (3)
     *  [C%O:1.0, C-H:1.0, H-O:1.0]
     * 
     *  BE 353.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R03959() throws Exception {

        String reactionID = "R03959";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R06989:Bond Cleaved and Formed (3)
     * [C%O:1.0, C-H:1.0, H-O:1.0]
     * 
     * BE 353.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R06989() throws Exception {

        String reactionID = "R06989";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R01432:Bond Cleaved and Formed (2)
     * [C%O:2.0, C-H:2.0]
     * 
     *  BE 706.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01432() throws Exception {

        String reactionID = "R01432";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R00307:Bond Cleaved and Formed (3) [C%O:2.0]
     *
     * BE 706.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00307() throws Exception {

        String reactionID = "R00307";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R03063:Bond Cleaved and Formed (3)
     *  [C-H:2.0, C-N:2.0, C=O:2.0]
     * 
     *   BE 2208.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R03063() throws Exception {

        String reactionID = "R03063";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *  
     *   MIN, fp 
     *   ID=R01015:Bond Cleaved and Formed (2)
     *   [C-H:2.0, H-O:2.0]
     * 
     *   BE 0.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01015() throws Exception {

        String reactionID = "R01015";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @COMPLEX Case TO DO IMP reaction to test ring re-arrangement Unmapped
     * atoms present
     *
     * MIN, fp ID=R03165:Bond Cleaved and Formed (5)
     *
     * [C%C:2.0, C-C:1.0, C-H:1.0, C-O:1.0, H-O:1.0]
     *
     * BE 1386.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R03165() throws Exception {
        ReactionMechanismTool testReactions = testReactions("R03165", KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R00344:Bond Cleaved and Formed (5) [C-C:1.0, C-H:1.0, C-O:1.0,
     * H-O:2.0, O-P:2.0]
     *
     * BE 1374.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00344() throws Exception {

        String reactionID = "R00344";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Ring exchange test IMP test
     *
     * MIN, fp ID=R04333:Bond Cleaved and Formed (2) [C-O:2.0, H-O:2.0]
     *
     * BE 716.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R04333() throws Exception {

        String reactionID = "R04333";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R01717:Bond Cleaved and Formed (1)
     *  [C-O:2.0]
     * 
     *  BE 716.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01717() throws Exception {

        String reactionID = "R01717";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *   ID=R04165:Bond Cleaved and Formed (2)
     *   [C-H:2.0, H-O:2.0]
     * 
     *   BE 0.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R04165() throws Exception {

        String reactionID = "R04165";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R01394:Bond Cleaved and Formed (2)
     *  [C-H:2.0, H-O:2.0]
     * 
     *  BE 0.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01394() throws Exception {

        String reactionID = "R01394";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R05069:Bond Cleaved and Formed (1)
     * [C-C:2.0]
     *
     * MIN, fp 
     * ID=R05069:Bond Order Change (1)
     * [C-O*C=O:2.0]
     *
     * BE 692.0, Fragment 0
     * @throws Exception
     */
    @Test
    public void R05069() throws Exception {

        String reactionID = "R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R07322:Bond Cleaved and Formed (2)
     * [C%C:5.0, C-H:2.0]
     * 
     * BE 1705.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R07322() throws Exception {

        String reactionID = "R07322";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * a single C-N bond and form a single C-C bond expected
     * 
     *
     * @throws Exception
     */
    @Test
    public void R03020() throws Exception {

        String reactionID = "R03020";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03020:Bond Cleaved and Formed (4)  C-C:1;  C-H:1;  C-N:1;  H-N:1; 
         * BE 651.0, Fragment 2
         */
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * water phophate bond interaction prefered
     * 
     *
     * @throws Exception
     */
    @Test
    public void R03332() throws Exception {

        String reactionID = "R03332";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03332:Bond Cleaved and Formed (2)  H-O:2;  O-P:2; 
         * BE 670.0, Fragment 2
         */
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Expected Stereo (R/S) changes
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01903() throws Exception {

        String reactionID = "R01903";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();

        assertEquals(1, stereoChangesWFingerprint.getFeatureCount());
    }

    /*
     * Stereo (E/Z)
     * 
     *
     * @throws Exception
     */
    @Test
    public void R04538() throws Exception {

        String reactionID = "R04538";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter fp = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        /*
         * Expected Solution
         * MIN, fp ID=R04538:Bond Cleaved and Formed (0)
         * ID=R04538:Bond Stereo Change (1)  C(E/Z):2; 
         * BE 0.0, Fragment 0
         */
        assertEquals(1, fp.getFeatureCount());
    }

    /*
     * O-P bond changes with water then ATP-ADP 
     * prefered big small with water R03187
     *
     * Expected Solution
     * MIN, fp ID=R03187:Bond Cleaved and Formed (5)  
     * [C%N:1;  C-O:1;  H-N:1;  H-O:3;  O-P:2] 
     * BE 1328.0, Fragment 3
     * 
     *
     * @throws Exception
     */
    @Test
    public void R03187() throws Exception {

        String reactionID = "R03187";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R00011() throws Exception {

        String reactionID = "R00011";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R00011:Bond Cleaved and Formed (2)  H-O:2;  O-O:1; 
         * BE 142.0, Fragment 1
         */
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R03364() throws Exception {

        String reactionID = "R03364";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03364:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
         * BE 665.0, Fragment 1
         */
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @COMPLEX CASE.
     * MAX, fp 
     * ID=R04558:Bond Cleaved and Formed (5)
     * [C-N:1.0, C-O:1.0, C=N:1.0, C=O:1.0, C@N:2.0]
     * MAX, fp 
     * ID=R04558:Bond Order Change (1)
     * [C-C*C@C:1.0]
     * BE 2687.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R04558() throws Exception {

        String reactionID = "R04558";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Should break phosphate C-N with O-C in RP00084 C00026_C00217 main
     * MIXTURE, fp 
     * ID=R01148:Bond Cleaved and Formed (3)
     * [C-H:2.0, C-N:2.0, C=O:2.0]
     * 
     * BE 2208.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01148() throws Exception {

        String reactionID = "R01148";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Expected bond changes
     *
     * R05645	C%O	1 
     * R05645	C(R/S)	10 
     * R05645	C-H	2 
     * R05645	C-O*C=O	1 
     * R05645	H-O	2
     * 
     *
     * @throws Exception
     */
    @Test
    public void R05645() throws Exception {

        String reactionID = "R05645";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R05645:Bond Cleaved and Formed (3)  C%O:1;  C-H:2;  H-O:2; 
         * BE 353.0, Fragment 0
         */
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R06361() throws Exception {

        String reactionID = "R06361";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R06361:Bond Cleaved and Formed (5)  C-N:2;  C-O:1;  C-S:1;  H-O:1;  H-S:1; 
         * BE 1240.0, Fragment 4
         */
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * IMP N -> O exchange
     *
     * @throws Exception
     */
    @Test
    public void R00093() throws Exception {

        String reactionID = "R00093";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R00093:Bond Cleaved and Formed (5)  C-H:2;  C-N:2;  C-O:1;  C=O:1;  H-O:1; 
         * BE 1767.0, Fragment 4 
         */
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R07635() throws Exception {

        String reactionID = "R07635";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R07635:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
         * BE 665.0, Fragment 
         */
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * System.out.println("CASE: Condition 11");
     *
     *
     * @throws Exception
     */
    @Test
    public void R05071() throws Exception {

        String reactionID = "R05071";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIXTURE, fp ID=R05071:Bond Cleaved and Formed (2)  C-C:2;  H-O:2; 
         * BE 692.0, Fragment 2
         */
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Water on product should cleave oxygen attached to the ring
     *
     *
     * @throws Exception
     */
    @Test
    public void R02707() throws Exception {

        String reactionID = "R02707";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R02707:Bond Cleaved and Formed (2)  C-O:2;  H-O:2; 
         * BE 716.0, Fragment 2
         */
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Phophate should be cleaved - attached to the ring
     *
     *
     * @throws Exception
     */
    @Test
    public void R00959() throws Exception {

        String reactionID = "R00959";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R00959:Bond Cleaved and Formed (2)  O-P:2; 
         * BE 670.0, Fragment 2
         */
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @Phophate should be cleaved
     * Expected Solution
     * MAX, fp ID=R00959:Bond Cleaved and Formed (2)  H-O:2; O-P:2; 
     * BE 670.0, Fragment 2
     */
    /**
     *
     * @throws Exception
     */
    @Test
    public void R01518() throws Exception {

        String reactionID = "R01518";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @25405	5.4.99.3 KEGG R05069
     *  
     * MIXTURE, fp 
     * ID=25405:Bond Cleaved and Formed (2)
     * [C-C:2.0, H-O:2.0]
     * BE 692.0, Fragment 0
     *   
     *
     * @throws Exception
     */
    @Test
    public void Rhea_25405() throws Exception {

        String reactionID = "25405";
        ReactionMechanismTool testReactions = testReactions(reactionID, RHEA_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MAX, fp 
     * ID=R07965:Bond Cleaved and Formed (6)
     * [C-H:2.0, C-N:1.0, C=O:1.0, H-N:1.0, H-O:2.0, O=O:1.0]
     * 
     * BE 1598.0, Fragment 0
     * Phosphate mappings
     * 
     *
     * @throws Exception
     */
    @Test
    public void R07965() throws Exception {

        String reactionID = "R07965";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(6, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R01569:Bond Cleaved and Formed (2)
     * [H-O:2.0, O-P:2.0]
     * 
     * BE 670.0, Fragment 0
     * Phosphate mappings
     * 
     *
     * @throws Exception
     */
    @Test
    public void R01569() throws Exception {

        String reactionID = "R01569";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        /*
         Check for the presence of O-P bond
         */
        if (!formedCleavedWFingerprint.getFeatures().contains(new Feature("O-P"))) {
            fail();
        }
    }

    /*
     * MIN, fp 
     * ID=R01569:Bond Cleaved and Formed (2)
     * [H-O:2.0, O-P:2.0]

     * BE 670.0, Fragment 0
     * Phosphate mappings
     * 
     *
     * @throws Exception
     */
    @Test
    public void R02918() throws Exception {

        String reactionID = "R02918";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        /*
         Check for the presence of O-P bond
         */
        if (!formedCleavedWFingerprint.getFeatures().contains(new Feature("O-P"))) {
            fail();
        }
    }

    /*
     * MIN, fp 
     * ID=R02555:Bond Cleaved and Formed (3)
     * [C%O:1.0, C-H:2.0, H-O:1.0]
     *
     * BE 353.0, Fragment 0
     * 
     *
     * @throws Exception
     */
    @Test
    public void R02555() throws Exception {

        String reactionID = "R02555";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R00630() throws Exception {

        String reactionID = "R00630";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R03979() throws Exception {

        String reactionID = "R03979";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R05421() throws Exception {

        String reactionID = "R05421";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R01145() throws Exception {

        String reactionID = "R01145";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void aromaticity() throws Exception {

        String reactionID = "aromaticity_check";
        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());

    }

    /*
     * R01529: Bond change should just report stereo changes
     *
     * @throws FileNotFoundException
     * @throws Exception
     */
    @Test
    public void R01529() throws FileNotFoundException, Exception {

        String reactionID = "R01529";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        assertEquals(1, stereoChangesWFingerprint.getFeatureCount());
    }

    /**
     *
     * @throws FileNotFoundException
     * @throws Exception
     */
    @Test
    public void R204() throws FileNotFoundException, Exception {

        String reactionID = "204";
        IBondChangeCalculator bondChanges = testReactions(reactionID, BRENDA_RXN_DIR, false, false)
                .getSelectedSolution().getBondChangeCalculator();
        IPatternFingerprinter formedCleavedWFingerprint = bondChanges
                .getFormedCleavedWFingerprint();
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());

    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void SIM() throws Exception {

        String reactionID1 = "R01194";
        IBondChangeCalculator testRCReactions1 = testReactions(reactionID1, KEGG_RXN_DIR, false, false)
                .getSelectedSolution().getBondChangeCalculator();
        String reactionID2 = "R01188";
        IBondChangeCalculator testRCReactions2 = testReactions(reactionID2, KEGG_RXN_DIR, false, false)
                .getSelectedSolution().getBondChangeCalculator();

        IPatternFingerprinter fp1 = new PatternFingerprinter();
        fp1.add(testRCReactions1.getFormedCleavedWFingerprint());
        fp1.add(testRCReactions1.getOrderChangesWFingerprint());
        fp1.add(testRCReactions1.getStereoChangesWFingerprint());

        IPatternFingerprinter fp2 = new PatternFingerprinter();
        fp2.add(testRCReactions2.getFormedCleavedWFingerprint());
        fp2.add(testRCReactions2.getOrderChangesWFingerprint());
        fp2.add(testRCReactions2.getStereoChangesWFingerprint());

        double similarityBC = getSimilarity(fp1, fp2);
        double similarityRC = getSimilarity(testRCReactions1.getReactionCenterWFingerprint(), testRCReactions2.getReactionCenterWFingerprint());
        out.println("Reaction Centre SIM: " + Math.round(similarityRC * 100) + "%");
        out.println("Bond Change SIM: " + Math.round(similarityBC * 100) + "%");
    }

    /*
     * Functional group regarrangment
     *
     * @throws Exception
     */
    @Test
    public void R07393() throws Exception {
        String reactionID = "R07393";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR, false, false);
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();

        assertEquals(2, orderChangesWFingerprint.getFeatureCount());
    }

    /*
     * TO DO 
     * Functional group regarrangment
     *
     *
     * @throws Exception
     */
    @Test
    public void R02996() throws Exception {
        String reactionID = "R02996";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp 
         ID=R02996:Bond Cleaved and Formed (2)
         [H-O:2.0, O-S:2.0]
         */
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Test for Non Ring RC should differ from RING RC although atoms are same
     * but aromatic and non aromatic bonds should differ
     *
     * @throws Exception
     */
    @Test
    public void NonRingRC() throws Exception {

        String reactionID = "R02718";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R02718:Bond Cleaved and Formed (3)  C-O:2;  H-O:2;  O-P:2; 
         * BE 1386.0, Fragment 4
         */
        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Test for Non Ring RC should differ from RING RC although atoms are same
     * but aromatic and non aromatic bonds should differ
     *
     * @throws Exception
     */
    @Test
    public void RingRC() throws Exception {

        String reactionID = "R01561";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R01561:Bond Cleaved and Formed (4)  C-N:1;  C-O:1;  H-N:1;  H-O:1; 
         * BE 663.0, Fragment 2
         */
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R08761:Bond Cleaved and Formed (4) [C-H:2.0, C-O:1.0, H-O:3.0,
     * O=O:1.0]
     *
     * BE 852.0, Fragment 0
     *
     * @throws java.lang.Exception
     */
    @Test
    public void R08761() throws Exception {

        String reactionID = "R08761";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R00045:Bond Cleaved and Formed (2)
     *  [H-O:8.0, O=O:1.0]
     * 
     *  BE 494.0, Fragment 0
     * keuekal bond changes mapping misplaced
     *
     *
     * @throws Exception
     */
    @Test
    public void R00045() throws Exception {

        String reactionID = "R00045";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        assertEquals(2, orderChangesWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R08765() throws Exception {

        String reactionID = "R08765";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R08765:Bond Cleaved and Formed (1)  C%C:2; 
         * BE 682.0, Fragment 0
         */
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }
}
