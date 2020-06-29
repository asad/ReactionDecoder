/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package uk.ac.bioinception;

import static org.junit.Assert.assertEquals;
import org.junit.Test;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.tools.TestUtility;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class SMILES2AAMTest extends TestUtility {

    @Test
    public void TestRHEA10006() throws Exception {

        String reactionSM = "N#CSCC1=CC=CC=C1>>S=C=NCC1=CC=CC=C1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "RHEA10006");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        //System.out.println("formedCleavedWFingerprint " + formedCleavedWFingerprint);
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void TestRHEA11004() throws Exception {

        String reactionSM = "[H+]."
                + "[H+]."
                + "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O."
                + "O=O."
                + "OC1=C(C=CC=C1)C([O-])=O"
                + ">>OC1=C(O)C=CC=C1."
                + "O=C=O."
                + "[H]O[H]."
                + "NC(=O)C1=CC=C[N+](=C1)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "RHEA11004");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        //System.out.println("formedCleavedWFingerprint " + formedCleavedWFingerprint);
        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @param cdkReaction
     * @param reactionName
     * @return
     * @throws InvalidSmilesException
     * @throws AssertionError
     * @throws Exception
     */
    public ReactionMechanismTool performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
        cdkReaction.setID(reactionName);
        /*
         RMT for the reaction mapping
         */
        boolean forceMapping = true;//Overrides any mapping present int the reaction
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers
        boolean complexMapping = true;//Rings
        StandardizeReaction standardizeReaction = new StandardizeReaction(); //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, complexMapping, standardizeReaction);
        return rmt;
    }
}
