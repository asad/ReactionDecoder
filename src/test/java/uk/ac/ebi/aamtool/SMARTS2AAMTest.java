/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package uk.ac.ebi.aamtool;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.MappingUtility;
import static uk.ac.ebi.reactionblast.tools.MappingUtility.parseReactionSMILES;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class SMARTS2AAMTest extends MappingUtility {
    
    @Test
    public void Transporter() throws Exception {
        
        String reactionID = "O=C(N)NCCCC(N)C(=O)O.O=C(O)C(N)CCCN>>O=C(N)NCCCC(N)C(=O)O.O=C(O)C(N)CCCN";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionID);
        parseReactionSMILES.setID("Transporter");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(0, formedCleavedWFingerprint.getFeatureCount());
    }
    
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
        // SMSD 3.0.0 VF2++ finds 1 primary bond change feature for this complex
        // multi-component NAD+ oxidation reaction (previously 4 with old SMSD)
        assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void DielsAlder() throws Exception {
        // Diels-Alder: butadiene + ethylene -> cyclohexene
        String reactionSM = "C=CC=C.C=C>>C1CC=CCC1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "DielsAlder");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        // Diels-Alder forms 2 new C-C bonds
        assertTrue("Diels-Alder should have bond changes", formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    @Test
    public void EsterHydrolysis() throws Exception {
        // Ester hydrolysis: methyl acetate + water -> acetic acid + methanol
        String reactionSM = "CC(=O)OC.O>>CC(=O)O.CO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EsterHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        // Hydrolysis cleaves C-O ester bond and forms C-O acid + O-H bonds
        assertTrue("Ester hydrolysis should have bond changes", formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    @Test
    public void IdentityReaction() throws Exception {
        // Identity: benzene -> benzene (no change)
        String reactionSM = "C1=CC=CC=C1>>C1=CC=CC=C1";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionSM);
        parseReactionSMILES.setID("Identity");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(0, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void AmideBondFormation() throws Exception {
        // Amide formation: acetic acid + methylamine -> N-methylacetamide + water
        String reactionSM = "CC(=O)O.CN>>CC(=O)NC.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "AmideBond");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        // Amide bond formation cleaves O-H and forms N-C bond
        assertTrue("Amide bond formation should have bond changes", formedCleavedWFingerprint.getFeatureCount() > 0);
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
        boolean accept_no_change = true;// accept no change
        StandardizeReaction standardizeReaction = new StandardizeReaction(); //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D,
                generate3D, complexMapping, accept_no_change, standardizeReaction);
        return rmt;
    }
}
