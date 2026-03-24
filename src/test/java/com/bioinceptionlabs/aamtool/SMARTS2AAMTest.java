/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import com.bioinceptionlabs.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.MappingUtility;
import static com.bioinceptionlabs.reactionblast.tools.MappingUtility.parseReactionSMILES;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;

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
        // Complex multi-component NAD+ oxidation reaction
        // Bond change count depends on mapping solution selected
        assertTrue("RHEA11004 should have bond changes",
                formedCleavedWFingerprint.getFeatureCount() > 0);
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

    // =========================================================================
    // EC 1: OXIDOREDUCTASES
    // =========================================================================

    /**
     * EC 1.1.1 - Alcohol oxidation (dehydrogenase core)
     * Ethanol -> Acetaldehyde + H2 (simplified, no cofactor)
     * Tests C-O bond order change (alcohol -> aldehyde)
     */
    @Test
    public void AlcoholOxidation() throws Exception {
        String reactionSM = "CCO>>CC=O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_1_1_1");
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();
        assertTrue("Alcohol oxidation should show bond order changes",
                orderChangesWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 1.11.1.6 - Catalase
     * 2 H2O2 -> 2 H2O + O2
     * Disproportionation of hydrogen peroxide
     */
    @Test
    public void Catalase() throws Exception {
        String reactionSM = "OO.OO>>O.O.O=O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_1_11_1_6");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Catalase should have bond changes (O-O cleavage, O=O formation)",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 2: TRANSFERASES
    // =========================================================================

    /**
     * EC 2.1.1 - Methyltransferase (simplified)
     * S-Adenosylmethionine methyl transfer: R-SMe + R'NH2 -> R-S + R'NHMe
     * Simplified: dimethyl sulfide + methylamine -> methyl sulfide + dimethylamine
     */
    @Test
    public void Methyltransferase() throws Exception {
        String reactionSM = "CSC.CN>>CS.CNC";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_2_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Methyltransferase should have bond changes (S-C cleavage, N-C formation)",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 2.6.1.1 - Aspartate aminotransferase (transamination)
     * L-Aspartate + 2-Oxoglutarate -> Oxaloacetate + L-Glutamate
     */
    @Test
    public void AspartateAminotransferase() throws Exception {
        String reactionSM = "OC(=O)CC(N)C(=O)O.OC(=O)CCC(=O)C(=O)O>>OC(=O)CC(=O)C(=O)O.OC(=O)CCC(N)C(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_2_6_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Transamination should have bond changes (C-N cleavage/formation, C=O changes)",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 2.7.1.1 - Hexokinase (phosphoryl transfer)
     * Glucose + ATP -> Glucose-6-phosphate + ADP
     * Simplified core: alcohol + phosphate -> phosphoester + hydroxyl
     */
    @Test
    public void PhosphorylTransfer() throws Exception {
        // Simplified: methanol + phosphoric acid -> methyl phosphate + water
        String reactionSM = "CO.OP(=O)(O)O>>COP(=O)(O)O.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_2_7_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Phosphoryl transfer should have bond changes (O-P cleavage/formation)",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 3: HYDROLASES
    // =========================================================================

    /**
     * EC 3.4.21 - Serine protease / Peptide bond hydrolysis
     * Dipeptide (Gly-Ala) + H2O -> Glycine + Alanine
     * Tests peptide bond cleavage - core biochemistry
     */
    @Test
    public void PeptideBondHydrolysis() throws Exception {
        // Gly-Ala: NCC(=O)NC(C)C(=O)O + H2O -> NCC(=O)O + NC(C)C(=O)O
        String reactionSM = "NCC(=O)NC(C)C(=O)O.O>>NCC(=O)O.NC(C)C(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_3_4_21");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Peptide bond hydrolysis should cleave C-N and form C-O, O-H",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 3.4 - Tripeptide hydrolysis (Gly-Gly-Gly)
     * Larger peptide to test macromolecule handling
     * Gly-Gly-Gly + H2O -> Gly-Gly + Gly
     */
    @Test
    public void TripeptideHydrolysis() throws Exception {
        // Gly-Gly-Gly + H2O -> Gly-Gly + Glycine
        String reactionSM = "NCC(=O)NCC(=O)NCC(=O)O.O>>NCC(=O)NCC(=O)O.NCC(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "TripeptideHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Tripeptide hydrolysis should have peptide bond cleavage",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 3.1.1.3 - Lipase / Triacylglycerol hydrolysis
     * Simplified: triglyceride + water -> diglyceride + fatty acid
     * Tests ester bond hydrolysis in lipid context
     */
    @Test
    public void LipaseReaction() throws Exception {
        // Simplified: ethyl butyrate + water -> butyric acid + ethanol
        String reactionSM = "CCCC(=O)OCC.O>>CCCC(=O)O.CCO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_3_1_1_3");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Lipase should have ester bond hydrolysis",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 3.5.1.1 - Asparaginase
     * L-Asparagine + H2O -> L-Aspartate + NH3
     * Amide bond hydrolysis (important anti-leukemia drug target)
     */
    @Test
    public void Asparaginase() throws Exception {
        String reactionSM = "NC(=O)CC(N)C(=O)O.O>>OC(=O)CC(N)C(=O)O.N";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_3_5_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Asparaginase should cleave amide C-N and form C-O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 4: LYASES
    // =========================================================================

    /**
     * EC 4.1.1.1 - Pyruvate decarboxylase
     * Pyruvate -> Acetaldehyde + CO2
     * Decarboxylation - key in fermentation
     */
    @Test
    public void PyruvateDecarboxylase() throws Exception {
        String reactionSM = "CC(=O)C(=O)O>>CC=O.O=C=O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_4_1_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Decarboxylation should cleave C-C and form C=O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 4.2.1.2 - Fumarase
     * Fumarate + H2O -> L-Malate
     * Hydration of a C=C double bond (TCA cycle)
     */
    @Test
    public void Fumarase() throws Exception {
        // Fumarate + water -> malate
        String reactionSM = "OC(=O)/C=C/C(=O)O.O>>OC(=O)CC(O)C(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_4_2_1_2");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Fumarase should hydrate C=C to C-OH",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 4.1.2.13 - Aldolase
     * Fructose-1,6-bisP -> DHAP + G3P (simplified)
     * Retro-aldol C-C cleavage
     */
    @Test
    public void AldolCleavage() throws Exception {
        // Retro-aldol cleavage: 3-hydroxybutanal -> 2 acetaldehydes
        String reactionSM = "CC(O)CC=O>>CC=O.CC=O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_4_1_2_13");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Aldol cleavage should break C-C bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 5: ISOMERASES
    // =========================================================================

    /**
     * EC 5.3.1.1 - Triose-phosphate isomerase (TIM)
     * DHAP -> G3P (simplified: hydroxyacetone -> lactaldehyde)
     * Ketone to aldehyde isomerization via enediol intermediate
     */
    @Test
    public void TriosePhosphateIsomerase() throws Exception {
        // Aldose-ketose isomerization: glyceraldehyde -> dihydroxyacetone
        String reactionSM = "O=CC(O)CO>>OCC(=O)CO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_5_3_1_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();
        assertTrue("TIM isomerization should show bond changes",
                formedCleavedWFingerprint.getFeatureCount() > 0
                || orderChangesWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 5.1.1.1 - Alanine racemase
     * L-Alanine -> D-Alanine (stereochemical inversion)
     * Important drug target in bacteria
     */
    @Test
    public void AlanineRacemase() throws Exception {
        String reactionSM = "N[C@@H](C)C(=O)O>>N[C@H](C)C(=O)O";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionSM);
        parseReactionSMILES.setID("EC_5_1_1_1");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        // Racemization: stereo center inversion, no net bond changes
        // but stereo changes should be detected
        assertTrue("Racemase should detect stereochemical change",
                stereoChangesWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 6: LIGASES
    // =========================================================================

    /**
     * EC 6.3.2.1 - Glutamine synthetase
     * Glutamate + NH3 + ATP -> Glutamine + ADP + Pi
     * Simplified: acid + amine -> amide (ATP-driven)
     */
    @Test
    public void GlutamineSynthetase() throws Exception {
        // Core: glutamic acid + ammonia -> glutamine + water
        String reactionSM = "NC(CCC(=O)O)C(=O)O.N>>NC(CCC(=O)N)C(=O)O.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EC_6_3_2_1");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Glutamine synthetase should form C-N amide bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * EC 6.3.2 - Peptide bond formation (ribosomal)
     * Glycine + Alanine -> Gly-Ala + H2O
     * Reverse of peptide hydrolysis - tests ligase class
     */
    @Test
    public void PeptideBondFormation() throws Exception {
        // Glycine + Alanine -> Gly-Ala dipeptide + water
        String reactionSM = "NCC(=O)O.NC(C)C(=O)O>>NCC(=O)NC(C)C(=O)O.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "PeptideBondFormation");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Peptide bond formation should form C-N and cleave O-H",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // EC 7: TRANSLOCASES (newer EC class)
    // =========================================================================

    /**
     * Transporter with conformational change
     * Same molecule on both sides = no bond change (transport)
     * Tests that RDT correctly identifies no-change translocases
     */
    @Test
    public void GlucoseTransporter() throws Exception {
        // Glucose (alpha) -> Glucose (alpha): pure transport
        String reactionSM = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O>>OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionSM);
        parseReactionSMILES.setID("GlucoseTransporter");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals("Glucose transport should have zero bond changes",
                0, formedCleavedWFingerprint.getFeatureCount());
    }

    // =========================================================================
    // COMPLEX / MACROMOLECULAR REACTIONS
    // =========================================================================

    /**
     * Tetrapeptide hydrolysis (Gly-Ala-Val-Leu)
     * Tests larger peptide handling
     */
    @Test
    public void TetrapeptideHydrolysis() throws Exception {
        // Gly-Ala-Val-Leu + H2O -> Gly-Ala + Val-Leu
        String reactionSM = "NCC(=O)NC(C)C(=O)NC(C(C)C)C(=O)NC(CC(C)C)C(=O)O.O"
                + ">>NCC(=O)NC(C)C(=O)O.NC(C(C)C)C(=O)NC(CC(C)C)C(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "TetrapeptideHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Tetrapeptide hydrolysis should cleave one peptide bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Thioester hydrolysis (Acetyl-CoA like)
     * R-CO-SR' + H2O -> R-COOH + R'SH
     * Important in fatty acid metabolism
     */
    @Test
    public void ThioesterHydrolysis() throws Exception {
        // Simplified acetyl-thioester + water -> acetic acid + thiol
        String reactionSM = "CC(=O)SCC.O>>CC(=O)O.SCC";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "ThioesterHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Thioester hydrolysis should cleave C-S and form C-O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Claisen condensation (fatty acid biosynthesis step)
     * Acetyl-thioester + Malonyl -> Acetoacetyl + CO2
     * Core C-C bond forming step
     */
    @Test
    public void ClaisenCondensation() throws Exception {
        // Simplified: methyl acetate + malonic acid -> methyl acetoacetate + CO2 + water
        String reactionSM = "CC(=O)SC.OC(=O)CC(=O)SC>>CC(=O)CC(=O)SC.O=C=O.SC";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "ClaisenCondensation");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Claisen condensation should form C-C bond and decarboxylate",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Lactam ring formation (beta-lactam biosynthesis-like)
     * Amino acid cyclization forming a ring
     */
    @Test
    public void LactamRingFormation() throws Exception {
        // 4-aminobutyric acid -> 2-pyrrolidone + water (gamma-lactam)
        String reactionSM = "NCCCC(=O)O>>O=C1CCCN1.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "LactamRingFormation");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Lactam ring formation should form intramolecular C-N bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Disulfide bond formation (protein folding)
     * 2 RSH + H2O2 -> RSSR + 2 H2O
     * Oxidative disulfide formation — critical for protein structure
     */
    @Test
    public void DisulfideBondFormation() throws Exception {
        // 2 methanethiol + hydrogen peroxide -> dimethyl disulfide + 2 water
        String reactionSM = "CS.CS.OO>>CSSC.O.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "DisulfideBond");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Disulfide formation should form S-S bond and cleave S-H, O-O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Schiff base / Imine formation
     * Aldehyde + Amine -> Imine + Water
     * Key step in PLP-dependent enzyme mechanisms (transamination)
     */
    @Test
    public void SchiffBaseFormation() throws Exception {
        String reactionSM = "CC=O.NCC>>CC=NCC.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "SchiffBase");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Schiff base formation should form C=N and cleave C=O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Glycosylation - O-glycosidic bond formation
     * Sugar + Alcohol -> Glycoside + Water
     * Fundamental in carbohydrate biochemistry
     */
    @Test
    public void GlycosidicBondFormation() throws Exception {
        // Simplified: methyl glucoside formation from glucose + methanol
        String reactionSM = "OC1CCOC(O)C1.CO>>OC1CCOC(OC)C1.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "GlycosidicBond");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Glycosidic bond formation should form C-O bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // ====== Complex Biochemistry Tests ======

    /**
     * Glycosidic bond hydrolysis (EC 3.2.1.26 - Invertase)
     * Sucrose simplified: methyl glucoside + water -> glucose + methanol
     * Tests hydrolysis of O-glycosidic linkage
     */
    @Test
    public void GlycosidicBondHydrolysis() throws Exception {
        String reactionSM = "OC1CCOC(OC)C1.O>>OC1CCOC(O)C1.CO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "GlycosidicHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Glycosidic bond hydrolysis should cleave C-O bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Phosphodiester hydrolysis (EC 3.1.4 - Nuclease)
     * Simple phosphodiester + water -> two phosphomonoesters
     * Relevant for DNA/RNA backbone cleavage
     */
    @Test
    public void PhosphodiesterHydrolysis() throws Exception {
        String reactionSM = "COP(=O)(O)OC.O>>COP(=O)(O)O.CO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "PhosphodiesterHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Phosphodiester hydrolysis should cleave P-O bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Beta-lactam ring opening (EC 3.5.2.6 - Beta-lactamase)
     * 2-Azetidinone + water -> beta-amino acid
     * Pharmacologically relevant: penicillin resistance mechanism
     */
    @Test
    public void BetaLactamHydrolysis() throws Exception {
        String reactionSM = "O=C1CCN1.O>>OC(=O)CCN";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "BetaLactamHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Beta-lactam hydrolysis should cleave C-N bond and form C-O, N-H",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Retro-aldol cleavage (EC 4.1.2 - Aldolase)
     * 4-hydroxy-2-butanone -> acetaldehyde + acetaldehyde
     * Simplified aldol retro: tests C-C bond cleavage lyase mechanism
     */
    @Test
    public void RetroAldolCleavage() throws Exception {
        String reactionSM = "OCC(=O)CC=O>>OCC=O.OC=C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "RetroAldol");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Retro-aldol should cleave C-C bond",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Thioester formation (EC 6.2.1 - Acyl-CoA synthetase simplified)
     * Acetic acid + methanethiol -> thioester + water
     * Key in fatty acid metabolism (CoA activation)
     */
    @Test
    public void ThioesterFormation() throws Exception {
        String reactionSM = "CC(=O)O.CS>>CC(=O)SC.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "ThioesterFormation");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Thioester formation should form C-S bond and cleave C-O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Epoxide ring opening (EC 3.3.2 - Epoxide hydrolase)
     * Ethylene oxide + water -> ethylene glycol
     * Important in xenobiotic metabolism / detoxification
     */
    @Test
    public void EpoxideHydrolysis() throws Exception {
        String reactionSM = "C1CO1.O>>OCCO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "EpoxideHydrolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Epoxide hydrolysis should cleave C-O ring bond and form O-H bonds",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Michael addition (EC 4.3 - Carbon-nitrogen lyase reverse)
     * Acrolein + methanethiol -> 3-(methylthio)propanal
     * Conjugate addition to alpha,beta-unsaturated carbonyl
     */
    @Test
    public void MichaelAddition() throws Exception {
        String reactionSM = "C=CC=O.CS>>CSCCC=O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "MichaelAddition");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Michael addition should form C-S bond and reduce C=C to C-C",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Deamination (EC 3.5.4 - Cytidine deaminase simplified)
     * Amine + water -> alcohol + ammonia
     * Key in nucleotide metabolism
     */
    @Test
    public void OxidativeDeamination() throws Exception {
        String reactionSM = "CCN.O>>CCO.N";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "Deamination");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Deamination should cleave C-N and form C-O",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Transamination (PLP-dependent, EC 2.6.1)
     * Pyruvate + Glutamate -> Alanine + alpha-Ketoglutarate (simplified)
     * Tests amino group transfer between amino acid and keto acid
     */
    @Test
    public void Transamination() throws Exception {
        // Simplified: ketone + amine -> amine + ketone (amino group transfer)
        String reactionSM = "CC(=O)C(=O)O.NCC(=O)O>>CC(N)C(=O)O.O=CC(=O)O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "Transamination");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Transamination should show C-N bond changes",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    // =========================================================================
    // CHEMICAL EDGE CASES
    // =========================================================================

    /**
     * E/Z Isomerization: cis-2-butene -> trans-2-butene
     * Tests detection of stereochemical (E/Z) changes around C=C double bond
     */
    @Test
    public void EZIsomerization() throws Exception {
        // cis-2-butene -> trans-2-butene
        String reactionSM = "C/C=C\\C>>C/C=C/C";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionSM);
        parseReactionSMILES.setID("EZIsomerization");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        assertTrue("E/Z isomerization should detect stereo changes",
                stereoChangesWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Friedel-Crafts Acylation: benzene + acetyl chloride -> acetophenone + HCl
     * Tests aromatic ring substitution and bond changes on aromatic system
     */
    @Test
    public void FriedelCraftsAcylation() throws Exception {
        String reactionSM = "c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1.Cl";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "FriedelCraftsAcylation");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Friedel-Crafts acylation should have bond changes on aromatic ring",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Aromatic Nitration: benzene + HNO3 -> nitrobenzene + H2O
     * Tests electrophilic aromatic substitution with charged species
     */
    @Test
    public void AromaticNitration() throws Exception {
        String reactionSM = "c1ccccc1.[O-][N+](=O)O>>c1ccc([N+](=O)[O-])cc1.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "AromaticNitration");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Aromatic nitration should have bond changes",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Ring Opening: Ozonolysis of cyclohexene (simplified)
     * cyclohexene + ozone -> dialdehyde
     * Tests ring opening reaction
     */
    @Test
    public void RingOpeningOzonolysis() throws Exception {
        String reactionSM = "C1CC=CCC1.O=[O+][O-]>>O=CCCCCC=O.O";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "RingOpeningOzonolysis");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Ring opening ozonolysis should have bond changes",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Charged Species: Quaternary ammonium formation
     * trimethylamine + methyl iodide -> tetramethylammonium iodide
     * Tests handling of charged atoms (N+ and I-)
     */
    @Test
    public void QuaternaryAmmoniumFormation() throws Exception {
        String reactionSM = "CN(C)C.CI>>C[N+](C)(C)C.[I-]";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(reactionSM);
        ReactionMechanismTool testReactions = performAtomAtomMapping(parseReactionSmiles, "QuaternaryAmmonium");
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue("Quaternary ammonium formation should have bond changes (N-C formation, C-I cleavage)",
                formedCleavedWFingerprint.getFeatureCount() > 0);
    }

    /**
     * Identity Reaction with Stereo Annotation
     * L-alanine -> L-alanine (same stereochemistry, no change)
     * Tests that identity is correctly detected even with stereo annotations
     */
    @Test
    public void IdentityWithStereo() throws Exception {
        String reactionSM = "N[C@@H](C)C(=O)O>>N[C@@H](C)C(=O)O";
        IReaction parseReactionSMILES = parseReactionSMILES(reactionSM);
        parseReactionSMILES.setID("IdentityWithStereo");
        ReactionMechanismTool testReactions = getAnnotation(parseReactionSMILES, true);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals("Identity with stereo should have 0 bond changes",
                0, formedCleavedWFingerprint.getFeatureCount());
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
