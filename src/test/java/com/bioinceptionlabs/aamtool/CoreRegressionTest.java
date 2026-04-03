/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import static com.bioinceptionlabs.reactionblast.tools.TestUtility.BRENDA_RXN_DIR;
import static com.bioinceptionlabs.reactionblast.tools.TestUtility.BUG_RXN_DIR;
import static com.bioinceptionlabs.reactionblast.tools.TestUtility.KEGG_RXN_DIR;
import static com.bioinceptionlabs.reactionblast.tools.TestUtility.RHEA_RXN_DIR;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.openscience.cdk.interfaces.IReaction;

import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.MappingUtility;

public class CoreRegressionTest extends MappingUtility {

    @Test
    public void identityReactionStaysNoChange() throws Exception {
        IReaction reaction = parseReactionSMILES("CC>>CC");
        reaction.setID("IdentityCore");

        ReactionMechanismTool mapped = getAnnotation(reaction, false);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());

        IPatternFingerprinter formedCleaved = mapped.getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(0, formedCleaved.getFeatureCount());
    }

    @Test
    public void dielsAlderRetainsBondChanges() throws Exception {
        IReaction reaction = parseReactionSMILES("C=CC=C.C=C>>C1CC=CCC1");
        reaction.setID("DielsAlderCore");

        ReactionMechanismTool mapped = getAnnotation(reaction, false);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());

        IPatternFingerprinter formedCleaved = mapped.getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertTrue(formedCleaved.getFeatureCount() > 0);
    }

    @Test
    public void representativeKeggReactionMaps() throws Exception {
        ReactionMechanismTool mapped = testReactions("R01081", KEGG_RXN_DIR);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());
    }

    @Test
    public void representativeRheaReactionMaps() throws Exception {
        ReactionMechanismTool mapped = testReactions("10006", RHEA_RXN_DIR);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());
    }

    @Test
    public void representativeBrendaReactionMaps() throws Exception {
        ReactionMechanismTool mapped = testReactions("ECBLAST_391", BRENDA_RXN_DIR);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());
    }

    @Test
    public void representativeBugCaseStillMaps() throws Exception {
        ReactionMechanismTool mapped = testReactions("R03200", BUG_RXN_DIR);
        assertNotNull(mapped);
        assertNotNull(mapped.getSelectedSolution());
    }
}
