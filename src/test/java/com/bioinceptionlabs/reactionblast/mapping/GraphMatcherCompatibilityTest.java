/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import org.junit.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class GraphMatcherCompatibilityTest {

    @Test
    public void mcsThreadRepairsLegacyAromaticFlags() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("c1ccccc1");
        IAtomContainer target = smilesParser.parseSmiles("c1ccccc1");
        query.setID("query");
        target.setID("target");

        for (IAtom atom : query.atoms()) {
            atom.setIsAromatic(false);
        }
        for (IBond bond : query.bonds()) {
            bond.setIsAromatic(true);
        }
        for (IAtom atom : target.atoms()) {
            atom.setIsAromatic(false);
        }
        for (IBond bond : target.bonds()) {
            bond.setIsAromatic(true);
        }

        GraphMatcher.MCSThread thread = new GraphMatcher.MCSThread(
                IMappingAlgorithm.MAX, 0, 0, query, target);
        GraphMatcher.MCSSolution solution = thread.call();

        assertNotNull(solution);
        assertEquals(6, solution.getAtomAtomMapping().getCount());
    }
}
