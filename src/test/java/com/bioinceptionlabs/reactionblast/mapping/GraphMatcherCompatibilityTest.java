/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.BestMatchContainer;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.HydrogenFreeFingerPrintContainer;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
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
                IMappingAlgorithm.MAX, 0, 0, query, target,
                "test-reaction", IMappingAlgorithm.MAX.name(), 1);
        GraphMatcher.MCSSolution solution = thread.call();

        assertNotNull(solution);
        assertEquals(6, solution.getAtomAtomMapping().getCount());
    }

    @Test
    public void mcsThreadShortCircuitsSingleAtomFragments() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("C");
        IAtomContainer target = smilesParser.parseSmiles("CC");
        query.setID("query");
        target.setID("target");
        relabelAtoms(query, 1);
        relabelAtoms(target, 10);

        GraphMatcher.MCSThread thread = new GraphMatcher.MCSThread(
                IMappingAlgorithm.MAX, 0, 0, query, target,
                "test-reaction", IMappingAlgorithm.MAX.name(), 1);
        GraphMatcher.MCSSolution solution = thread.call();

        assertNotNull(solution);
        assertEquals(1, solution.getAtomAtomMapping().getCount());
    }

    @Test
    public void mcsThreadReturnsEmptyWhenFragmentsShareNoAtoms() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = smilesParser.parseSmiles("N");
        IAtomContainer target = smilesParser.parseSmiles("O");
        query.setID("query");
        target.setID("target");
        relabelAtoms(query, 1);
        relabelAtoms(target, 10);

        GraphMatcher.MCSThread thread = new GraphMatcher.MCSThread(
                IMappingAlgorithm.MAX, 0, 0, query, target,
                "test-reaction", IMappingAlgorithm.MAX.name(), 1);
        GraphMatcher.MCSSolution solution = thread.call();

        assertNotNull(solution);
        assertEquals(0, solution.getAtomAtomMapping().getCount());
    }

    @Test
    public void structuralPairKeyIgnoresOccurrenceIdsButSeparatesModes() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer queryA = smilesParser.parseSmiles("CCO");
        IAtomContainer targetA = smilesParser.parseSmiles("CCO");
        IAtomContainer queryB = smilesParser.parseSmiles("CCO");
        IAtomContainer targetB = smilesParser.parseSmiles("CCO");

        relabelAtoms(queryA, 1);
        relabelAtoms(targetA, 10);
        relabelAtoms(queryB, 100);
        relabelAtoms(targetB, 200);

        String maxKeyA = MappingKeyUtil.buildPairKey(
                queryA, targetA, IMappingAlgorithm.MAX.name(),
                false, true, false, false);
        String maxKeyB = MappingKeyUtil.buildPairKey(
                queryB, targetB, IMappingAlgorithm.MAX.name(),
                false, true, false, false);
        String ringsKey = MappingKeyUtil.buildPairKey(
                queryB, targetB, IMappingAlgorithm.RINGS.name(),
                false, false, true, true);

        assertEquals(maxKeyA, maxKeyB);
        assertNotEquals(maxKeyA, ringsKey);
    }

    @Test
    public void replicateMappingRebindsTemplateOntoDifferentOccurrenceAtomIds() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer representativeQuery = smilesParser.parseSmiles("CCO");
        IAtomContainer representativeTarget = smilesParser.parseSmiles("CCO");
        IAtomContainer occurrenceQuery = smilesParser.parseSmiles("CCO");
        IAtomContainer occurrenceTarget = smilesParser.parseSmiles("CCO");

        relabelAtoms(representativeQuery, 1);
        relabelAtoms(representativeTarget, 10);
        relabelAtoms(occurrenceQuery, 101);
        relabelAtoms(occurrenceTarget, 201);

        representativeQuery.setID("R0");
        representativeTarget.setID("P0");
        occurrenceQuery.setID("R1");
        occurrenceTarget.setID("P1");

        AtomAtomMapping template = new AtomAtomMapping(representativeQuery, representativeTarget);
        for (int atomIndex = 0; atomIndex < representativeQuery.getAtomCount(); atomIndex++) {
            template.put(representativeQuery.getAtom(atomIndex), representativeTarget.getAtom(atomIndex));
        }

        GraphMatcher.MCSSolution templateSolution = new GraphMatcher.MCSSolution(
                0, 0, representativeQuery, representativeTarget, template);

        ReactionContainer reactionContainer = new ReactionContainer();
        reactionContainer.putEduct(0, representativeQuery);
        reactionContainer.putProduct(0, representativeTarget);
        reactionContainer.putEduct(1, occurrenceQuery);
        reactionContainer.putProduct(1, occurrenceTarget);

        Holder holder = new Holder(
                IMappingAlgorithm.MAX,
                "reaction",
                Arrays.asList("R0", "R1"),
                Arrays.asList("P0", "P1"),
                reactionContainer,
                new BestMatchContainer(),
                new HydrogenFreeFingerPrintContainer());

        GraphMatcher.MCSSolution replayed = GraphMatcher.replicateMappingOnContainers(
                holder,
                new GraphMatcher.Combination(1, 1),
                templateSolution);

        assertNotNull(replayed);
        assertEquals(3, replayed.getAtomAtomMapping().getCount());

        Set<String> mappedQueryIds = new HashSet<>();
        Set<String> mappedTargetIds = new HashSet<>();
        replayed.getAtomAtomMapping().getMappingsByAtoms().forEach((queryAtom, targetAtom) -> {
            mappedQueryIds.add(queryAtom.getID());
            mappedTargetIds.add(targetAtom.getID());
        });

        assertEquals(new HashSet<>(Arrays.asList("101", "102", "103")), mappedQueryIds);
        assertEquals(new HashSet<>(Arrays.asList("201", "202", "203")), mappedTargetIds);
    }

    private void relabelAtoms(IAtomContainer container, int startId) {
        int atomId = startId;
        for (IAtom atom : container.atoms()) {
            atom.setID(Integer.toString(atomId++));
        }
    }
}
