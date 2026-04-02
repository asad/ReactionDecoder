package com.bioinceptionlabs.generic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.bioinception.smsd.core.SearchEngine;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.BaseMapping.Algorithm;

public class SMSDCompatibilityFlagsTest {

    private final SmilesParser smilesParser
            = new SmilesParser(SilentChemObjectBuilder.getInstance());

    @Test
    public void legacyMcsConstructorRemainsConnectedByDefault() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CCO");
        IAtomContainer target = smilesParser.parseSmiles("CCO");

        Isomorphism iso = new Isomorphism(
                query,
                target,
                Algorithm.VFLibMCS,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false));

        assertEquals(3, iso.getFirstAtomMapping().getCount());
    }

    @Test
    public void explicitMcsFlagsCanExcludeTargetAtoms() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CCO");
        IAtomContainer target = smilesParser.parseSmiles("CCO");

        SearchEngine.McsOptions options = new SearchEngine.McsOptions();
        options.excludedTargetAtoms = new LinkedHashSet<>();
        options.excludedTargetAtoms.add(2);

        Isomorphism iso = new Isomorphism(
                query,
                target,
                Algorithm.VFLibMCS,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false),
                options);

        java.lang.reflect.Field field = Isomorphism.class.getDeclaredField("mcsOptions");
        field.setAccessible(true);
        SearchEngine.McsOptions effective = (SearchEngine.McsOptions) field.get(iso);

        assertNotNull(effective);
        assertNotNull(effective.excludedTargetAtoms);
        assertTrue(effective.excludedTargetAtoms.contains(2));
        assertEquals(3, iso.getFirstAtomMapping().getCount());
    }

    @Test
    public void substructureFlagsCanLimitEnumeratedMappings() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CC");
        IAtomContainer target = smilesParser.parseSmiles("CCCC");

        Substructure substructure = new Substructure(
                query,
                target,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false),
                true,
                2,
                5_000L);

        assertTrue(substructure.isSubgraph());
        assertEquals(2, substructure.getAllAtomMapping().size());
        assertNotNull(substructure.getLastSearchResult());
        assertTrue(substructure.getLastSearchResult().exists);
    }

    @Test
    public void equivalentSubstructureMatchesUseDeterministicFirstMapping() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CC");
        IAtomContainer target = smilesParser.parseSmiles("CCC");

        Substructure substructure = new Substructure(
                query,
                target,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false),
                true,
                2,
                5_000L);
        substructure.setChemFilters(true, true, true);

        Map<Integer, Integer> mapping = substructure.getFirstAtomMapping().getMappingsByIndex();
        assertEquals(2, substructure.getAllAtomMapping().size());
        assertEquals(Integer.valueOf(0), mapping.get(0));
        assertEquals(Integer.valueOf(1), mapping.get(1));
    }

    @Test
    public void originalAtomRanksGuideEquivalentMappingSelection() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CC");
        IAtomContainer target = smilesParser.parseSmiles("CCC");

        query.getAtom(0).setProperty("OLD_RANK", 1);
        query.getAtom(1).setProperty("OLD_RANK", 2);
        target.getAtom(0).setProperty("OLD_RANK", 30);
        target.getAtom(1).setProperty("OLD_RANK", 10);
        target.getAtom(2).setProperty("OLD_RANK", 20);

        Substructure substructure = new Substructure(
                query,
                target,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false),
                true,
                2,
                5_000L);
        substructure.setChemFilters(true, true, true);

        Map<Integer, Integer> mapping = substructure.getFirstAtomMapping().getMappingsByIndex();
        assertEquals(2, substructure.getAllAtomMapping().size());
        Set<Integer> mappedTargets = new HashSet<>(mapping.values());
        assertEquals(Set.of(1, 2), mappedTargets);
    }

    @Test
    public void sourceAtomIdsGuideEquivalentMappingSelection() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CC");
        IAtomContainer target = smilesParser.parseSmiles("CCC");

        query.getAtom(0).setProperty("sourceAtomId", "R:0:0");
        query.getAtom(1).setProperty("sourceAtomId", "R:0:1");
        target.getAtom(0).setProperty("sourceAtomId", "P:0:2");
        target.getAtom(1).setProperty("sourceAtomId", "P:0:0");
        target.getAtom(2).setProperty("sourceAtomId", "P:0:1");

        Substructure substructure = new Substructure(
                query,
                target,
                AtomBondMatcher.atomMatcher(false, false),
                AtomBondMatcher.bondMatcher(false, false),
                true,
                2,
                5_000L);
        substructure.setChemFilters(true, true, true);

        Map<Integer, Integer> mapping = substructure.getFirstAtomMapping().getMappingsByIndex();
        assertEquals(2, substructure.getAllAtomMapping().size());
        Set<Integer> mappedTargets = new HashSet<>(mapping.values());
        assertEquals(Set.of(1, 2), mappedTargets);
    }
}
