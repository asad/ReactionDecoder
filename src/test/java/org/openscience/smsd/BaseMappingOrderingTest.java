package org.openscience.smsd;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.BaseMapping.Algorithm;

public class BaseMappingOrderingTest {

    private final SmilesParser smilesParser
            = new SmilesParser(SilentChemObjectBuilder.getInstance());

    @Test
    public void bondCompleteMappingStaysAheadOfSymmetricShortcut() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("C1CCCC1");
        IAtomContainer target = buildSubstitutedCyclopentane(query);

        AtomAtomMapping invalid = new AtomAtomMapping(query, target);
        invalid.put(query.getAtom(0), target.getAtom(5));
        invalid.put(query.getAtom(1), target.getAtom(0));
        invalid.put(query.getAtom(2), target.getAtom(1));
        invalid.put(query.getAtom(3), target.getAtom(2));
        invalid.put(query.getAtom(4), target.getAtom(3));

        AtomAtomMapping valid = new AtomAtomMapping(query, target);
        valid.put(query.getAtom(0), target.getAtom(0));
        valid.put(query.getAtom(1), target.getAtom(1));
        valid.put(query.getAtom(2), target.getAtom(2));
        valid.put(query.getAtom(3), target.getAtom(3));
        valid.put(query.getAtom(4), target.getAtom(4));

        ManualBaseMapping mapping = new ManualBaseMapping(query, target);
        mapping.add(invalid);
        mapping.add(valid);

        assertEquals(query.getBondCount(), mapping.getAllBondMaps().get(0).size());
        assertEquals(valid.getMappingsByIndex(), mapping.getFirstAtomMapping().getMappingsByIndex());
    }

    @Test
    public void isSubgraphChecksAllEquivalentMappings() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("C1CCCC1");
        IAtomContainer target = buildSubstitutedCyclopentane(query);

        Isomorphism iso = new Isomorphism(
                query,
                target,
                Algorithm.VFLibMCS,
                AtomMatcher.forElement(),
                BondMatcher.forOrder());

        iso.clearMaps();

        AtomAtomMapping invalid = new AtomAtomMapping(query, target);
        invalid.put(query.getAtom(0), target.getAtom(5));
        invalid.put(query.getAtom(1), target.getAtom(0));
        invalid.put(query.getAtom(2), target.getAtom(1));
        invalid.put(query.getAtom(3), target.getAtom(2));
        invalid.put(query.getAtom(4), target.getAtom(3));

        AtomAtomMapping valid = new AtomAtomMapping(query, target);
        valid.put(query.getAtom(0), target.getAtom(0));
        valid.put(query.getAtom(1), target.getAtom(1));
        valid.put(query.getAtom(2), target.getAtom(2));
        valid.put(query.getAtom(3), target.getAtom(3));
        valid.put(query.getAtom(4), target.getAtom(4));

        iso.getMCSList().add(invalid);
        iso.getMCSList().add(valid);

        assertTrue(iso.isSubgraph());
    }

    @Test
    public void equivalentMappingsPreferLowerSourceDisplacement() throws Exception {
        IAtomContainer query = smilesParser.parseSmiles("CC");
        IAtomContainer target = smilesParser.parseSmiles("CCC");

        query.getAtom(0).setProperty("sourceAtomId", "R:0:0");
        query.getAtom(1).setProperty("sourceAtomId", "R:0:1");
        target.getAtom(0).setProperty("sourceAtomId", "P:0:2");
        target.getAtom(1).setProperty("sourceAtomId", "P:0:0");
        target.getAtom(2).setProperty("sourceAtomId", "P:0:1");

        AtomAtomMapping farther = new AtomAtomMapping(query, target);
        farther.put(query.getAtom(0), target.getAtom(0));
        farther.put(query.getAtom(1), target.getAtom(1));

        AtomAtomMapping closer = new AtomAtomMapping(query, target);
        closer.put(query.getAtom(0), target.getAtom(1));
        closer.put(query.getAtom(1), target.getAtom(2));

        ManualBaseMapping mapping = new ManualBaseMapping(query, target);
        mapping.add(farther);
        mapping.add(closer);

        assertEquals(closer.getMappingsByIndex(), mapping.getFirstAtomMapping().getMappingsByIndex());
    }

    private static final class ManualBaseMapping extends BaseMapping {

        private ManualBaseMapping(IAtomContainer query, IAtomContainer target) {
            super(query, target, AtomMatcher.forElement(), BondMatcher.forOrder());
        }

        private void add(AtomAtomMapping mapping) {
            getMCSList().add(mapping);
        }
    }

    private IAtomContainer buildSubstitutedCyclopentane(IAtomContainer query)
            throws CloneNotSupportedException {
        IAtomContainer target = query.clone();
        IAtom sideChainCarbon = target.getBuilder().newInstance(IAtom.class, "C");
        target.addAtom(sideChainCarbon);
        target.addBond(0, target.getAtomCount() - 1, IBond.Order.SINGLE);
        return target;
    }
}
