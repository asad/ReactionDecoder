/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.mcs;

import java.util.List;
import static org.junit.Assert.assertEquals;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MCSTest {

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);

        org.openscience.smsd.algorithm.mcsplus1.MCSPlus mcs
                = new org.openscience.smsd.algorithm.mcsplus1.MCSPlus(query, target, atomMatcher, bondMatcher);
        mcs.search_cliques();
        List<List<Integer>> overlaps = mcs.getFinalMappings();
        assertEquals(overlaps.size(), 2);
        assertEquals(7, (mcs.getFinalMappings().iterator().next().size() / 2));
    }

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws Exception
     */
    @Test
    public void testSet_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("*C1OC(CO)C(OP=O)C1O");
        IAtomContainer target = sp.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);
        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);
        AtomAtomMapping firstAtomMapping = mcs.getFirstAtomMapping();
//        System.out.println("mcs " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
        assertEquals(9, (firstAtomMapping.getCount()));
    }

    @Test
    public void test_VFMCS_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("*C1OC(CO)C(OP=O)C1O");
        IAtomContainer target = sp.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");

//        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
//        IAtomContainer target = sp.parseSmiles("OC1=CC=CC=C1");
//        IAtomContainer query = sp.parseSmiles("CC1=C(CCC([O-])=O)C2=[N+]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N+]5[Fe-]3(n14)n1c(=C6)c(C)c(CCC([O-])=O)c1=C2");
//        IAtomContainer target = sp.parseSmiles("C1=2N3C(C=C4[N+]5=C(C=C6N7C8=CC9=[N+](C(=C1)C(=C9CCC([O-])=O)C)[Fe-2]573OC(CCC=%10C=%11C=C%12C(=C(C%13=CC%14=[N+]%15C(=CC=%16N%17C(=C(C%16C=C)C)C=C([N+]%11[Fe-2]%15%17(N%13%12)OC(CCC8=C6C)=O)C%10C)C(=C%14C=C)C)C)CCC([O-])=O)=O)C(=C4C)C=C)=C(C2C)C=C");
//
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

        Isomorphism mcs
                = new Isomorphism(query, target, Algorithm.CDKMCS, atomMatcher, bondMatcher);
//                = new Isomorphism(query, target, Algorithm.VFLibMCS, false, false, false);

//        System.out.println("mcs.getFirstAtomMapping() " + mcs.getFirstAtomMapping());
//        List<AtomAtomMapping> allAtomMapping = mcs.getAllAtomMapping();
        // assertEquals(2, allAtomMapping.size());
//        System.out.println("MCS " + (mcs.getFirstAtomMapping().getCommonFragmentAsSMILES()));
        assertEquals(9, mcs.getFirstAtomMapping().getCount());
    }

    @Test
    public void test_Pattern_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OC1=CC=CC=C1");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        AtomMatcher am;
        BondMatcher bm;
        boolean matchBonds = true;
        boolean matchAtom = true;

        if (matchAtom) {
            am = AtomMatcher.forQuery();
        } else {
            am = AtomMatcher.forElement();
        }

        if (matchBonds) {
            bm = BondMatcher.forOrder();
        } else {
            bm = BondMatcher.forAny();
        }

        boolean flagSubGraph = false;
        if (query.getAtomCount() >= target.getAtomCount()) {

            org.openscience.cdk.isomorphism.Pattern pattern = VentoFoggia.findSubstructure(target, am, bm); // create pattern
            Mappings limit = pattern.matchAll(query).limit(1);
            flagSubGraph = limit.count() > 0;
        }
        if (query.getAtomCount() < target.getAtomCount()) {
            org.openscience.cdk.isomorphism.Pattern pattern
                    = VentoFoggia.findSubstructure(query, AtomMatcher.forElement(), BondMatcher.forOrder()); // create pattern
            Mappings limit = pattern.matchAll(target).limit(1);
            flagSubGraph = limit.count() > 0;
        }
        assertEquals(true, flagSubGraph);

    }

    @Test
    public void test_Substructure_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OC1=CC=CC=C1");

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.Substructure sub = new org.openscience.smsd.Substructure(query, target, atomMatcher, bondMatcher, false);
        assertEquals(true, sub.isSubgraph());

    }

    @Test
    public void test_MCSPlusBK_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

//        IAtomContainer query = sp.parseSmiles("cccc");
//        IAtomContainer target = sp.parseSmiles("nccc");
//        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
//        IAtomContainer target = sp.parseSmiles("C1=CC=CC=C1");
        IAtomContainer query = sp.parseSmiles("*C1OC(CO)C(OP=O)C1O");
        IAtomContainer target = sp.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");
//        IAtomContainer query = sp.parseSmiles("CC1=C(CCC([O-])=O)C2=[N+]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N+]5[Fe-]3(n14)n1c(=C6)c(C)c(CCC([O-])=O)c1=C2");
//        IAtomContainer target = sp.parseSmiles("C1=2N3C(C=C4[N+]5=C(C=C6N7C8=CC9=[N+](C(=C1)C(=C9CCC([O-])=O)C)[Fe-2]573OC(CCC=%10C=%11C=C%12C(=C(C%13=CC%14=[N+]%15C(=CC=%16N%17C(=C(C%16C=C)C)C=C([N+]%11[Fe-2]%15%17(N%13%12)OC(CCC8=C6C)=O)C%10C)C(=C%14C=C)C)C)CCC([O-])=O)=O)C(=C4C)C=C)=C(C2C)C=C");
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);

        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);

//        System.out.println("MCS " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
        assertEquals("C(C1C(C(CO1)O)O)O", mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());

    }

//    @Test
//    public void test_BIG_count() throws Exception {
//        ////////System.out.println("3");
//        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
//
//        IAtomContainer query = sp.parseSmiles("C\\C1=C2\\N\\C(=C/C3=N/C(=C\\C4=C(CC(O)=O)C(CCC(O)=O)=C(CC5=C(CCC(O)=O)C(CC(O)=O)=C1N5)N4)/[C@@H](CCC(O)=O)[C@]3(C)CC(O)=O)[C@@H](CCC(O)=O)[C@]2(C)CC(O)=O");
//        IAtomContainer target = sp.parseSmiles("C[C@]1(CC(O)=O)[C@H](CCC(O)=O)\\C2=C\\C3=C(CC(O)=O)C(CCC(O)=O)=C(CC4=C(CCC(O)=O)C(CC(O)=O)=C(N4)[C@](C)(O)[C@@]45N\\C(=C/C1=N2)[C@@H](CCC(O)=O)[C@]4(C)CC(=O)O5)N3");
//        MoleculeInitializer.initializeMolecule(query);
//        MoleculeInitializer.initializeMolecule(target);
//
//        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
//        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);
//
////        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
////                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);
//        Isomorphism mcs
//                = new Isomorphism(query, target, Algorithm.VFLibMCS, atomMatcher, bondMatcher);
//        String create = new SmilesGenerator(SmiFlavor.Canonical).create(mcs.getFirstAtomMapping().getMapCommonFragmentOnQuery());
////        System.out.println("MCS " + create);
//        assertEquals(63, mcs.getFirstAtomMapping().getCount());
//
//    }
    @Test
    public void test_Ring_Match() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("OC1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OC(=C)\\C=C/C=C");
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

//        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
//                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);
        Isomorphism mcs
                = new Isomorphism(query, target, Algorithm.VFLibMCS, atomMatcher, bondMatcher);
        String create = new SmilesGenerator(SmiFlavor.Canonical).create(mcs.getFirstAtomMapping().getMapCommonFragmentOnQuery());
//        System.out.println("MCS " + create);
        assertEquals(0, mcs.getFirstAtomMapping().getCount());

    }

    @Test
    public void test_MCSPlusBig_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("O=C1C=C2CCC3C(CCC4(C)C(C(=O)CO)CCC34)C2(C)CC1");
        IAtomContainer target = sp.parseSmiles("O=CC(O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC21C");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

        try {
//            System.out.println("CALLING Substructure");
            Substructure sub
                    = new Substructure(query, target, atomMatcher, bondMatcher, true);
//            System.out.println("Subgraph " + sub.isSubgraph());
            //System.out.println("MCS " + (mcs.getFirstAtomMapping().getCommonFragmentAsSMILES()));
            assertEquals(Boolean.TRUE, sub.isSubgraph());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void MCSCheck() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("O=C(O)CCC(N)C(=O)O");
        IAtomContainer target = sp.parseSmiles("O=C(O)C(N)CCC(=O)N");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, false);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

        Isomorphism mcs
                = new Isomorphism(query, target, Algorithm.CDKMCS, atomMatcher, bondMatcher);
//        System.out.println("mcs " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
        assertEquals(9, mcs.getFirstAtomMapping().getCount());

    }

    @Test
    public void SubStructureCheck() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //.[O][C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)[C@@H2]O
        //.[O][C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)[C@@H2]O
        IAtomContainer query = sp.parseSmiles("OC[C@@H]1C[C@H](O)[C@@H](O)C(O)O1");
        IAtomContainer target1 = sp.parseSmiles("OC[C@H]1OC[C@H](O)[C@@H](O)[C@@H]1O");
        IAtomContainer target2 = sp.parseSmiles("OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target1);
        MoleculeInitializer.initializeMolecule(target2);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

        try {

            Substructure mcs
                    = new Substructure(query, target1, atomMatcher, bondMatcher, true);
//            System.out.println("sub " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
            assertEquals(0, mcs.getFirstAtomMapping().getCount());
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {

            Substructure mcs
                    = new Substructure(query, target2, atomMatcher, bondMatcher, true);
//            System.out.println("sub " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
            assertEquals(11, mcs.getFirstAtomMapping().getCount());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void MCSRingCheck() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("OC[C@H]1OC[C@H](O)[C@@H](O)[C@@H]1O");
        IAtomContainer target = sp.parseSmiles("OC[C@H]1O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O");

        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

        Isomorphism mcs
                = new Isomorphism(query, target, Algorithm.DEFAULT, atomMatcher, bondMatcher);

//            org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
//                    = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, false, false, false);
//        System.out.println("mcs " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
        assertEquals(11, mcs.getFirstAtomMapping().getCount());

    }

    @Test
    public void test_MCS_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("O[C@@H]1C=CC(=C[C@H]1OC(=C)C(O)=O)C(O)=O");
        IAtomContainer target = sp.parseSmiles("O[C@H]1C=C[C@](CC(=O)C(O)=O)(C=C1)C(O)=O");
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, true);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);
//        Isomorphism mcs
//                = new Isomorphism(query, target, Algorithm.VFLibMCS, atomMatcher, bondMatcher);
//        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
//        System.out.println("MCS " + mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());
//        System.out.println("R MCS " + smilesGenerator.create(mcs.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
//        System.out.println("P MCS " + smilesGenerator.create(mcs.getFirstAtomMapping().getMapCommonFragmentOnTarget()));

//        SmilesGenerator smileGenerator = new SmilesGenerator(
//                SmiFlavor.Unique
//                | SmiFlavor.UseAromaticSymbols
//                | SmiFlavor.AtomAtomMap
//                | SmiFlavor.Stereo
//        );
//        String smile = smileGenerator.create(mcs.getFirstAtomMapping().getCommonFragment());
//        System.out.println("SM " + smile);
        assertEquals("[cH]-1=[cH]-[c@H](-[cH2]-c-c1C(=O)O)O", mcs.getFirstAtomMapping().getCommonFragmentAsSMILES());

    }

    @Test
    public void test_RingMCS_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("O[C@@H]1COC(O)[C@H](O)[C@H]1O");
        IAtomContainer target = sp.parseSmiles("OCC1(O)OC[C@@H](O)[C@@H]1O");
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        org.openscience.smsd.algorithm.matchers.AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
        org.openscience.smsd.algorithm.matchers.BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, true);

        org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper mcs
                = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(query, target, atomMatcher, bondMatcher);
//        Isomorphism mcs
//                = new Isomorphism(query, target, Algorithm.CDKMCS, atomMatcher, bondMatcher);
        String moiety = mcs.getFirstAtomMapping().getCommonFragmentAsSMILES();
//        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.AtomAtomMap);
//        System.out.println("MCS " + moiety);
//        System.out.println("R MCS " + smilesGenerator.create(mcs.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
//        System.out.println("P MCS " + smilesGenerator.create(mcs.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
//        System.out.println("MCS " + mcs.getFirstAtomMapping().toString());
        assertEquals("C([C@@H]([C@H]([C@@H](CO)O)O)O)O", moiety);

    }
}
