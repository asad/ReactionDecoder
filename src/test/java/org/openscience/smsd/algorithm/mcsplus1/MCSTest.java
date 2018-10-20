/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.mcsplus1;

import java.util.List;
import static org.junit.Assert.assertEquals;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
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
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        IAtomContainer target = sp.parseSmiles("Nc1ccccc1");

        org.openscience.smsd.algorithm.mcsplus.mcsplus1.MCSPlus mcs
                = new org.openscience.smsd.algorithm.mcsplus.mcsplus1.MCSPlus(query, target, true, true, true);
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
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("*C1OC(CO)C(OP=O)C1O");
        IAtomContainer target = sp.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");

        org.openscience.smsd.algorithm.mcsplus.mcsplus1.MCSPlus mcs
                = new org.openscience.smsd.algorithm.mcsplus.mcsplus1.MCSPlus(query, target, true, false, false);
        mcs.search_cliques();
        List<List<Integer>> overlaps = mcs.getFinalMappings();
        assertEquals(overlaps.size(), 1);
        assertEquals(9, (mcs.getFinalMappings().iterator().next().size() / 2));
    }

    /**
     * Test of set method, of class MCSPlusHandler.
     *
     * @throws Exception
     */
    @Test
    public void test_MCS_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //CHEBI:58504
        IAtomContainer query = sp.parseSmiles("CC1=C2N3[C@H]([C@H](CC(N)=O)[C@@]2"
                + "(C)CCC([O-])=O)[C@]2(C)[N+]4=C([C@@H](CCC(N)=O)[C@]2(C)CC(N)=O)"
                + "C(C)=C2[N+]5=C(C=C6[N+](=C1[C@@H](CCC(N)=O)C6(C)C)[Co--]345"
                + "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc3c(N)ncnc13)[C@@H]"
                + "(CCC(N)=O)[C@]2(C)CC(N)=O");
        //CHEBI:2480
        IAtomContainer target = sp.parseSmiles("[H][C@@]12[C@H](CC(N)=O)[C@@](C)"
                + "(CCC(=O)NC[C@@H](C)O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]"
                + "([C@]1(C)[C@@](C)(CC(N)=O)[C@@H]8CCC(N)=O)[Co--]57(C[C@H]1O[C@H]"
                + "([C@H](O)[C@@H]1O)n1cnc5c(N)ncnc15)N23)[C@@](C)(CC(N)=O)"
                + "[C@@H]6CCC(N)=O)C(C)(C)[C@@H]4CCC(N)=O");

        Isomorphism mcs
                //= new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
                = new Isomorphism(query, target, Algorithm.DEFAULT, false, false, false);

        List<AtomAtomMapping> allAtomMapping = mcs.getAllAtomMapping();
        assertEquals(81, mcs.getFirstAtomMapping().getCount());
    }

    @Test
    public void test_VFMCS_count() throws Exception {
        ////////System.out.println("3");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = sp.parseSmiles("*C1OC(CO)C(OP=O)C1O");
//        IAtomContainer target = sp.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");

        IAtomContainer query = sp.parseSmiles("[N].[N]=C[N].[CH].[O].[NH2]");
        IAtomContainer target = sp.parseSmiles("[C]C(=O)N.[N].[N]C=[N].[NH2]");

        Isomorphism mcs
                //= new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
                = new Isomorphism(query, target, Algorithm.VFLibMCS, false, false, false);

        System.out.println("mcs.getFirstAtomMapping() " + mcs.getFirstAtomMapping());
        List<AtomAtomMapping> allAtomMapping = mcs.getAllAtomMapping();
        // assertEquals(2, allAtomMapping.size());
        assertEquals(7, mcs.getFirstAtomMapping().getCount());
    }

}
