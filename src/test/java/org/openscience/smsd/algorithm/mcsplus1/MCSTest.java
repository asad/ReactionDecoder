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

        org.openscience.smsd.algorithm.mcsplus1.MCSPlus mcs
                = new org.openscience.smsd.algorithm.mcsplus1.MCSPlus(query, target, true, true, true);
        mcs.search_cliques();
        List<List<Integer>> overlaps = mcs.getFinalMappings();
        assertEquals(overlaps.size(), 2);
        assertEquals((mcs.getFinalMappings().iterator().next().size() / 2), 7);
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

        org.openscience.smsd.algorithm.mcsplus1.MCSPlus mcs
                    = new org.openscience.smsd.algorithm.mcsplus1.MCSPlus(query, target, true, false, false);
            mcs.search_cliques();
        List<List<Integer>> overlaps = mcs.getFinalMappings();
        assertEquals(overlaps.size(), 1);
        assertEquals((mcs.getFinalMappings().iterator().next().size() / 2), 9);
    }

}
