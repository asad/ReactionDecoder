/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.mcsplus1;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 *
 * @author arahman
 */
public class MCSTest {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException, Exception {
//        String querySM = "OCc1[nH]c(Cc2[nH]c(Cc3[nH]c(Cc4[nH]cc(CCC(O)=O)c4CC(O)=O)c(CCC(O)=O)c3CC(O)=O)c(CCC(O)=O)c2CC(O)=O)c(CCC(O)=O)c1CC(O)=O";
//        String targetSM = "OC(=O)CCc1c2Cc3[nH]c(Cc4[nH]c(Cc5[nH]c(Cc([nH]2)c1CC(O)=O)c(CCC(O)=O)c5CC(O)=O)c(CCC(O)=O)c4CC(O)=O)c(CC(O)=O)c3CCC(O)=O";
        String querySM = "*C1OC(CO)C(OP=O)C1O";
        String targetSM = "*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        try {
            IAtomContainer file1 = smilesParser.parseSmiles(querySM);
            IAtomContainer file2 = smilesParser.parseSmiles(targetSM);
            Isomorphism smsd = new Isomorphism(file1, file2, Algorithm.MCSPlus, true, false, false);
            System.out.println("smsd " + smsd.getAllAtomMapping().size());
            System.out.println("smsd first " + smsd.getAllAtomMapping().iterator().next().getCount());
        } catch (Exception ex) {
            Logger.getLogger(MCSTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("=====================================================");
        try {
            IAtomContainer file1 = smilesParser.parseSmiles(querySM);
            IAtomContainer file2 = smilesParser.parseSmiles(targetSM);

            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(file1);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(file2);
            ExtAtomContainerManipulator.aromatizeCDK(file1);
            ExtAtomContainerManipulator.aromatizeCDK(file2);

            org.openscience.smsd.algorithm.mcsplus1.MCSPlus mcs
                    = new org.openscience.smsd.algorithm.mcsplus1.MCSPlus(file1, file2, true, false, false);
            mcs.search_cliques();
            System.out.println("mcs.final_MAPPINGS " + mcs.getFinalMappings().iterator().next().size() / 2);
            List<List<Integer>> overlaps = mcs.getFinalMappings();
            System.out.println("Total Mappings " + overlaps.size());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(MCSTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    static IAtomContainer getMol(String MolFile) throws Exception {
        FileInputStream ReadMolecule;
        ReadMolecule = new FileInputStream(MolFile);
        MDLV2000Reader MolRead = new MDLV2000Reader(new InputStreamReader(ReadMolecule));
        return (IAtomContainer) MolRead.read(new AtomContainer());
    }

}
