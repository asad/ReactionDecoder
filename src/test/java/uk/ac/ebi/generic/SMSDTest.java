/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.generic;

import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author asad
 */
public class SMSDTest {

    /**
     * @param args the command line arguments
     * @throws java.lang.CloneNotSupportedException
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static void main(String[] args) throws CloneNotSupportedException, CDKException {
        try {
            // TODO code application logic here
            SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer ac1 = smilesParser.parseSmiles("*C1OC(CO)C(OP=O)C1O");
            IAtomContainer ac2 = smilesParser.parseSmiles("*[P](=O)(O)OCC1O[CH]C(OC(=O)C(N)CC2=CC=C(O)C=C2)C1O");

//            IAtomContainer ac1 = smilesParser.parseSmiles("OCc1[nH]c(Cc2[nH]c(Cc3[nH]c(Cc4[nH]cc(CCC(O)=O)c4CC(O)=O)c(CCC(O)=O)c3CC(O)=O)c(CCC(O)=O)c2CC(O)=O)c(CCC(O)=O)c1CC(O)=O");
//            IAtomContainer ac2 = smilesParser.parseSmiles("OC(=O)CCc1c2Cc3[nH]c(Cc4[nH]c(Cc5[nH]c(Cc([nH]2)c1CC(O)=O)c(CCC(O)=O)c5CC(O)=O)c(CCC(O)=O)c4CC(O)=O)c(CC(O)=O)c3CCC(O)=O");
            Isomorphism smsd = new Isomorphism(ac1, ac2, Algorithm.MCSPlus, false, false, true);
            smsd.setChemFilters(true, true, true);
            AtomAtomMapping aam = smsd.getAllAtomMapping().iterator().next();
            System.out.println("Mapping " + aam.getCommonFragmentAsSMILES());
            System.out.println("size " + aam.getCount());
            for (Map.Entry<IAtom, IAtom> m : aam.getMappingsByAtoms().entrySet()) {
                System.out.println(m.getKey().getID() + "" + m.getValue().getID());
            }

        } catch (InvalidSmilesException ex) {
            Logger.getLogger(SMSDTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
