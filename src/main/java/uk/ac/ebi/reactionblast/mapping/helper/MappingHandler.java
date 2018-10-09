/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.mapping.helper;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;

/**
 * @Author: Syed Asad Rahman <asad @ ebi.ac.uk>
 * @Date: 2009/06/3
 * @Revision: 1.10
 */
public class MappingHandler extends BasicDebugger {

    private static final Logger LOG = getLogger(MappingHandler.class.getName());

    /**
     *
     * @param MappedReaction
     */
    public static void cleanMapping(IReaction MappedReaction) {
        int count = MappedReaction.getMappingCount();
        for (int i = count; i > 0; i--) {
            MappedReaction.removeMapping(i);
        }

        for (int eMol = 0; eMol < MappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = MappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {

//                IAtom atomE = eMolecule.getAtom(eAtom);
//                System.out.println("Atom: " + atomE.getSymbol());
                IAtom atomEMap = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
//                System.out.println("AtomCopy: " + atomEMap.getSymbol());
                String atomLabel = Integer.toString(-1);
                atomEMap.setFlag(MAPPED, false);
                atomEMap.setID(atomLabel);
            }
        }

        for (int pMol = 0; pMol < MappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = MappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                IAtom atomPMap = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                String atomLabel = Integer.toString(-1);
                atomPMap.setFlag(MAPPED, false);
                atomPMap.setID(atomLabel);
            }
        }
    }

    /**
     *
     * @param expLabReaction
     * @param MappedReaction
     * @param counter
     * @return
     */
    protected static synchronized int setMappingFlags(IReaction expLabReaction, IReaction MappedReaction, int counter) {
        IAtomContainerSet expEductSet = expLabReaction.getReactants();
        IAtomContainerSet expProductSet = expLabReaction.getProducts();

//        System.out.println("Mapping size " + expLabReaction.getMappingCount());
        for (IMapping map : expLabReaction.mappings()) {

            IAtom I_Atom = (IAtom) map.getChemObject(0);
            IAtom J_Atom = (IAtom) map.getChemObject(1);

//            System.out.println("Mapped Atom ID " + I_Atom.getID().trim() + "  Atom ID " + J_Atom.getID().trim());
            if (I_Atom != null && J_Atom != null) {

                /*
                *******************************
                * Mapping the Reactants ******************************
                 */
                boolean eFlag = false;
                IAtom firstAtom = null;
                IAtom secondAtom = null;
                for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                    IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                    for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                        if (I_Atom.getID().trim().equalsIgnoreCase(eMolecule.getAtom(eAtom).getID().trim())) {

                            String atomLabel = Integer.toString(counter);
//                            System.out.println("_atomMappings.get(i).getSymbol().trim() " + I_Atom.getSymbol().trim() + " eMolecule.getAtom(eAtom).getID().trim() " + eMolecule.getAtom(eAtom).getSymbol().trim());

                            firstAtom = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                            firstAtom.setID(atomLabel);
                            firstAtom.setFlag(MAPPED, true);
                            eFlag = true;
                            break;
                        }
                    }

                    if (eFlag) {
                        break;
                    }

                }
                /*
                *******************************
                * Mapping the Products ******************************
                 */
                boolean pFlag = false;
                for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                    IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                    for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                        if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {
//                            System.out.println("Hi Matched product");
//                            System.out.println("ID:" + J_Atom.getID().trim());

                            String atomLabel = Integer.toString(counter);

                            secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                            secondAtom.setID(atomLabel);
                            secondAtom.setFlag(MAPPED, true);
                            IMapping mappingObject = MappedReaction.getBuilder().newInstance(IMapping.class, firstAtom, secondAtom);
                            MappedReaction.addMapping(mappingObject);
                            counter++;
                            pFlag = true;
                            break;
                        }
                    }

                    if (pFlag) {
                        break;
                    }
                }
            }
        }
        return counter;
    }

    /**
     *
     * @param MappedReaction
     * @param ReactionWithUniqueSTOICHIOMETRY
     * @param coreMappedReaction
     * @param counter
     * @return
     */
    protected static int setMappingFlags(IReaction MappedReaction, IReaction ReactionWithUniqueSTOICHIOMETRY, IReaction coreMappedReaction, int counter) {
        IAtomContainerSet expEductSet = ReactionWithUniqueSTOICHIOMETRY.getReactants();
        IAtomContainerSet expProductSet = ReactionWithUniqueSTOICHIOMETRY.getProducts();

//        System.out.println("Mapping size " + coreMappedReaction.getMappingCount());
        for (IMapping map : coreMappedReaction.mappings()) {

            IAtom I_Atom = (IAtom) map.getChemObject(0);
            IAtom J_Atom = (IAtom) map.getChemObject(1);
//
//            System.out.println("Mapped Atom ID " + I_Atom.getID().trim() + "  Atom ID " + J_Atom.getID().trim());
//            System.out.println("Mapped Atom Symbol " + I_Atom.getSymbol() + "  Atom Symbol " + J_Atom.getSymbol());

            if (I_Atom != null && J_Atom != null) {

                /*
                *******************************
                * Mapping the Reactants ******************************
                 */
                boolean eFlag = false;
                IAtom firstAtom = null;
                IAtom secondAtom = null;
                for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                    IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                    for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                        if (I_Atom.getID().trim().equalsIgnoreCase(eMolecule.getAtom(eAtom).getID().trim())) {

//                            System.out.println("Mapped Atom ID " + I_Atom.getID().trim() + " eMoleculeID " + eMolecule.getAtom(eAtom).getID().trim());
//                            System.out.println("Mapped Atom Symbol " + I_Atom.getSymbol() + " eMolecule Symbol " + eMolecule.getAtom(eAtom).getSymbol());
//                            System.out.println("Hi Matched product");
//                            System.out.println("ID:" + I_Atom.getID().trim());
                            String atomLabel = Integer.toString(counter);
//                            System.out.println("_atomMappings.get(i).getSymbol().trim() " + I_Atom.getSymbol().trim() + " eMolecule.getAtom(eAtom).getID().trim() " + eMolecule.getAtom(eAtom).getSymbol().trim());

                            firstAtom = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                            firstAtom.setID(atomLabel);
                            firstAtom.setFlag(MAPPED, true);
                            eFlag = true;
                            break;
                        }
                    }

                    if (eFlag) {
                        break;
                    }

                }
                /*
                *******************************
                * Mapping the Products ******************************
                 */
                boolean pFlag = false;
                for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                    IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                    for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                        if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {
//                            System.out.println("Hi Matched product");
//                            System.out.println("ID:" + J_Atom.getID().trim());

                            String atomLabel = Integer.toString(counter);
//                            System.out.println("_atomMappings.get(j).getSymbol().trim() " + J_Atom.getSymbol().trim()
//                                    + " eMolecule.getAtom(pAtom).getID().trim() " + pMolecule.getAtom(pAtom).getSymbol().trim());

                            secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                            secondAtom.setID(atomLabel);
                            secondAtom.setFlag(MAPPED, true);
                            counter++;
                            pFlag = true;
                            break;
                        }
                    }

                    if (pFlag) {
                        break;
                    }

                }
            }

        }

//        System.out.println("expLabReaction in setMapping: ");
//        printReaction(expLabReaction);
//        System.out.println("MappedReaction in setMapping: ");
//        printReaction(MappedReaction);
//        System.out.println("Mapping Counter " + counter);
        return counter;
    }
}
