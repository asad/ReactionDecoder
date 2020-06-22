/* Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.tools;

import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.helper.MoleculeInitializer;

/**
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Utility {

    static final String NEW_LINE = getProperty("line.separator");

    //print matrix of type MARCS on shell
    /**
     *
     * @param MCGregor_Matrix
     * @param bondnum_A
     * @param i_bonds_A
     * @param c_bonds_A
     * @param bondnum_B
     * @param i_bonds_B
     * @param c_bonds_B
     * @return
     */
    static int print_matrix(List<Integer> MCGregor_Matrix,
            int bondnum_A, List<Integer> i_bonds_A, List<String> c_bonds_A,
            int bondnum_B, List<Integer> i_bonds_B, List<String> c_bonds_B) {

        System.out.println("bondnum_A " + bondnum_A);
        System.out.println("bondnum_B " + bondnum_B);

        System.out.println("c_bonds_A " + c_bonds_A.size());
        print_list(c_bonds_A);
        System.out.println("i_bonds_A " + i_bonds_A.size());
        print_list(i_bonds_A);
        System.out.println("c_bonds_B " + c_bonds_B.size());
        print_list(c_bonds_B);
        System.out.println("i_bonds_B " + i_bonds_B.size());
        print_list(i_bonds_B);

        System.out.print("matrix: " + NEW_LINE + "-" + "    ");
        for (int a = 0; a < bondnum_B; a++) {
            System.out.print(" " + c_bonds_B.get((a * 4) + 0) + c_bonds_B.get((a * 4) + 1));
        }
        System.out.print(NEW_LINE + "     ");
        for (int a = 0; a < bondnum_B; a++) {
            System.out.print(" " + i_bonds_B.get((a * 3) + 0) + i_bonds_B.get((a * 3) + 1));
        }
        System.out.println("");
        for (int a = 0; a < bondnum_A; a++) {
            System.out.print(c_bonds_A.get((a * 4) + 0) + "" + c_bonds_A.get((a * 4) + 1));
            System.out.print(" " + i_bonds_A.get((a * 3) + 0) + i_bonds_A.get((a * 3) + 1));
            for (int b = 0; b < bondnum_B; b++) {
                System.out.print("   " + MCGregor_Matrix.get((a * bondnum_B) + b));
            }
            System.out.println("");

        }
        System.out.println("");

        return 0;
    }

    static void print_list(List list) {
        list.stream().forEach((o) -> {
            System.out.print(o + " ");
        });
        System.out.println("");

    }

    public static List<Integer> getBubbleSort(List<Integer> unSortedVector) {
        List<Integer> sortedVector = new ArrayList<>(unSortedVector);
        int j;
        boolean flag = true;   // set flag to true to begin first pass
        int temp;   //holding variable

        while (flag) {
            flag = false;    //set flag to false awaiting a possible swap
            for (j = 1; j < sortedVector.size() - 1; j++) {
                if (sortedVector.get(j) > sortedVector.get(j + 1)) // change to > for ascending sort
                {
                    temp = sortedVector.get(j);
                    //swap elements
                    sortedVector.set(j, sortedVector.get(j + 1));
                    sortedVector.set(j + 1, temp);
                    flag = true;              //shows a swap occurred  
                }
            }
        }
//        System.out.println("Bubble Sort: " + sortedVector);
        return sortedVector;
    }

    /**
     * If either is a subgraph
     *
     * @param ac1
     * @param ac2
     * @param either
     * @return
     * @throws CDKException
     */
    public static boolean isMatch(IAtomContainer ac1, IAtomContainer ac2, boolean either) throws CDKException {

        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac1);
        MoleculeInitializer.initializeMolecule(ac1);
        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac2);
        MoleculeInitializer.initializeMolecule(ac2);

        AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, true);
        BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);

        if (ac1.getAtomCount() <= ac2.getAtomCount()) {
            Substructure pattern = new Substructure(ac1, ac2, atomMatcher, bondMatcher, false); // create pattern
            return pattern.isSubgraph();
        }
        if (either && ac1.getAtomCount() >= ac2.getAtomCount()) {
            Substructure pattern = new Substructure(ac2, ac1, atomMatcher, bondMatcher, false); // create pattern
            return pattern.isSubgraph();
        }
        return false;
    }

    /**
     * ac1 is subgraph of ac2
     *
     * @param source
     * @param target
     * @param matchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @param matchRingSize
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static Map<IAtom, IAtom> findSubgraph(
            IAtomContainer source, IAtomContainer target,
            boolean matchAtomType, boolean matchBonds, boolean shouldMatchRings,
            boolean matchRingSize) throws CDKException {

        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(source);
        MoleculeInitializer.initializeMolecule(source);

        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);
        MoleculeInitializer.initializeMolecule(target);

        AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(matchAtomType, matchRingSize);
        BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(matchBonds, shouldMatchRings);

        Substructure s;
        if (source.getAtomCount() <= target.getAtomCount()) {
            try {
                s = new Substructure(source, target, atomMatcher, bondMatcher, false);
                s.setChemFilters(true, true, true);
                return s.getFirstAtomMapping().getMappingsByAtoms();
            } catch (CDKException ex) {
                Logger.getLogger(uk.ac.ebi.reactionblast.mechanism.helper.Utility.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return new HashMap<>();
    }
}
