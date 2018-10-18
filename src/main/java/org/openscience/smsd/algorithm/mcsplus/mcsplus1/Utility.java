/* Copyright (C) 2009-2018  Syed Asad Rahman <asad at ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus.mcsplus1;

import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.smsd.algorithm.matchers.DefaultMatcher;

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

    static List<Integer> getBubbleSort(List<Integer> unSortedVector) {
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
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return
     */
    protected boolean isMatchFeasible(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {

        if (!shouldMatchBonds && !shouldMatchRings && !matchAtomType) {
            return true;
        }

        if (bondA1 instanceof IQueryBond) {
            if (((IQueryBond) bondA1).matches(bondA2)) {
                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
                return atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0));
            }
            return false;
        } else {
            /*
             This one also matches atom type, not just symbols
             */
            return DefaultMatcher.matches(bondA1, bondA2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        }
    }

}
