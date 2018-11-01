/* Copyright (C) 2009-2018   Syed Asad Rahman <asad at ebi.ac.uk>
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
 *
 */
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 * Checks if atom is matching between query and target molecules.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class AtomBondMatcher {

//    /**
//     *
//     * @param queryBondMatcher
//     * @param targetBond
//     * @return true if condition meet else false
//     */
//    public static boolean isBondMatch(BondMatcher queryBondMatcher, IBond targetBond) {
//        return queryBondMatcher.matches(targetBond);
//    }
//
//    /**
//     *
//     * @param bondA1
//     * @param bondA2
//     * @param shouldMatchRings
//     * @param matchAtomTypes
//     * @return
//     */
//    public static boolean isAtomMatch(IBond bondA1, IBond bondA2, boolean shouldMatchRings, boolean matchAtomTypes) {
//
//        AtomMatcher atomMatcher1;
//        AtomMatcher atomMatcher2;
//        atomMatcher1 = atomMatcher(bondA1.getAtom(0), shouldMatchRings, matchAtomTypes);
//        atomMatcher2 = atomMatcher(bondA1.getAtom(1), shouldMatchRings, matchAtomTypes);
//        return ((isAtomMatch(atomMatcher1, bondA2.getAtom(0)) && isAtomMatch(atomMatcher2, bondA2.getAtom(1)))
//                || (atomMatcher1.matches(bondA2.getAtom(1)) && atomMatcher2.matches(bondA2.getAtom(0))));
//
//    }
//
//    /**
//     *
//     * @param a
//     * @param shouldMatchRings
//     * @param matchAtomTypes
//     * @return
//     */
//    public static AtomMatcher atomMatcher(IAtom a, boolean shouldMatchRings, boolean matchAtomTypes) {
//        AtomMatcher atomMatcher;
//
//        if (matchAtomTypes) {
//            atomMatcher = new DefaultAtomTypeMatcher(a, shouldMatchRings);
//        } else {
//            atomMatcher = new DefaultAtomMatcher(a, shouldMatchRings);
//        }
//        return atomMatcher;
//    }
//
//    /**
//     *
//     * @param a
//     * @param b
//     * @return
//     */
//    public static boolean isAtomMatch(AtomMatcher a, IAtom b) {
//        return a.matches(b);
//    }
//
//    /**
//     *
//     * @param bondA1
//     * @param bondA2
//     * @param matchBond
//     * @param shouldMatchRings
//     * @param matchAtomTypes (atom type also matched and symbol matched)
//     * @return true if condition meet else false
//     */
//    public static boolean matches(
//            IBond bondA1,
//            IBond bondA2,
//            boolean matchBond,
//            boolean shouldMatchRings,
//            boolean matchAtomTypes) {
//
//        if (bondA1 instanceof IQueryBond) {
//            if (((IQueryBond) bondA1).matches(bondA2)) {
//                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
//                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
//                return atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
//                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0));
//            }
//            return false;
//        }
//
//        if (!isAtomMatch(bondA1, bondA2, shouldMatchRings, matchAtomTypes)) {
//            return false;
//        }
//        return isBondMatch(new DefaultBondMatcher(bondA1, matchBond, shouldMatchRings, matchAtomTypes), bondA2);
//    }
    /**
     *
     * @param bondA1
     * @param bondA2
     * @param matchBond
     * @param shouldMatchRings
     * @return
     */
    public static boolean matches(
            IBond bondA1,
            IBond bondA2,
            boolean matchBond,
            boolean shouldMatchRings) {

        BondMatcher bm = BondMatcher.forAny();

        if (matchBond) {
            bm = BondMatcher.forOrder();
        }

        if (shouldMatchRings) {
            bm = BondMatcher.forStrictOrder();
        }

        return bm.matches(bondA1, bondA2);
    }

    /**
     *
     * @param a1
     * @param a2
     * @param shouldMatchRings
     * @param matchAtomTypes
     * @return
     */
    public static boolean matches(
            IAtom a1,
            IAtom a2,
            boolean shouldMatchRings,
            boolean matchAtomTypes) {

        AtomMatcher am = AtomMatcher.forElement();

        if (shouldMatchRings) {
            am = AtomMatcher.forRingMatcher();
        }

        if (matchAtomTypes) {
            am = AtomMatcher.forAtomTypeMatcher();
        }

        return am.matches(a1, a1);
    }
}
