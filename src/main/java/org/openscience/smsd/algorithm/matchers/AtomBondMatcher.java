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
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * Checks if atom is matching between query and target molecules.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class AtomBondMatcher {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(AtomBondMatcher.class);
    private final static boolean DEBUG = false;

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return
     */
    public static boolean matchAtomAndBond(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {
        if (DEBUG) {
            System.out.println("\n\n matchAtomAndBond ");
        }
        AtomMatcher atomMatcher = atomMatcher(shouldMatchRings, matchAtomType);
        boolean atomMatch1 = matches(bondA1.getBegin(), bondA2.getBegin(), atomMatcher)
                && matches(bondA1.getEnd(), bondA2.getEnd(), atomMatcher);
        boolean atomMatch2 = matches(bondA1.getBegin(), bondA2.getEnd(), atomMatcher)
                && matches(bondA1.getEnd(), bondA2.getBegin(), atomMatcher);
        BondMatcher bondMatcher = bondMatcher(shouldMatchBonds, shouldMatchRings);
        boolean bondMatch = matches(bondA1, bondA2, bondMatcher);

        if (DEBUG) {
            System.out.println(" bondA1 a0:" + bondA1.getBegin().getSymbol()
                    + " a1:" + bondA1.getEnd().getSymbol());
            System.out.println(" bondB1 b0:" + bondA2.getBegin().getSymbol()
                    + " b1:" + bondA2.getEnd().getSymbol());

            System.out.println(" atomMatch1 " + atomMatch1
                    + "| atomMatch2 " + atomMatch2
                    + ", bondMatch " + bondMatch);
        }
        return (atomMatch1 || atomMatch2) && bondMatch;
    }

    /**
     *
     * @param b1
     * @param b2
     * @param atomMatcher
     * @param bondMatcher
     * @return
     */
    public static boolean matchAtomAndBond(
            IBond b1,
            IBond b2,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher) {
        if (DEBUG) {
            System.out.println("\n\n matchAtomAndBond ");
        }
        
        boolean atomMatch1 = matches(b1.getBegin(), b2.getBegin(), atomMatcher)
                && matches(b1.getEnd(), b2.getEnd(), atomMatcher);
        boolean atomMatch2 = matches(b1.getBegin(), b2.getEnd(), atomMatcher)
                && matches(b1.getEnd(), b2.getBegin(), atomMatcher);
        boolean bondMatch = matches(b1, b2, bondMatcher);

        if (DEBUG) {
            System.out.println(" bondA1 a0:" + b1.getBegin().getSymbol()
                    + " a1:" + b1.getEnd().getSymbol());
            System.out.println(" bondB1 b0:" + b2.getBegin().getSymbol()
                    + " b1:" + b2.getEnd().getSymbol());

            System.out.println(" atomMatch1 " + atomMatch1
                    + "| atomMatch2 " + atomMatch2
                    + ", bondMatch " + bondMatch);
        }
        return (atomMatch1 || atomMatch2) && bondMatch;
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param bm Bond Matcher
     * @return
     */
    public static boolean matches(
            IBond bondA1,
            IBond bondA2,
            BondMatcher bm) {
        return bm.matches(bondA1, bondA2);
    }

    /**
     *
     * @param a1
     * @param a2
     * @param am Atom Matcher
     * @return
     */
    public static boolean matches(
            IAtom a1,
            IAtom a2,
            AtomMatcher am) {
        return am.matches(a1, a2);

    }

    /**
     * Get Atom Matcher
     *
     * @param shouldMatchRings
     * @param matchAtomTypes
     * @return
     */
    public static AtomMatcher atomMatcher(
            boolean shouldMatchRings,
            boolean matchAtomTypes) {
        AtomMatcher am = AtomMatcher.forElement();

        if (shouldMatchRings) {
            if (DEBUG) {
                System.out.println("shouldMatchRings " + shouldMatchRings);
            }
            am = AtomMatcher.forRingMatcher();
        }

        if (matchAtomTypes) {
            if (DEBUG) {
                System.out.println("matchAtomTypes " + matchAtomTypes);
            }
            am = AtomMatcher.forAtomTypeMatcher();
        }
        return am;
    }

    /**
     * Get Bond Matcher
     *
     * @param matchBond
     * @param shouldMatchRings
     * @return
     */
    public static BondMatcher bondMatcher(
            boolean matchBond,
            boolean shouldMatchRings) {
        BondMatcher bm = BondMatcher.forAny();

        if (matchBond) {
            if (DEBUG) {
                System.out.println("matchBond " + matchBond);
            }
            bm = BondMatcher.forOrder();
        }

        if (shouldMatchRings) {
            if (DEBUG) {
                System.out.println("shouldMatchRings " + shouldMatchRings);
            }
            bm = BondMatcher.forStrictOrder();
        }
        return bm;
    }

    /**
     * Query Atom Matcher
     *
     * @return
     */
    public static AtomMatcher queryAtomMatcher() {
        return AtomMatcher.forQuery();
    }

    /**
     * Query Bond Matcher
     *
     * @return
     */
    public static BondMatcher queryBondMatcher() {
        return BondMatcher.forQuery();
    }
}
