/* Copyright (C) 2009-2020   Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
package org.openscience.smsd;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * Checks if atom is matching between query and target molecules.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class AtomBondMatcher {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(AtomBondMatcher.class);

    /**
     *
     * @param b1
     * @param b2
     * @param atomMatcher
     * @param bondMatcher
     * @param undirected
     * @return
     */
    public static boolean matchAtomAndBond(
            IBond b1,
            IBond b2,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean undirected) {
        LOGGER.debug("matchAtomAndBond");

        boolean atomMatch = matches(b1.getBegin(), b2.getBegin(), atomMatcher)
                && matches(b1.getEnd(), b2.getEnd(), atomMatcher);
        boolean bondMatch = matches(b1, b2, bondMatcher);

        if (undirected) {
            atomMatch |= matches(b1.getBegin(), b2.getEnd(), atomMatcher)
                    && matches(b1.getEnd(), b2.getBegin(), atomMatcher);
        }

        LOGGER.debug(" bondA1 a0:" + b1.getBegin().getSymbol()
                + " a1:" + b1.getEnd().getSymbol());
        LOGGER.debug(" bondB1 b0:" + b2.getBegin().getSymbol()
                + " b1:" + b2.getEnd().getSymbol());
        LOGGER.debug(" atomMatch " + atomMatch
                + ", bondMatch " + bondMatch);
        return atomMatch && bondMatch;
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
            boolean matchAtomTypes,
            boolean shouldMatchRings) {

        AtomMatcher am = AtomMatcher.forElement();

        if (matchAtomTypes) {
            LOGGER.debug("matchAtomTypes " + matchAtomTypes);
            am = AtomMatcher.forAtomTypeMatcher();
        }

        if (shouldMatchRings && !matchAtomTypes) {
            LOGGER.debug("shouldMatchRings " + shouldMatchRings);
            am = AtomMatcher.forRingMatcher();
        }

        if (shouldMatchRings && matchAtomTypes) {
            LOGGER.debug("matchAtomTypes " + matchAtomTypes);
            am = AtomMatcher.forRingAtomTypeMatcher();
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
            LOGGER.debug("Order Match Choosen " + matchBond);
            bm = BondMatcher.forOrder();
        }

        if (shouldMatchRings) {
            LOGGER.debug("Ring Match Choosen " + shouldMatchRings);
            bm = BondMatcher.forRing();
        }

        if (matchBond && shouldMatchRings) {
            LOGGER.debug("Order & Ring Match Choosen " + shouldMatchRings);
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

    // ==================== Inner class AtomMatcher ====================

    /**
     * CDK class adapted SMSD
     *
     * @author John May
     * @author Syed Asad Rahman
     */
    public static abstract class AtomMatcher {

        @Override
        public abstract String toString();

        /**
         * Are the semantics of {@code atom1} compatible with {@code atom2}.
         *
         * @param atom1 an atom from a query container
         * @param atom2 an atom from the target container
         * @return the atom1 can be paired with atom2
         */
        public abstract boolean matches(IAtom atom1, IAtom atom2);

        // ---- Shared helper methods ----

        protected static int atomicNumber(IAtom atom) {
            Integer elem = atom.getAtomicNumber();
            if (elem != null) {
                return elem;
            }
            if (atom instanceof IPseudoAtom) {
                return 0;
            }
            throw new NullPointerException("an atom had unset atomic number");
        }

        protected static boolean matchAtomType(IAtom atom1, IAtom atom2) {
            String rAtom = atom1.getAtomTypeName() == null
                    ? atom1.getSymbol() : atom1.getAtomTypeName();
            String tAtom = atom2.getAtomTypeName() == null
                    ? atom2.getSymbol() : atom2.getAtomTypeName();
            return rAtom.equals(tAtom);
        }

        protected static boolean isRingSizeMatch(IAtom atom1, IAtom atom2) {
            if (atom1.isInRing() && atom2.isInRing()) {
                List<Integer> ringsizesQ = atom1.getProperty(CDKConstants.RING_SIZES);
                List<Integer> ringsizesT = atom2.getProperty(CDKConstants.RING_SIZES);
                if (ringsizesQ == null || ringsizesT == null) {
                    return false;
                } else {
                    return ringsizesT.containsAll(ringsizesQ)
                            || ringsizesQ.containsAll(ringsizesT);
                }
            }
            return !atom1.isAromatic() && !atom2.isAromatic();
        }

        protected static boolean matchCharge(IAtom atom1, IAtom atom2) {
            Integer c1 = atom1.getFormalCharge();
            Integer c2 = atom2.getFormalCharge();
            if (c1 == null) c1 = 0;
            if (c2 == null) c2 = 0;
            return c1.equals(c2);
        }

        protected static boolean matchIsotope(IAtom atom1, IAtom atom2) {
            Integer m1 = atom1.getMassNumber();
            Integer m2 = atom2.getMassNumber();
            if (m1 == null || m2 == null) return true;
            return m1.equals(m2);
        }

        // ---- Factory methods ----

        public static AtomMatcher forAny() {
            return new AnyMatcher();
        }

        public static AtomMatcher forElement() {
            return new ElementMatcher();
        }

        public static AtomMatcher forQuery() {
            return new QueryAtomMatcher();
        }

        public static boolean matchSymbol(IAtom atom1, IAtom atom2) {
            if (atom1.getAtomicNumber() != null && atom2.getAtomicNumber() != null) {
                return atom1.getAtomicNumber().equals(atom2.getAtomicNumber());
            }
            String s1 = atom1.getSymbol();
            String s2 = atom2.getSymbol();
            return s1 != null && s1.equals(s2);
        }

        public static AtomMatcher forRingAtomTypeMatcher() {
            return new RingAtomTypeMatcher();
        }

        public static AtomMatcher forAtomTypeMatcher() {
            return new AtomTypeElementMatcher();
        }

        public static AtomMatcher forRingMatcher() {
            return new RingElementMatcher();
        }

        // ---- Inner matcher classes ----

        private static final class AnyMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return true;
            }

            @Override
            public String toString() {
                return "AnyMatcher";
            }
        }

        private static final class QueryAtomMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return ((IQueryAtom) atom1).matches(atom2);
            }

            @Override
            public String toString() {
                return "QueryMatcher";
            }
        }

        private static final class ElementMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return atomicNumber(atom1) == atomicNumber(atom2)
                        && matchCharge(atom1, atom2)
                        && matchIsotope(atom1, atom2);
            }

            @Override
            public String toString() {
                return "ElementMatcher";
            }
        }

        private static final class RingElementMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return atomicNumber(atom1) == atomicNumber(atom2)
                        && isRingSizeMatch(atom1, atom2)
                        && matchCharge(atom1, atom2);
            }

            @Override
            public String toString() {
                return "RingElementMatcher";
            }
        }

        private static final class AtomTypeElementMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return atomicNumber(atom1) == atomicNumber(atom2)
                        && matchAtomType(atom1, atom2)
                        && matchCharge(atom1, atom2);
            }

            @Override
            public String toString() {
                return "AtomTypeElementMatcher";
            }
        }

        private static final class RingAtomTypeMatcher extends AtomMatcher {
            @Override
            public boolean matches(IAtom atom1, IAtom atom2) {
                return atomicNumber(atom1) == atomicNumber(atom2)
                        && matchAtomType(atom1, atom2)
                        && isRingSizeMatch(atom1, atom2)
                        && matchCharge(atom1, atom2);
            }

            @Override
            public String toString() {
                return "RingAtomTypeMatcher";
            }
        }
    }

    // ==================== Inner class BondMatcher ====================

    /**
     * CDK class adapted SMSD
     *
     * @author John May
     * @author Syed Asad Rahman
     */
    public static abstract class BondMatcher {

        @Override
        public abstract String toString();

        /**
         * Determines if {@code bond1} is compatible with {@code bond2}.
         *
         * @param bond1 a bond from the query structure
         * @param bond2 a bond from the target structure
         * @return the bonds are compatible
         */
        public abstract boolean matches(IBond bond1, IBond bond2);

        public static BondMatcher forAny() {
            return new AnyBondMatcher();
        }

        public static BondMatcher forStrictOrder() {
            return new StrictOrderMatcher();
        }

        public static BondMatcher forOrder() {
            return new OrderMatcher();
        }

        public static BondMatcher forRing() {
            return new RingMatcher();
        }

        public static BondMatcher forQuery() {
            return new QueryBondMatcher();
        }

        private static final class OrderMatcher extends BondMatcher {
            @Override
            public boolean matches(IBond bond1, IBond bond2) {
                return bond1.isAromatic() && bond2.isAromatic()
                        || bond1.getOrder() == bond2.getOrder();
            }

            @Override
            public String toString() {
                return "OrderMatcher";
            }
        }

        private static final class RingMatcher extends BondMatcher {
            @Override
            public boolean matches(IBond bond1, IBond bond2) {
                return (bond1.isAromatic() == bond2.isAromatic())
                        || (!bond1.isAromatic() && !bond2.isAromatic());
            }

            @Override
            public String toString() {
                return "RingMatcher";
            }
        }

        private static final class StrictOrderMatcher extends BondMatcher {
            @Override
            public boolean matches(IBond bond1, IBond bond2) {
                return bond1.isAromatic() == bond2.isAromatic()
                        && (bond1.getOrder() == bond2.getOrder()
                        || bond1.isAromatic() && bond2.isAromatic());
            }

            @Override
            public String toString() {
                return "StrictOrderMatcher";
            }
        }

        private static final class AnyBondMatcher extends BondMatcher {
            @Override
            public boolean matches(IBond bond1, IBond bond2) {
                return true;
            }

            @Override
            public String toString() {
                return "AnyMatcher";
            }
        }

        private static final class QueryBondMatcher extends BondMatcher {
            @Override
            public boolean matches(IBond bond1, IBond bond2) {
                return ((IQueryBond) bond1).matches(bond2);
            }

            @Override
            public String toString() {
                return "QueryMatcher";
            }
        }
    }
}
