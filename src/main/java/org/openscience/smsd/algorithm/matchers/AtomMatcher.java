/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.matchers;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;

/**
 * CDK class adapted SMSD
 *
 * @author John May
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public abstract class AtomMatcher {

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

    /**
     * Atoms are always compatible.
     *
     * @return a matcher for which all atoms match
     */
    public static AtomMatcher forAny() {
        return new AnyMatcher();
    }

    /**
     * Atoms are compatible if they are the same element.
     *
     * @return a matcher which checks element compatibility
     */
    public static AtomMatcher forElement() {
        return new ElementMatcher();
    }

    /**
     * Atoms are compatible if the second atom ({@code atom2}) is accepted by
     * the {@link IQueryAtom}, {@code atom1}.
     *
     * @return a matcher which checks query atom compatibility
     */
    public static AtomMatcher forQuery() {
        return new QueryMatcher();
    }

    /**
     * Atom should match ring size and is atom type or non ring atoms with same
     * atom type
     *
     */
    public static AtomMatcher forRingAtomTypeMatcher() {
        return new RingAtomTypeMatcher();
    }

    /**
     * Returns true if atom atom type matches.
     *
     * @return
     */
    public static AtomMatcher forAtomTypeMatcher() {
        return new forAtomTypeElementMatcher();
    }

    /**
     * Returns true if atom is part of a ring system (same size). It also
     * returns true if both are not part of the ring system.
     *
     * @return
     */
    public static AtomMatcher forRingMatcher() {
        return new forRingElementMatcher();
    }

    /**
     * A matcher defines all atoms as compatible.
     */
    private static final class AnyMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return true;
        }

        @Override
        public String toString() {
            return "AnyMatcher";
        }

    }

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class QueryMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return ((IQueryAtom) atom1).matches(atom2);
        }

        @Override
        public String toString() {
            return "QueryMatcher";
        }
    }

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class ElementMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return atomicNumber(atom1) == atomicNumber(atom2);
        }

        /**
         * Null safe atomic number access.
         *
         * @param atom an atom
         * @return the atomic number
         */
        private int atomicNumber(IAtom atom) {
            Integer elem = atom.getAtomicNumber();
            if (elem != null) {
                return elem;
            }
            if (atom instanceof IPseudoAtom) {
                return 0;
            }
            throw new NullPointerException("an atom had unset atomic number");
        }

        @Override
        public String toString() {
            return "ElementMatcher";
        }
    }

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class forRingElementMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return atomicNumber(atom1) == atomicNumber(atom2)
                    && isRingSizeMatch(atom1, atom2);
        }

        /**
         * Null safe atomic number access.
         *
         * @param atom an atom
         * @return the atomic number
         */
        private int atomicNumber(IAtom atom) {
            Integer elem = atom.getAtomicNumber();
            if (elem != null) {
                return elem;
            }
            if (atom instanceof IPseudoAtom) {
                return 0;
            }
            throw new NullPointerException("an atom had unset atomic number");
        }

        private boolean isRingSizeMatch(IAtom atom1, IAtom atom2) {
            if (atom1.isInRing() & atom2.isInRing()) {
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

        @Override
        public String toString() {
            return "forRingElementMatcher";
        }
    }

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class forAtomTypeElementMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return atomicNumber(atom1) == atomicNumber(atom2)
                    && matchAtomType(atom1, atom2);
        }

        /**
         * Null safe atomic number access.
         *
         * @param atom an atom
         * @return the atomic number
         */
        private int atomicNumber(IAtom atom) {
            Integer elem = atom.getAtomicNumber();
            if (elem != null) {
                return elem;
            }
            if (atom instanceof IPseudoAtom) {
                return 0;
            }
            throw new NullPointerException("an atom had unset atomic number");
        }

        private boolean matchAtomType(IAtom atom1, IAtom atom2) {
            String rAtom = atom1.getAtomTypeName() == null
                    ? atom1.getSymbol() : atom1.getAtomTypeName();
            String tAtom = atom2.getAtomTypeName() == null
                    ? atom2.getSymbol() : atom2.getAtomTypeName();
            return rAtom.equals(tAtom);
        }

        @Override
        public String toString() {
            return "forAtomTypeElementMatcher";
        }
    }

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class RingAtomTypeMatcher extends AtomMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom atom1, IAtom atom2) {
            return atomicNumber(atom1) == atomicNumber(atom2)
                    && matchAtomType(atom1, atom2)
                    && isRingSizeMatch(atom1, atom2);
        }

        /**
         * Null safe atomic number access.
         *
         * @param atom an atom
         * @return the atomic number
         */
        private int atomicNumber(IAtom atom) {
            Integer elem = atom.getAtomicNumber();
            if (elem != null) {
                return elem;
            }
            if (atom instanceof IPseudoAtom) {
                return 0;
            }
            throw new NullPointerException("an atom had unset atomic number");
        }

        private boolean matchAtomType(IAtom atom1, IAtom atom2) {
            String rAtom = atom1.getAtomTypeName() == null
                    ? atom1.getSymbol() : atom1.getAtomTypeName();
            String tAtom = atom2.getAtomTypeName() == null
                    ? atom2.getSymbol() : atom2.getAtomTypeName();
            return rAtom.equals(tAtom);
        }

        private boolean isRingSizeMatch(IAtom atom1, IAtom atom2) {
            if (atom1.isInRing() & atom2.isInRing()) {
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

        @Override
        public String toString() {
            return "RingAtomTypeMatcher";
        }

    }
}
