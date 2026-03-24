/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd;

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

    // ---- Shared helper methods (extracted from inner classes) ----

    /**
     * Null safe atomic number access.
     *
     * @param atom an atom
     * @return the atomic number, or 0 for pseudo-atoms
     * @throws NullPointerException if the atom has no atomic number and is not a pseudo-atom
     */
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

    /**
     * Check whether two atoms have the same atom type name (falling back to symbol).
     *
     * @param atom1 first atom
     * @param atom2 second atom
     * @return true if atom type names (or symbols) are equal
     */
    protected static boolean matchAtomType(IAtom atom1, IAtom atom2) {
        String rAtom = atom1.getAtomTypeName() == null
                ? atom1.getSymbol() : atom1.getAtomTypeName();
        String tAtom = atom2.getAtomTypeName() == null
                ? atom2.getSymbol() : atom2.getAtomTypeName();
        return rAtom.equals(tAtom);
    }

    /**
     * Check whether two atoms belong to compatible ring systems.
     * Returns true if both are in rings that share at least one common ring size,
     * or if neither atom is aromatic (non-ring atoms).
     *
     * @param atom1 first atom
     * @param atom2 second atom
     * @return true if ring sizes are compatible
     */
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

    /**
     * Check whether two atoms have the same formal charge.
     * Null charges are treated as zero.
     *
     * @param atom1 first atom
     * @param atom2 second atom
     * @return true if formal charges are equal
     */
    protected static boolean matchCharge(IAtom atom1, IAtom atom2) {
        Integer c1 = atom1.getFormalCharge();
        Integer c2 = atom2.getFormalCharge();
        if (c1 == null) c1 = 0;
        if (c2 == null) c2 = 0;
        return c1.equals(c2);
    }

    /**
     * Check whether two atoms have compatible isotope (mass number).
     * Only requires a match if BOTH atoms have an explicit mass number set;
     * if either is null the atoms are considered compatible.
     *
     * @param atom1 first atom
     * @param atom2 second atom
     * @return true if isotopes are compatible
     */
    protected static boolean matchIsotope(IAtom atom1, IAtom atom2) {
        Integer m1 = atom1.getMassNumber();
        Integer m2 = atom2.getMassNumber();
        if (m1 == null || m2 == null) return true;
        return m1.equals(m2);
    }

    // ---- Factory methods ----

    /**
     * Atoms are always compatible.
     *
     * @return a matcher for which all atoms match
     */
    public static AtomMatcher forAny() {
        return new AnyMatcher();
    }

    /**
     * Atoms are compatible if they are the same element, have the same charge,
     * and compatible isotopes.
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
     * Simple symbol-based match for single atom mapping.
     *
     * @param atom1 first atom
     * @param atom2 second atom
     * @return true if atoms have the same element symbol
     */
    public static boolean matchSymbol(IAtom atom1, IAtom atom2) {
        if (atom1.getAtomicNumber() != null && atom2.getAtomicNumber() != null) {
            return atom1.getAtomicNumber().equals(atom2.getAtomicNumber());
        }
        String s1 = atom1.getSymbol();
        String s2 = atom2.getSymbol();
        return s1 != null && s1.equals(s2);
    }

    /**
     * Atom should match ring size, atom type, and charge; non-ring atoms
     * must share the same atom type and charge.
     *
     * @return a matcher for ring-aware atom type matching
     */
    public static AtomMatcher forRingAtomTypeMatcher() {
        return new RingAtomTypeMatcher();
    }

    /**
     * Returns a matcher that checks element, atom type, and charge.
     *
     * @return a matcher for atom type matching
     */
    public static AtomMatcher forAtomTypeMatcher() {
        return new AtomTypeElementMatcher();
    }

    /**
     * Returns a matcher that checks element, ring size, and charge.
     * Returns true if both atoms are in compatible ring systems or both
     * are non-aromatic.
     *
     * @return a matcher for ring-aware element matching
     */
    public static AtomMatcher forRingMatcher() {
        return new RingElementMatcher();
    }

    // ---- Inner matcher classes ----

    /**
     * A matcher defines all atoms as compatible.
     */
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

    /**
     * A matcher to use when all atoms are {@link IQueryAtom}s. {@code atom1} is
     * cast to a query atom and matched against {@code atom2}.
     */
    private static final class QueryMatcher extends AtomMatcher {

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
     * Matches atoms by element (atomic number), formal charge, and isotope.
     */
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

    /**
     * Matches atoms by element, ring size compatibility, and charge.
     */
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

    /**
     * Matches atoms by element, atom type name, and charge.
     */
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

    /**
     * Matches atoms by element, atom type name, ring size compatibility, and charge.
     */
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
