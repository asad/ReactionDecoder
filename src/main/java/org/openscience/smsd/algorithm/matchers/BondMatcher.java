/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.matchers;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * @author John May
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public abstract class BondMatcher extends org.openscience.cdk.isomorphism.BondMatcher {

    /**
     * Determines if {@code bond1} is compatible with {@code bond2}.
     *
     * @param bond1 a bond from the query structure
     * @param bond2 a bond from the target structure
     * @return the bonds are compatible
     */
    @Override
    public abstract boolean matches(IBond bond1, IBond bond2);

    /**
     * All bonds are compatible.
     *
     * @return a bond matcher
     */
    public static BondMatcher forAny() {
        return new AnyMatcher();
    }

    /**
     * Bonds are compatible if they are both aromatic or their orders are equal
     * and they are non-aromatic. Under this matcher a single/double bond will
     * not match a single/double bond which is aromatic.
     *
     * @return a bond matcher
     */
    public static BondMatcher forStrictOrder() {
        return new StrictOrderMatcher();
    }

    /**
     * Bonds are compatible if they are both aromatic or their orders are equal.
     * This matcher allows a single/double bond to match a single/double
     * aromatic bond.
     *
     * @return a bond matcher
     */
    public static BondMatcher forOrder() {
        return new OrderMatcher();
    }

    /**
     * Bonds are compatible if the first {@code bond1} (an {@link IQueryBond})
     * matches the second, {@code bond2}.
     *
     * @return a bond matcher
     */
    public static BondMatcher forQuery() {
        return new QueryMatcher();
    }

    /**
     * Bonds are compatible if they are both aromatic or their orders are equal.
     */
    private static final class OrderMatcher extends BondMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IBond bond1, IBond bond2) {
            return bond1.getFlag(ISAROMATIC) && bond2.getFlag(ISAROMATIC)
                    || bond1.getOrder() == bond2.getOrder();
        }
    }

    /**
     * Bonds are compatible if they are both aromatic or their orders are equal
     * and they are non-aromatic. In this matcher a single or double bond will
     * not match a single or double bond which is part of an aromatic system.
     */
    private static final class StrictOrderMatcher extends BondMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IBond bond1, IBond bond2) {
            return bond1.getFlag(ISAROMATIC) == bond2.getFlag(ISAROMATIC)
                    && (bond1.getOrder() == bond2.getOrder()
                    || bond1.getFlag(ISAROMATIC) && bond2.getFlag(ISAROMATIC));
        }
    }

    /**
     * All bonds are considered compatible.
     */
    private static final class AnyMatcher extends BondMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IBond bond1, IBond bond2) {
            return true;
        }
    }

    /**
     * Bonds are compatible if the first {@code bond1} (an {@link IQueryBond})
     * matches the second, {@code bond2}.
     */
    private static final class QueryMatcher extends BondMatcher {

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IBond bond1, IBond bond2) {
            return ((IQueryBond) bond1).matches(bond2);
        }
    }
}
