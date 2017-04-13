/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.algorithm.matchers;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * Checks if a bond is matching between query and target molecules.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class DefaultBondMatcher implements BondMatcher {

    static final long serialVersionUID = -7861469841127328812L;
    private final IBond queryBond;
    private final boolean shouldMatchBonds;
    private final boolean matchAtomTypes;
    private final boolean shouldMatchRings;

    /**
     * Constructor
     */
    public DefaultBondMatcher() {
        this.queryBond = null;
        shouldMatchBonds = false;
        matchAtomTypes = false;
        shouldMatchRings = false;
    }

    /**
     * Constructor
     *
     * @param queryBond query GraphMolecule
     * @param shouldMatchBonds bond match flag
     * @param shouldMatchRings
     * @param matchAtomTypes
     */
    public DefaultBondMatcher(IBond queryBond,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomTypes) {
        super();
        this.queryBond = queryBond;
        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomTypes = matchAtomTypes;
    }

    /**
     * {@inheritDoc}
     *
     * @param targetBond target bond
     * @return true if bonds match
     */
    @Override
    public boolean matches(IBond targetBond) {

        if (this.queryBond != null && queryBond instanceof IQueryBond) {
            return ((IQueryBond) queryBond).matches(targetBond);
        } else if ((queryBond != null && targetBond != null)
                && isBondMatchFlag() && isBondTypeMatch(targetBond)) {
            return true;
        } else if ((queryBond != null && targetBond != null)
                && !isBondMatchFlag() && isShouldMatchRings()) {
            if (queryBond.getFlag(CDKConstants.ISAROMATIC)
                    && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
                return true;
            } else if (!queryBond.getFlag(CDKConstants.ISAROMATIC)
                    && !targetBond.getFlag(CDKConstants.ISAROMATIC)) {
                return true;
            }
        } else if ((queryBond != null && targetBond != null)
                && !isBondMatchFlag() && !isShouldMatchRings()) {
            return true;
        }
        return false;
    }

    /**
     * Return true if a bond is matched between query and target
     *
     * @param targetBond
     * @return
     */
    private boolean isBondTypeMatch(IBond targetBond) {

        if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && queryBond.getOrder().equals(targetBond.getOrder())) {
            return true;
        }

        if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }

        return !matchAtomTypes
                && queryBond.getFlag(CDKConstants.ISINRING)
                && targetBond.getFlag(CDKConstants.ISINRING)
                && (queryBond.getOrder() == IBond.Order.UNSET
                || targetBond.getOrder() == IBond.Order.UNSET);
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isBondMatchFlag() {
        return shouldMatchBonds;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isShouldMatchRings() {
        return shouldMatchRings;
    }
}
