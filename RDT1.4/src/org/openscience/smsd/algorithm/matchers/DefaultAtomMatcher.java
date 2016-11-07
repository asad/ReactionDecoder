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
 * 
 */
package org.openscience.smsd.algorithm.matchers;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;

/**
 * Checks if atom is matching between query and target molecules.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class DefaultAtomMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private final String symbol;
    private final IAtom qAtom;
    private final boolean shouldMatchRings;

    /**
     * Constructor
     *
     * @param qAtom
     * @param symbol
     * @param shouldMatchRings
     */
    public DefaultAtomMatcher(IAtom qAtom,
            String symbol,
            boolean shouldMatchRings) {
        this.qAtom = qAtom;
        this.symbol = symbol;
        this.shouldMatchRings = shouldMatchRings;
    }

    /**
     * Constructor
     *
     * @param atom query atom
     * @param shouldMatchRings ring matching flag
     */
    public DefaultAtomMatcher(IAtom atom, boolean shouldMatchRings) {
        this(atom, atom.getSymbol(), shouldMatchRings);
    }

    private boolean matchSymbol(IAtom atom) {
        if (getAtomSymbol() == null) {
            return false;
        }
        return getAtomSymbol().equals(atom.getSymbol());
    }

    /**
     * {@inheritDoc}
     *
     * @param targetAtom
     * @return true if condition meet else false
     */
    @Override
    public boolean matches(IAtom targetAtom) {
        if (targetAtom instanceof IQueryAtom) {
            return ((IQueryAtom) targetAtom).matches(getQueryAtom());
        } else if (getQueryAtom() != null && getQueryAtom() instanceof IQueryAtom) {
            return ((IQueryAtom) getQueryAtom()).matches(targetAtom);
        } else {
            if (!matchSymbol(targetAtom)) {
                return false;
            }
            if (isMatchRings() && isRingAtom(getQueryAtom()) && !isRingAtom(targetAtom)) {
                return false;
            }
            if (isMatchRings() && !isRingAtom(getQueryAtom()) && isRingAtom(targetAtom)) {
                return false;
            }
            if (isMatchRings()
                    && isRingAtom(getQueryAtom())
                    && isRingAtom(targetAtom)
                    && !isRingSizeMatch(targetAtom)) {
                return false;
            }
        }
        return true;
    }

    private boolean isRingSizeMatch(IAtom atom) {
        List<Integer> ringsizesQ = qAtom.getProperty(CDKConstants.RING_SIZES);
        List<Integer> ringsizesT = atom.getProperty(CDKConstants.RING_SIZES);
        if (ringsizesQ != null && ringsizesT != null) {
            for (int i : ringsizesQ) {
                if (ringsizesT.contains(i)) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean isAliphaticAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISALIPHATIC);
    }

    private boolean isRingAtom(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING);
    }

    /**
     * @return the qAtom
     */
    @Override
    public IAtom getQueryAtom() {
        return qAtom;
    }

    /**
     * @return the symbol
     */
    public String getAtomSymbol() {
        return symbol;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }
}
