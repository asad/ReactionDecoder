/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.tools;

import java.io.IOException;
import java.util.Comparator;

import static org.openscience.cdk.config.Isotopes.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getSingleBondEquivalentSum;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * GraphAtomContainer Comparator
 */
public class AtomContainerSetComparator implements Comparator<IAtomContainer> {

    /**
     * Configure LoggingTool
     */
    private final ILoggingTool LOGGER
            = createLoggingTool(AtomContainerSetComparator.class);

    /**
     * Creates a new instance of AtomContainerComparator
     */
    public AtomContainerSetComparator() {
    }

    /*
     * <p>Compares two IAtomContainers for order with the following criteria with decreasing priority:</p>
     * <ul>
     *   <li>Compare atom count
     *   <li>Compare molecular weight (heavy atoms only)
     *   <li>Compare bond count
     *   <li>Compare sum of bond orders (heavy atoms only)
     * </ul>
     * <p>If no difference can be found with the above criteria, the IAtomContainers are
     * considered equal.</p>
     * <p>Returns a negative integer, zero, or a positive integer as the first argument is less than,
     * equal to, or greater than the second.</p>
     * <p>This method is null safe.</p>
     *
     * @param o1 the first IAtomContainer
     * @param o2 the second IAtomContainer
     * @return a negative integer, zero, or a positive integer as the first argument is less than, equal
     *         to, or greater than the second.
     */
    /**
     *
     * @param o1
     * @param o2
     * @return
     */
    @Override
    public int compare(IAtomContainer o1, IAtomContainer o2) {
        // Check for nulls
        if (o1 == null && o2 == null) {
            return 0;
        }
        if (o1 == null) {
            return -1;
        }
        if (o2 == null) {
            return 1;
        }

        // Check for correct instances
        if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
            return 0;
        }
        if (!(o1 instanceof IAtomContainer)) {
            return -1;
        }
        if (!(o2 instanceof IAtomContainer)) {
            return 1;
        }

        // Check for correct instances
        if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
            return 0;
        }
        if (!(o1 instanceof IAtomContainer)) {
            return -1;
        }
        if (!(o2 instanceof IAtomContainer)) {
            return 1;
        }

        IAtomContainer atomContainer1 = o1;
        IAtomContainer atomContainer2 = o2;

        // 1. Compare atom count
        if (atomContainer1.getAtomCount() > atomContainer2.getAtomCount()) {
            return -1;
        } else if (atomContainer1.getAtomCount() < atomContainer2.getAtomCount()) {
            return 1;
        } else {
            // 2. Atom count equal, compare molecular weight (heavy atoms only)
            double mw1;
            double mw2;
            try {
                mw1 = getMolecularWeight(atomContainer1);
                mw2 = getMolecularWeight(atomContainer2);
            } catch (CDKException e) {
                LOGGER.warn("Exception in molecular mass calculation.");
                return 0;
            }
            if (mw1 > mw2) {
                return -1;
            } else if (mw1 < mw2) {
                return 1;
            } else {
                // 3. Molecular weight equal, compare bond count
                if (atomContainer1.getBondCount() > atomContainer2.getBondCount()) {
                    return -1;
                } else if (atomContainer1.getBondCount() < atomContainer2.getBondCount()) {
                    return 1;
                } else {
                    // 4. Bond count equal, compare sum of bond orders (heavy atoms only)
                    double bondOrderSum1 = getSingleBondEquivalentSum(atomContainer1);
                    double bondOrderSum2 = getSingleBondEquivalentSum(atomContainer2);
                    if (bondOrderSum1 > bondOrderSum2) {
                        return -1;
                    } else if (bondOrderSum1 < bondOrderSum2) {
                        return 1;
                    }
                }

            }
        }
        // AtomContainers are equal in terms of this comparator
        return 0;
    }

    /**
     * Returns the molecular weight (exact mass) of the major isotopes of all
     * heavy atoms of the given IAtomContainer.
     *
     * @param atomContainer an IAtomContainer to calculate the molecular weight
     * for
     * @throws org.openscience.cdk.exception.CDKException if an error occurs
     * with the IsotopeFactory
     * @return the molecular weight (exact mass) of the major isotopes of all
     * heavy atoms of the given IAtomContainer
     */
    private double getMolecularWeight(IAtomContainer atomContainer) throws CDKException {
        double mw = 0.0;
        try {
            for (IAtom atom : atomContainer.atoms()) {
                if (!atom.getSymbol().equals("H") && !atom.getSymbol().equals("R")) {
                    try {
                        try {
                            IIsotope majorIsotope = getInstance().getMajorIsotope(atom.getSymbol());
                            mw += majorIsotope.getExactMass();
                        } catch (NullPointerException e) {
                            mw += getInstance().getMajorIsotope("Ra").getExactMass();
                            LOGGER.warn("Isotopes not defined in the CDK " + atom.getSymbol());
                        }
                    } catch (IOException e) {
                        LOGGER.warn("Molecular weight calculation failed for atom " + atom.getSymbol());
                    }
                } else if (atom.getSymbol().equals("R")) {
                    mw += getInstance().getMajorIsotope("C").getExactMass();
                }
            }
        } catch (IOException e) {
            LOGGER.warn("Molecular weight calculation failed for atleast one atom ");
        }
        return mw;
    }
}
