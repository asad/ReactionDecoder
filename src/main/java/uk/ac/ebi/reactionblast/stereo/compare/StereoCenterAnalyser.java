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
package uk.ac.ebi.reactionblast.stereo.compare;

import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.cip.CIPTool;
import static org.openscience.cdk.geometry.cip.CIPTool.HYDROGEN;
import static org.openscience.cdk.geometry.cip.CIPTool.checkIfAllLigandsAreDifferent;
import static org.openscience.cdk.geometry.cip.CIPTool.defineLigand;
import static org.openscience.cdk.geometry.cip.CIPTool.order;
import org.openscience.cdk.geometry.cip.ILigand;
import org.openscience.cdk.geometry.cip.VisitedAtoms;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.SaturationChecker;

/**
 * Small utility class to check an atom to see if it is a tetrahedral stereo
 * center - that is, it has 4 different neighbours. It uses the {@link CIPTool}
 * to determine this.
 *
 * @author maclean
 *
 */
public class StereoCenterAnalyser {

    /**
     * Check an atom to see if it has a potential tetrahedral stereo center.
     * This can only be true if:
     * <ol>
     * <li>It has 4 neighbours OR 3 neighbours and a single implicit
     * hydrogen</li>
     * <li>These four neighbours are different according to CIP rules</li>
     * </ol>
     * If these conditions are met, it returns true.
     *
     * @param atom the central atom of the stereocenter
     * @param atomContainer the atom container the atom is in
     * @return true if all conditions for a stereocenter are met
     */
    public static synchronized boolean hasPotentialStereoCenter(IAtom atom, IAtomContainer atomContainer) {
        List<IAtom> neighbours = atomContainer.getConnectedAtomsList(atom);
        int numberOfNeighbours = neighbours.size();
        boolean hasImplicitHydrogen = false;
        if (numberOfNeighbours == 4) {
            hasImplicitHydrogen = false;
        } else if (numberOfNeighbours == 3) {
            Integer implicitCount = atom.getImplicitHydrogenCount();
            if (implicitCount != null && implicitCount == 1) {
                hasImplicitHydrogen = true;
            } else {
                SaturationChecker checker = new SaturationChecker();
                try {
                    if (checker.calculateNumberOfImplicitHydrogens(
                            atom, atomContainer) == 1) {
                        hasImplicitHydrogen = true;
                    }
                } catch (CDKException e) {
                    e.printStackTrace();
                }
            }
            if (!hasImplicitHydrogen) {
                return false;
            }
        } else if (numberOfNeighbours > 4) {
            return false;   // not tetrahedral, anyway
        } else if (numberOfNeighbours < 3) {
            return false;   // definitely not chiral
        }
        ILigand[] ligands = new ILigand[4];
        int index = 0;
        VisitedAtoms bitSet = new VisitedAtoms();
        int chiralAtomIndex = atomContainer.indexOf(atom);
        for (IAtom neighbour : neighbours) {
            int ligandAtomIndex = atomContainer.indexOf(neighbour);
            ligands[index] = defineLigand(
                    atomContainer, bitSet, chiralAtomIndex, ligandAtomIndex);
            index++;
        }
        if (hasImplicitHydrogen) {
            ligands[index] = defineLigand(atomContainer, bitSet, chiralAtomIndex, HYDROGEN);
        }
        order(ligands);
        return checkIfAllLigandsAreDifferent(ligands);
    }

    private StereoCenterAnalyser() {
    }
}
