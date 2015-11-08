/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

package uk.ac.ebi.reactionblast.stereo.tools;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.geometry.cip.CIPTool.CIP_CHIRALITY;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import uk.ac.ebi.reactionblast.stereo.compare.ChiralityTool;
import uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoLifter;

/**
 * Tool for comparing chiralities.
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class Chirality2DTool implements ChiralityTool {

    @Override
    public Map<IAtom, IStereoAndConformation> getTetrahedralChiralities(IAtomContainer atomContainer) {
        return getTetrahedralChiralities(atomContainer, false);
    }

    /**
     * 
     * @param atomContainer
     * @param getNoneAssesments
     * @return
     */
    @Override
    public Map<IAtom, IStereoAndConformation> getTetrahedralChiralities(
            IAtomContainer atomContainer, boolean getNoneAssesments) {
        Map<IAtom, IStereoAndConformation> chiralities = new HashMap<>();
        WedgeStereoLifter lifter = new WedgeStereoLifter();
        for (IAtom atom : atomContainer.atoms()) {
            IStereoAndConformation chirality = getChirality2D(lifter, atom, atomContainer);
            if (getNoneAssesments || chirality != IStereoAndConformation.NONE) {
                chiralities.put(atom, chirality);
            }
        }
        return chiralities;
    }

    /**
     * 
     * @param lifter
     * @param atom
     * @param atomContainer
     * @return
     */
    public static IStereoAndConformation getChirality2D(
            WedgeStereoLifter lifter, IAtom atom, IAtomContainer atomContainer) {
        IStereoElement stereoElement = lifter.lift(atom, atomContainer);
        return getChirality2D(stereoElement, atomContainer);
    }

    /**
     * 
     * @param stereoElement
     * @param atomContainer
     * @return
     */
    public static IStereoAndConformation getChirality2D(
            IStereoElement stereoElement, IAtomContainer atomContainer) {
        if (stereoElement instanceof ITetrahedralChirality) {
            CIP_CHIRALITY chiral = CIPTool.getCIPChirality(
                    atomContainer, (ITetrahedralChirality) stereoElement);
            switch (chiral) {
                case NONE: return IStereoAndConformation.NONE;
                case R: return IStereoAndConformation.R;
                case S: return IStereoAndConformation.S;
                default: return IStereoAndConformation.NONE;
            }
        } else {
            return IStereoAndConformation.NONE;
        }
    }
    private static final Logger LOG = Logger.getLogger(Chirality2DTool.class.getName());
}
