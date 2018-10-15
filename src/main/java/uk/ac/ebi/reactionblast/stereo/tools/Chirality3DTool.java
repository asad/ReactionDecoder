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
package uk.ac.ebi.reactionblast.stereo.tools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.geometry.cip.CIPTool.CIP_CHIRALITY;
import static org.openscience.cdk.geometry.cip.CIPTool.CIP_CHIRALITY.NONE;
import static org.openscience.cdk.geometry.cip.CIPTool.getCIPChirality;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.stereo.StereoTool;
import static org.openscience.cdk.stereo.StereoTool.getStereo;
import org.openscience.cdk.stereo.TetrahedralChirality;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import uk.ac.ebi.reactionblast.stereo.compare.ChiralityTool;

/**
 * Takes a molecule with 3D coordinates, and uses the {@link StereoTool} and the
 * {@link CIPTool} to determine chirality of the stereo centers.
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class Chirality3DTool implements ChiralityTool {

    /**
     * Get R/S Chirality assignments for an atom container that should have 3D
     * coordinates as Point3d set in the atoms. Note that assignments of
     * IStereoAndConformation.NONE are not returned by default.
     *
     * @param atomContainer the atom container to perform assignments on
     * @return a map of chiral atoms to assignments
     */
    @Override
    public Map<IAtom, IStereoAndConformation> getTetrahedralChiralities(IAtomContainer atomContainer) {
        return getTetrahedralChiralities(atomContainer, false);
    }

    /**
     * Get R/S Chirality assignments for an atom container that should have 3D
     * coordinates as Point3d set in the atoms. If getNoneAssignments is set to
     * true, atoms with 4 neighbors that are not chiral will be mapped to
     * IStereoAndConformation.NONE.
     *
     * @param atomContainer the atom container to perform assignments on
     * @param getNoneAssigments if true, map non-chiral tetrahedral centers to
     * NONE
     * @return a map of chiral atoms to assignments
     */
    @Override
    public Map<IAtom, IStereoAndConformation> getTetrahedralChiralities(IAtomContainer atomContainer, boolean getNoneAssigments) {
        Map<IAtom, IStereoAndConformation> chiralMap = new HashMap<>();
        for (IAtom atom : atomContainer.atoms()) {
            List<IAtom> neighbours = atomContainer.getConnectedAtomsList(atom);
            if (neighbours.size() == 4) {

                IAtom n1 = neighbours.get(0);
                IAtom n2 = neighbours.get(1);
                IAtom n3 = neighbours.get(2);
                IAtom n4 = neighbours.get(3);
                Stereo stereo = getStereo(n1, n2, n3, n4);
                IAtom[] ligands = new IAtom[]{n1, n2, n3, n4};
                ITetrahedralChirality stereoCenter
                        = new TetrahedralChirality(atom, ligands, stereo);
                CIP_CHIRALITY chirality = getCIPChirality(atomContainer, stereoCenter);
                if (getNoneAssigments || chirality != CIP_CHIRALITY.NONE) {
                    switch (chirality) {
                        case NONE:
                            chiralMap.put(atom, IStereoAndConformation.NONE);
                        case R:
                            chiralMap.put(atom, IStereoAndConformation.R);
                        case S:
                            chiralMap.put(atom, IStereoAndConformation.S);
                        default:
                            chiralMap.put(atom, IStereoAndConformation.NONE);
                    }
                }
            }
        }

        return chiralMap;
    }
}
