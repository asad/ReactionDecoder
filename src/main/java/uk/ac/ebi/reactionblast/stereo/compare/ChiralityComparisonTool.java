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

import static java.lang.System.out;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.Isomorphism;
import static org.openscience.smsd.interfaces.Algorithm.DEFAULT;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import uk.ac.ebi.reactionblast.stereo.tools.Chirality2DTool;
import uk.ac.ebi.reactionblast.stereo.tools.Chirality3DTool;

/**
 * Tool for comparing chiralities.
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class ChiralityComparisonTool {

    /**
     *
     * @param atomContainerA
     * @param atomContainerB
     * @throws Exception
     */
    public static void compare(IAtomContainer atomContainerA, IAtomContainer atomContainerB) throws Exception {
        Map<IAtom, IStereoAndConformation> chiralityMapA;
        if (has3DCoordinates(atomContainerA)) {
            chiralityMapA = new Chirality3DTool().getTetrahedralChiralities(atomContainerA);
        } else {
            chiralityMapA = new Chirality2DTool().getTetrahedralChiralities(atomContainerA);
        }

        Map<IAtom, IStereoAndConformation> chiralityMapB;
        if (has3DCoordinates(atomContainerB)) {
            chiralityMapB = new Chirality3DTool().getTetrahedralChiralities(atomContainerB);
        } else {
            chiralityMapB = new Chirality2DTool().getTetrahedralChiralities(atomContainerB);
        }

        Isomorphism isomorphism = new Isomorphism(atomContainerA, atomContainerB, DEFAULT, true, false, false);
        Map<IAtom, IAtom> atomMap = isomorphism.getFirstAtomMapping().getMappingsByAtoms();
        atomMap.keySet().stream().forEach((atomA) -> {
            IAtom atomB = atomMap.get(atomA);
            boolean isStereoA = chiralityMapA.containsKey(atomA);
            boolean isStereoB = chiralityMapB.containsKey(atomB);
            if (isStereoA && isStereoB) {
                IStereoAndConformation stereoA = chiralityMapA.get(atomA);
                IStereoAndConformation stereoB = chiralityMapB.get(atomB);
                out.println(stereoA + " " + stereoB);
            }
        });
    }

    private static boolean has3DCoordinates(IAtomContainer atomContainerA) {
        // XXX - check all atoms?
        for (IAtom atom : atomContainerA.atoms()) {
            return atom.getPoint3d() != null;
        }
        return false;
    }

    private ChiralityComparisonTool() {
    }
}
