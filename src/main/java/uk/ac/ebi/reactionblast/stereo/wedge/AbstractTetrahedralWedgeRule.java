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
package uk.ac.ebi.reactionblast.stereo.wedge;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.stereo.TetrahedralChirality;

/**
 *
 * @author Gilleain Torrance
 */
public abstract class AbstractTetrahedralWedgeRule extends WedgeRule {

    /**
     *
     * @param centralAtom
     * @param atomContainer
     * @param angleMap
     * @return
     */
    @Override
    public IStereoElement execute(IAtom centralAtom,
            IAtomContainer atomContainer,
            SortedMap<Double, IBond> angleMap) {
        int[] permutation = getMatchPermutation();
        List<IBond> bonds = new ArrayList<>(angleMap.values());
        IAtom[] ligandAtoms = new IAtom[4];
        for (int index = 0; index < 4; index++) {
            IBond bond = bonds.get(permutation[index]);
            ligandAtoms[index] = bond.getOther(centralAtom);
        }
        return new TetrahedralChirality(centralAtom, ligandAtoms, getStereo());
    }

    /**
     *
     * @return
     */
    public abstract ITetrahedralChirality.Stereo getStereo();
}
