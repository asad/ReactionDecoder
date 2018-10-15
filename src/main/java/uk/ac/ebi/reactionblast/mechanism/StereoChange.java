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
package uk.ac.ebi.reactionblast.mechanism;

import java.io.Serializable;

import org.openscience.cdk.interfaces.IAtom;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class StereoChange implements Serializable {

    private static final long serialVersionUID = 6778787889667901L;
    private final IStereoAndConformation rAtomStereo;
    private final IStereoAndConformation pAtomStereo;
    private final IAtom rAtom;
    private final IAtom pAtom;

    /**
     *
     * @param rAtomStereo
     * @param pAtomStereo
     * @param rAtom
     * @param pAtom
     */
    public StereoChange(IStereoAndConformation rAtomStereo, IStereoAndConformation pAtomStereo, IAtom rAtom, IAtom pAtom) {
        this.rAtomStereo = rAtomStereo;
        this.pAtomStereo = pAtomStereo;
        this.rAtom = rAtom;
        this.pAtom = pAtom;
    }

    @Override
    public String toString() {
        return "StereoChange{" + "rAtomStereo=" + rAtomStereo + ", pAtomStereo=" + pAtomStereo + ", rAtom="
                + rAtom.getSymbol() + rAtom.getID() + ", pAtom=" + pAtom.getSymbol() + pAtom.getID() + '}';
    }

    /**
     * @return the rAtomStereo
     */
    public IStereoAndConformation getReactantAtomStereo() {
        return rAtomStereo;
    }

    /**
     * @return the pAtomStereo
     */
    public IStereoAndConformation getProductAtomStereo() {
        return pAtomStereo;
    }

    /**
     * @return the rAtom
     */
    public IAtom getReactantAtom() {
        return rAtom;
    }

    /**
     * @return the pAtom
     */
    public IAtom getProductAtom() {
        return pAtom;
    }
}
