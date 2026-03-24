/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.Serializable;

import org.openscience.cdk.interfaces.IAtom;
import com.bioinceptionlabs.reactionblast.mechanism.IStereoAndConformation;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class AtomStereoChangeInformation implements Serializable {

    private static final long serialVersionUID = 1896986585959789L;
    private final IAtom reactantAtom;
    private final IAtom productAtom;
    private boolean stereoChange = false;
    private IStereoAndConformation atomStereoR;
    private IStereoAndConformation atomStereoP;

    /**
     *
     * @param rAtom
     * @param pAtom
     */
    public AtomStereoChangeInformation(IAtom rAtom, IAtom pAtom) {
        this.reactantAtom = rAtom;
        this.productAtom = pAtom;
        setStereoChange(true);
    }

    /**
     *
     * @param atomE
     * @param atomP
     * @param aStereoR
     * @param aStereoP
     */
    public AtomStereoChangeInformation(IAtom atomE, IAtom atomP, IStereoAndConformation aStereoR, IStereoAndConformation aStereoP) {
        this(atomE, atomP);
        this.atomStereoR = aStereoR;
        this.atomStereoP = aStereoP;
    }

    /**
     * @return the reactantAtom
     */
    public IAtom getReactantAtom() {
        return reactantAtom;
    }

    /**
     * @return the productAtom
     */
    public IAtom getProductAtom() {
        return productAtom;
    }

    /**
     * @return the stereoChange
     */
    public boolean isStereoChange() {
        return stereoChange;
    }

    /**
     * @param stereoChange the stereoChange to set
     */
    private void setStereoChange(boolean stereoChange) {
        this.stereoChange = stereoChange;
    }

    /**
     * @return the atomStereo
     */
    public IStereoAndConformation getReactantAtomStereo() {
        return atomStereoR;
    }

    /**
     * @return the atomStereo
     */
    public IStereoAndConformation getProductAtomStereo() {
        return atomStereoP;
    }

    /**
     * @param atomStereoR
     * @param atomStereoP
     */
    public void setAtomStereo(IStereoAndConformation atomStereoR, IStereoAndConformation atomStereoP) {
        this.atomStereoR = atomStereoR;
        this.atomStereoP = atomStereoP;
    }
}
