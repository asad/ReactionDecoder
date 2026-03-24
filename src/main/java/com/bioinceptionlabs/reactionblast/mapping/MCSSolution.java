/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping;

import java.io.Serializable;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class MCSSolution implements Serializable {

    private static final long serialVersionUID = 0xc678991ddf0L;
    private final IAtomContainer queryContainer;
    private final int targetPosition;
    private final IAtomContainer targetContainer;
    private final AtomAtomMapping atomatomMapping;
    private final int queryPosition;
    private Integer stereoScore;
    private Integer fragmentSize;
    private Double energy;

    /**
     *
     * @param queryPosition
     * @param targetPosition
     * @param queryContainer
     * @param targetContainer
     * @param aam
     */
    public MCSSolution(int queryPosition, int targetPosition,
            IAtomContainer queryContainer, IAtomContainer targetContainer, AtomAtomMapping aam) {
        this.queryContainer = queryContainer;
        this.targetPosition = targetPosition;
        this.targetContainer = targetContainer;
        this.atomatomMapping = aam;
        this.queryPosition = queryPosition;
        this.energy = null;
        this.fragmentSize = null;
        this.stereoScore = null;
    }

    /**
     * @return the stereoScore, or 0 if null
     */
    public Integer getStereoScore() {
        return stereoScore != null ? stereoScore : 0;
    }

    /**
     * @param stereoScore the stereoScore to set
     */
    public void setStereoScore(Integer stereoScore) {
        this.stereoScore = stereoScore;
    }

    /**
     * @return the fragmentSize, or 0 if null
     */
    public Integer getFragmentSize() {
        return fragmentSize != null ? fragmentSize : 0;
    }

    /**
     * @param fragmentSize the fragmentSize to set
     */
    public void setFragmentSize(Integer fragmentSize) {
        this.fragmentSize = fragmentSize;
    }

    /**
     * @return the energy, or 0.0 if null
     */
    public Double getEnergy() {
        return energy != null ? energy : 0.0;
    }

    /**
     * @param energy the energy to set
     */
    public void setEnergy(Double energy) {
        this.energy = energy;
    }

    /**
     * @return the queryContainer
     */
    public IAtomContainer getQueryContainer() {
        return queryContainer;
    }

    /**
     * @return the targetContainer
     */
    public IAtomContainer getTargetContainer() {
        return targetContainer;
    }

    /**
     * @return the atomatomMapping
     */
    public AtomAtomMapping getAtomAtomMapping() {
        return atomatomMapping;
    }

    /**
     * @return the targetPosition
     */
    public int getTargetPosition() {
        return targetPosition;
    }

    /**
     * @return the queryPosition
     */
    public int getQueryPosition() {
        return queryPosition;
    }
}
