/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.blocks;

import static java.lang.String.format;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author asad
 */
public class DefinedMapping {

    private final IAtomContainer rAtomContainer;

    private final IAtomContainer pAtomContainer;

    private final int rIndex;

    private final int pIndex;

    private boolean visited = false;

    private final int index;

    /**
     *
     * @param rIndex
     * @param pIndex
     * @param index
     * @param rAtomContainer
     * @param pAtomContainer
     */
    public DefinedMapping(int rIndex, int pIndex, int index,
            IAtomContainer rAtomContainer, IAtomContainer pAtomContainer) {
        this.rIndex = rIndex;
        this.pIndex = pIndex;
        this.index = index;
        this.rAtomContainer = rAtomContainer;
        this.pAtomContainer = pAtomContainer;
    }

    /**
     *
     * @return
     */
    public IAtom getRAtom() {
        return getrAtomContainer().getAtom(getrIndex());
    }

    /**
     *
     * @return
     */
    public List<IAtom> getRAtomNeighbours() {
        return getrAtomContainer().getConnectedAtomsList(getRAtom());
    }

    /**
     *
     * @return
     */
    public IAtom getPAtom() {
        return getpAtomContainer().getAtom(getpIndex());
    }

    /**
     *
     * @return
     */
    public List<IAtom> getPAtomNeighbours() {
        return getpAtomContainer().getConnectedAtomsList(getPAtom());
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof DefinedMapping) {
            DefinedMapping other = (DefinedMapping) o;
            return getrAtomContainer() == other.getrAtomContainer()
                    && getpAtomContainer() == other.getpAtomContainer()
                    && getrIndex() == other.getrIndex()
                    && getpIndex() == other.getpIndex();
        }
        return false;
    }

    @Override
    public int hashCode() {
        return getrAtomContainer().hashCode() ^ getpAtomContainer().hashCode() * getrIndex() * getpIndex();
    }

    @Override
    public String toString() {
        return format("rI:%s pI:%s i:%s rID:%s pID:%s", getrIndex(), getpIndex(), getIndex(), getrAtomContainer().getID(), getpAtomContainer().getID());
    }

    /**
     * @return the visited
     */
    public boolean isVisited() {
        return visited;
    }

    /**
     * @param visited the visited to set
     */
    public void setVisited(boolean visited) {
        this.visited = visited;
    }

    /**
     * @return the rAtomContainer
     */
    public IAtomContainer getrAtomContainer() {
        return rAtomContainer;
    }

    /**
     * @return the pAtomContainer
     */
    public IAtomContainer getpAtomContainer() {
        return pAtomContainer;
    }

    /**
     * @return the rIndex
     */
    public int getrIndex() {
        return rIndex;
    }

    /**
     * @return the pIndex
     */
    public int getpIndex() {
        return pIndex;
    }

    /**
     * @return the index
     */
    public int getIndex() {
        return index;
    }
}
