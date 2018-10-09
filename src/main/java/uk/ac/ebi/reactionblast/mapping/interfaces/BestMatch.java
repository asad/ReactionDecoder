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
package uk.ac.ebi.reactionblast.mapping.interfaces;

import java.io.IOException;
import org.openscience.smsd.AtomAtomMapping;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class BestMatch {

    /**
     *
     */
    public BestMatch() {
    }

    /**
     *
     * @throws java.io.IOException
     */
    public abstract void Clear() throws IOException;

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @throws java.io.IOException
     */
    public abstract void Erase(int indexI, int indexJ) throws IOException;

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     * @throws IOException
     */
    public abstract boolean containsKey(int indexI, int indexJ) throws IOException;

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     * @throws IOException
     */
    public abstract AtomAtomMapping getAtomMatch(int indexI, int indexJ) throws IOException;

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    public abstract double getBondEnergy(int indexI, int indexJ);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    public abstract double getGraphSimilarity(int indexI, int indexJ);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    public abstract double getStereoScore(int indexI, int indexJ);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    public abstract int getTotalFragmentCount(int indexI, int indexJ);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param matchingAtoms
     * @throws IOException
     */
    public abstract void putBestMapping(int indexI, int indexJ, AtomAtomMapping matchingAtoms) throws IOException;

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param value
     */
    public abstract void setBondEnergy(int indexI, int indexJ, double value);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param _simVal
     */
    public abstract void setGraphSimilarity(int indexI, int indexJ, double _simVal);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param stereoVal
     */
    public abstract void setStereoScore(int indexI, int indexJ, double stereoVal);

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param value
     */
    public abstract void setTotalFragmentCount(int indexI, int indexJ, Integer value);
}
