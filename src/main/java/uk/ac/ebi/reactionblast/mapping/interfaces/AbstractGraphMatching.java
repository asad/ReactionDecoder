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

import java.util.BitSet;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class AbstractGraphMatching {

    /**
     *
     * @param holder
     * @param substrateIndex
     * @param productIndex
     * @throws java.lang.Exception
     */
    public static void setMCSUpdationFlags(Holder holder, int substrateIndex, int productIndex) throws Exception {
        ReactionContainer reactionStructureInformation = holder.getReactionContainer();
        reactionStructureInformation.setEductModified(substrateIndex, true);
        reactionStructureInformation.setProductModified(productIndex, true);
    }

    /**
     * @return the matchedPart
     */
    public abstract IAtomContainer getMatchedPart();

    /**
     *
     * @return
     */
    public abstract IAtomContainer getRemainingEduct();

    /**
     *
     * @return
     */
    public abstract IAtomContainer getRemainingProduct();

    /**
     *
     * @param holder Data holder
     * @param removeHydrogen
     * @param I
     * @param J
     * @param eductFP
     * @param prodFP
     * @return
     */
    public abstract boolean mcsMatch(Holder holder, boolean removeHydrogen, Integer I, Integer J, BitSet eductFP, BitSet prodFP);

    /**
     * Removed matched part from the reaction
     *
     * @param reaction
     * @return
     */
    public abstract int removeMatchedAtomsAndUpdateAAM(IReaction reaction);
}
