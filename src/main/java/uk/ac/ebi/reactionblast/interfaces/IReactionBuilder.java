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
package uk.ac.ebi.reactionblast.interfaces;

import org.openscience.cdk.interfaces.IReactionSet;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public interface IReactionBuilder {

    /**
     *
     * @return get IReaction objects
     */
    IReactionSet getReactions();

    /**
     *process the reaction process
     */
    void process();

    /**
     *write the fingerprints in a file
     */
    void writeCompoundFingerPrints();

    /**
     *write the InChi into a file
     */
    void writeCompoundInChI();

    /**
     * write the reaction in a file
     */
    void writeReactions();
}
