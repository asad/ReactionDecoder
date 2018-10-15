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
package uk.ac.ebi.reactionblast.tools;

import static java.lang.System.currentTimeMillis;

import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mapping.container.CDKReactionBuilder;
import static uk.ac.ebi.reactionblast.mapping.helper.MappingHandler.cleanMapping;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class StandardizeReaction implements IStandardizer {

    /**
     *
     * @param reaction
     * @return New Standardized reaction Object
     * @throws Exception
     */
    @Override
    public synchronized IReaction standardize(IReaction reaction) throws Exception {
        String ReactionID = reaction.getID();
        cleanMapping(reaction);

        if (ReactionID == null) {
            ReactionID = Long.toString(currentTimeMillis());
            reaction.setID(ReactionID);
        }
        CDKReactionBuilder rBuilder = new CDKReactionBuilder();
        IReaction standardizedReaction = rBuilder.standardize(reaction);
        return standardizedReaction;
    }
}
