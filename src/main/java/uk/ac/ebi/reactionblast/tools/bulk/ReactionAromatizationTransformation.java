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
package uk.ac.ebi.reactionblast.tools.bulk;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;
import static uk.ac.ebi.reactionblast.interfaces.ITransformation.TargetType.REACTION;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeDayLight;

/**
 *
 * @author asad
 */
public class ReactionAromatizationTransformation implements
        ITransformation<IReaction> {

    /**
     *
     * @return
     */
    @Override
    public ITransformation.TargetType getTargetType() {
        return REACTION;
    }

    /**
     *
     * @param reaction
     * @return
     */
    @Override
    public IReaction transform(IReaction reaction) {
        for (IAtomContainer atomContainer : getAllAtomContainers(reaction)) {
            try {
                aromatizeDayLight(atomContainer);
            } catch (CDKException e) {
                e.printStackTrace();
            }
        }
        return reaction;
    }
}
