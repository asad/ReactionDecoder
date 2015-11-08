/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;

/**
 * Just sets the atom ID to the mapping ID.
 * 
 * @author maclean
 *
 */
public class ReactionAtomIDTransformation implements ITransformation<IReaction> {

    @Override
    public ITransformation.TargetType getTargetType() {
        return TargetType.REACTION;
    }

    @Override
    public IReaction transform(IReaction reaction) {
        for (IAtomContainer atomContainer : ReactionManipulator.getAllAtomContainers(reaction)) {
            for (IAtom atom : atomContainer.atoms()) {
                Object prop = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
                if (prop != null) {
                    int mappingID = (Integer) prop;
                    atom.setID(String.valueOf(mappingID));
                }
            }
        }

        return reaction;
    }
    private static final Logger LOG = Logger.getLogger(ReactionAtomIDTransformation.class.getName());
}