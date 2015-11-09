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

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

public class ReactionImplicitHydrogenTransformation implements
        ITransformation<IReaction> {

    @Override
    public ITransformation.TargetType getTargetType() {
        return TargetType.REACTION;
    }

    @Override
    public IReaction transform(IReaction reaction) {
        CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(
                DefaultChemObjectBuilder.getInstance());
        for (IAtomContainer atomContainer : ReactionManipulator.getAllAtomContainers(reaction)) {
            try {
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
                adder.addImplicitHydrogens(atomContainer);
            } catch (CDKException e) {
                e.printStackTrace();
            }
        }
        return reaction;
    }
}