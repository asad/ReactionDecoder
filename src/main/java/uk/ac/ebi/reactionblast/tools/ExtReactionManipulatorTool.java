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

import java.util.HashMap;
import java.util.Map;


import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ExtReactionManipulatorTool extends ReactionManipulator {

   
    /**
     *
     * @param reaction
     * @return deep clone of the reactions with mol IDs set and reaction ids set
     * plus flags copied
     * @throws CloneNotSupportedException
     */
    public static synchronized IReaction deepClone(IReaction reaction) throws CloneNotSupportedException {
        IReaction clone = new Reaction();
        // clone the reactants, products and agents

        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer acClone = new AtomContainer(ac).clone();
            /*Set IDs as CDK clone doesn't*/
            for (int i = 0; i < ac.getAtomCount(); i++) {
                acClone.getAtom(i).setID(ac.getAtom(i).getID());
            }
            acClone.setID(ac.getID());
            acClone.addProperties(ac.getProperties());
            clone.getReactants().addAtomContainer(acClone);
        }

        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer acClone = new AtomContainer(ac).clone();
            /*Set IDs as CDK clone doesn't*/
            for (int i = 0; i < ac.getAtomCount(); i++) {
                acClone.getAtom(i).setID(ac.getAtom(i).getID());
            }
            acClone.setID(ac.getID());
            acClone.addProperties(ac.getProperties());
            clone.getProducts().addAtomContainer(acClone);
        }

        for (IAtomContainer ac : reaction.getAgents().atomContainers()) {
            IAtomContainer acClone = new AtomContainer(ac).clone();
            acClone.setID(ac.getID());
            acClone.addProperties(ac.getProperties());
            clone.getAgents().addAtomContainer(acClone);
        }

        // create a Map of corresponding atoms for molecules (key: original Atom, 
        // value: clone Atom)
        Map<IChemObject, IChemObject> atomatom = new HashMap<>();
        for (int i = 0; i < reaction.getReactants().getAtomContainerCount(); ++i) {
            IAtomContainer mol = reaction.getReactants().getAtomContainer(i);
            IAtomContainer mol2 = clone.getReactants().getAtomContainer(i);
            for (int j = 0; j < mol.getAtomCount(); ++j) {
                atomatom.put(mol.getAtom(j), mol2.getAtom(j));
            }
        }
        for (int i = 0; i < reaction.getProducts().getAtomContainerCount(); ++i) {
            IAtomContainer mol = reaction.getProducts().getAtomContainer(i);
            IAtomContainer mol2 = clone.getProducts().getAtomContainer(i);
            for (int j = 0; j < mol.getAtomCount(); ++j) {
                atomatom.put(mol.getAtom(j), mol2.getAtom(j));
            }
        }
        //Add mapping to the clone
        for (IMapping mapping : reaction.mappings()) {
            clone.addMapping(new Mapping(atomatom.get(mapping.getChemObject(0)), atomatom.get(mapping.getChemObject(1))));

        }
        clone.setID(reaction.getID());
        return clone;
    }

    /**
     *
     * @param reaction
     * @return a new mol with explicit Hydrogens
     * @throws CloneNotSupportedException
     */
    public static IReaction addExplicitH(IReaction reaction) throws CloneNotSupportedException {
        IReaction r = reaction.getBuilder().newInstance(IReaction.class);
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer addExplicitH = ExtAtomContainerManipulator.addExplicitH(ac);
            r.addReactant(addExplicitH, reaction.getReactantCoefficient(ac));
        }

        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer addExplicitH = ExtAtomContainerManipulator.addExplicitH(ac);
            r.addProduct(addExplicitH, reaction.getProductCoefficient(ac));
        }

        r.setDirection(reaction.getDirection());
        r.setID(reaction.getID() == null ? "" : reaction.getID());

        return r;
    }
}
