/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.helper;

import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public class IOReaction {

    private static final Logger LOG = Logger.getLogger(IOReaction.class.getName());

    private static void printReaction(IReaction reaction) {
        IAtomContainerSet Educt = reaction.getReactants();
        IAtomContainerSet Product = reaction.getProducts();

        System.out.println("*******************************");

        System.out.println("Educt Mol Count: " + Educt.getAtomContainerCount());

        System.out.println("*******************************");

        for (int j = 0; j < Educt.getAtomContainerCount(); j++) {

            IAtomContainer M = Educt.getAtomContainer(j);
            System.out.println("Mol ID: " + M.getID());
            System.out.println("Stoic: " + reaction.getReactantCoefficient(M));
            System.out.println("Split Mol Atom Count: " + M.getAtomCount());

            printAtoms(M);

        }

        System.out.println("*******************************");

        System.out.println("Product Mol Count: " + Product.getAtomContainerCount());

        System.out.println("*******************************");

        for (int j = 0; j < Product.getAtomContainerCount(); j++) {

            IAtomContainer M = Product.getAtomContainer(j);
            System.out.println("Mol ID: " + M.getID());
            System.out.println("Stoic: " + reaction.getProductCoefficient(M));
            System.out.println("Split Mol Atom Count: " + M.getAtomCount());

            printAtoms(M);

        }

        System.out.println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
    }

    private static void printAtoms(IAtomContainer mol) {
        System.out.print("Atom: ");
        for (IAtom a : mol.atoms()) {

            System.out.print(a.getSymbol());

        }
        System.out.println();
        System.out.println();
    }
}
