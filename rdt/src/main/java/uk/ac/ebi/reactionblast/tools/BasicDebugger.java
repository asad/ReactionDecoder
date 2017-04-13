/*
 * Copyright (C) 2007-2017 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

import static java.lang.System.out;
import java.util.Map;
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
public abstract class BasicDebugger {

    /**
     *
     * @param mappings
     */
    public void printAtomAtomMapping(Map<IAtom, IAtom> mappings) {
        out.println();
        mappings.entrySet().stream().map((m) -> {
            out.println("e:" + m.getKey().getID());
            return m;
        }).forEach((m) -> {
            out.println("p:" + m.getValue().getID());
        });
    }

    /**
     *
     * @param reaction
     */
    protected void printReaction(IReaction reaction) {
        IAtomContainerSet Educt = reaction.getReactants();
        IAtomContainerSet Product = reaction.getProducts();

        out.println("*******************************");

        out.println("Educt Mol Count: " + Educt.getAtomContainerCount());

        out.println("*******************************");

        for (int j = 0; j < Educt.getAtomContainerCount(); j++) {

            IAtomContainer M = Educt.getAtomContainer(j);
            out.println("Mol ID: " + M.getID());
            out.println("SingleElectron: " + M.getSingleElectronCount());
            out.println("Stoic: " + reaction.getReactantCoefficient(M));
            out.println("Split Mol Atom Count: " + M.getAtomCount());
            printAtoms(M);
        }

        out.println("*******************************");

        out.println("Product Mol Count: " + Product.getAtomContainerCount());

        out.println("*******************************");

        for (int j = 0; j < Product.getAtomContainerCount(); j++) {

            IAtomContainer M = Product.getAtomContainer(j);
            out.println("Mol ID: " + M.getID());
            out.println("SingleElectron: " + M.getSingleElectronCount());
            out.println("Stoic: " + reaction.getProductCoefficient(M));
            out.println("Split Mol Atom Count: " + M.getAtomCount());
            printAtoms(M);

        }

        out.println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
    }

    /**
     * Print Atoms in molecules
     *
     * @param mol
     */
    protected void printAtoms(IAtomContainer mol) {
        out.print("Atom: ");
        for (IAtom a : mol.atoms()) {

            out.print(a.getSymbol());
            if (a.getID() != null) {
                out.print("[" + a.getID() + "]");
            }

        }
        out.println();
        out.println();
    }

    /**
     * Prints atoms in molecules
     *
     * @param molecule
     */
    protected void printMolecule(IAtomContainer molecule) {

        out.println("AtomContainer " + molecule.getID() + ": " + molecule.getAtomCount());

        for (int i = 0; i < molecule.getAtomCount(); i++) {

            out.print(molecule.getAtom(i).getSymbol() + " : " + molecule.getAtom(i).getID() + ",  ");
        }

        out.println();

    }
}
