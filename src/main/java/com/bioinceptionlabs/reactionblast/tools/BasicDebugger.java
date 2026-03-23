/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.tools;

import static java.lang.System.getProperty;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 */
public abstract class BasicDebugger {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(BasicDebugger.class);

    protected static final String NEW_LINE = getProperty("line.separator");

    /**
     *
     * @param mappings
     */
    public static void printAtomAtomMapping(Map<IAtom, IAtom> mappings) {
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        mappings.entrySet().stream().map((m) -> {
            sb.append("e:").append(m.getKey().getID()).append(NEW_LINE);
            return m;
        }).forEach((m) -> {
            sb.append("p:").append(m.getValue().getID()).append(NEW_LINE);
        });
        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param reaction
     */
    protected static void printReaction(IReaction reaction) {
        IAtomContainerSet Educt = reaction.getReactants();
        IAtomContainerSet Product = reaction.getProducts();

        StringBuilder sb = new StringBuilder();
        sb.append("*******************************").append(NEW_LINE);
        sb.append("Educt Mol Count: ").append(Educt.getAtomContainerCount()).append(NEW_LINE);
        sb.append("*******************************").append(NEW_LINE);

        for (int j = 0; j < Educt.getAtomContainerCount(); j++) {

            IAtomContainer M = Educt.getAtomContainer(j);
            sb.append("Mol ID: ").append(M.getID()).append(NEW_LINE);
            sb.append("SingleElectron: ").append(M.getSingleElectronCount()).append(NEW_LINE);
            sb.append("Stoic: ").append(reaction.getReactantCoefficient(M)).append(NEW_LINE);
            sb.append("Split Mol Atom Count: ").append(M.getAtomCount()).append(NEW_LINE);
            appendAtoms(sb, M);
        }

        sb.append("*******************************").append(NEW_LINE);
        sb.append("Product Mol Count: ").append(Product.getAtomContainerCount()).append(NEW_LINE);
        sb.append("*******************************").append(NEW_LINE);

        for (int j = 0; j < Product.getAtomContainerCount(); j++) {

            IAtomContainer M = Product.getAtomContainer(j);
            sb.append("Mol ID: ").append(M.getID()).append(NEW_LINE);
            sb.append("SingleElectron: ").append(M.getSingleElectronCount()).append(NEW_LINE);
            sb.append("Stoic: ").append(reaction.getProductCoefficient(M)).append(NEW_LINE);
            sb.append("Split Mol Atom Count: ").append(M.getAtomCount()).append(NEW_LINE);
            appendAtoms(sb, M);

        }

        sb.append(NEW_LINE).append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%").append(NEW_LINE).append(NEW_LINE);
        LOGGER.debug(sb.toString());
    }

    private static void appendAtoms(StringBuilder sb, IAtomContainer mol) {
        sb.append("Atom: ");
        for (IAtom a : mol.atoms()) {
            sb.append(a.getSymbol());
            if (a.getID() != null) {
                sb.append("[").append(a.getID()).append("]");
            }
        }
        sb.append(NEW_LINE).append(NEW_LINE);
    }

    /**
     * Print Atoms in molecules
     *
     * @param mol
     */
    protected static void printAtoms(IAtomContainer mol) {
        StringBuilder sb = new StringBuilder();
        sb.append("Atom: ");
        for (IAtom a : mol.atoms()) {

            sb.append(a.getSymbol());
            if (a.getID() != null) {
                sb.append("[").append(a.getID()).append("]");
            }

        }
        LOGGER.debug(sb.toString());
    }

    /**
     * Prints atoms in molecules
     *
     * @param molecule
     */
    protected static void printMolecule(IAtomContainer molecule) {

        StringBuilder sb = new StringBuilder();
        sb.append("AtomContainer ").append(molecule.getID()).append(": ").append(molecule.getAtomCount()).append(NEW_LINE);

        for (int i = 0; i < molecule.getAtomCount(); i++) {

            sb.append(molecule.getAtom(i).getSymbol()).append(" : ").append(molecule.getAtom(i).getID()).append(",  ");
        }

        LOGGER.debug(sb.toString());

    }
}
