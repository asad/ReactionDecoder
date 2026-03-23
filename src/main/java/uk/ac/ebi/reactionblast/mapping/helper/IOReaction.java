/*
 * Copyright (C) 2003-2020 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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

import static java.lang.System.getProperty;

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
public class IOReaction {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(IOReaction.class);

    static final String NEW_LINE = getProperty("line.separator");

    private static void printReaction(IReaction reaction) {
        IAtomContainerSet Educt = reaction.getReactants();
        IAtomContainerSet Product = reaction.getProducts();

        StringBuilder sb = new StringBuilder();
        sb.append("*******************************").append(NEW_LINE);
        sb.append("Educt Mol Count: ").append(Educt.getAtomContainerCount()).append(NEW_LINE);
        sb.append("*******************************").append(NEW_LINE);

        for (int j = 0; j < Educt.getAtomContainerCount(); j++) {

            IAtomContainer M = Educt.getAtomContainer(j);
            sb.append("Mol ID: ").append(M.getID()).append(NEW_LINE);
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
        }
        sb.append(NEW_LINE).append(NEW_LINE);
    }

    private IOReaction() {
    }
}
