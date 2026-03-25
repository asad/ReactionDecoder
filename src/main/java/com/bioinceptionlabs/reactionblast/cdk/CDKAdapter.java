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
package com.bioinceptionlabs.reactionblast.cdk;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;

import com.bioinceptionlabs.reactionblast.model.AtomNode;
import com.bioinceptionlabs.reactionblast.model.BondEdge;
import com.bioinceptionlabs.reactionblast.model.MolecularGraph;
import com.bioinceptionlabs.reactionblast.model.ReactionGraph;

/**
 * Bidirectional conversion between CDK types and graph model types.
 * Use during the Strangler Fig migration: existing CDK code wraps/unwraps
 * through this adapter as it's gradually refactored to use graph model types.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class CDKAdapter {

    private CDKAdapter() {}

    // ---- CDK → Graph Model ----

    public static MolecularGraph fromCDK(IAtomContainer mol) {
        return new CDKMolecularGraph(mol);
    }

    public static ReactionGraph fromCDK(IReaction rxn) {
        return new CDKReactionGraph(rxn);
    }

    public static AtomNode fromCDK(IAtom atom) {
        return new CDKAtomNode(atom);
    }

    public static BondEdge fromCDK(IBond bond) {
        return new CDKBondEdge(bond);
    }

    // ---- Graph Model → CDK ----

    public static IAtomContainer toCDK(MolecularGraph graph) {
        if (graph instanceof CDKMolecularGraph) {
            return ((CDKMolecularGraph) graph).getCDKContainer();
        }
        throw new IllegalArgumentException(
                "Cannot convert non-CDK MolecularGraph to IAtomContainer. "
                + "Use CDKToolkit for CDK-based operations.");
    }

    public static IReaction toCDK(ReactionGraph rxn) {
        if (rxn instanceof CDKReactionGraph) {
            return ((CDKReactionGraph) rxn).getCDKReaction();
        }
        throw new IllegalArgumentException(
                "Cannot convert non-CDK ReactionGraph to IReaction. "
                + "Use CDKToolkit for CDK-based operations.");
    }

    public static IAtom toCDK(AtomNode node) {
        if (node instanceof CDKAtomNode) {
            return ((CDKAtomNode) node).getCDKAtom();
        }
        throw new IllegalArgumentException(
                "Cannot convert non-CDK AtomNode to IAtom.");
    }

    public static IBond toCDK(BondEdge edge) {
        if (edge instanceof CDKBondEdge) {
            return ((CDKBondEdge) edge).getCDKBond();
        }
        throw new IllegalArgumentException(
                "Cannot convert non-CDK BondEdge to IBond.");
    }

    // ---- Type checking ----

    public static boolean isCDK(MolecularGraph graph) {
        return graph instanceof CDKMolecularGraph;
    }

    public static boolean isCDK(ReactionGraph rxn) {
        return rxn instanceof CDKReactionGraph;
    }
}
