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
package com.bioinceptionlabs.reactionblast.model;

import java.util.Map;

/**
 * Toolkit adapter interface + global registry.
 * Implement this for CDK, RDKit, OpenBabel, etc.
 *
 * Usage:
 * <pre>
 *   ChemToolkit.register(new CDKToolkit());            // once at startup
 *   ReactionGraph rxn = ChemToolkit.get().parseReactionSmiles("CC>>CC");
 * </pre>
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface ChemToolkit {

    // ---- Parsing ----

    ReactionGraph parseReactionSmiles(String smiles);

    MolecularGraph parseMoleculeSmiles(String smiles);

    // ---- Serialization ----

    String toSmiles(MolecularGraph mol);

    String toSmiles(ReactionGraph rxn);

    String toCanonicalSmiles(MolecularGraph mol);

    // ---- Perception ----

    void perceiveAtomTypes(MolecularGraph mol);

    void perceiveAromaticity(MolecularGraph mol);

    void addImplicitHydrogens(MolecularGraph mol);

    // ---- Substructure / MCS ----

    boolean isSubstructure(MolecularGraph query, MolecularGraph target);

    Map<AtomNode, AtomNode> findMCS(MolecularGraph mol1, MolecularGraph mol2);

    // ---- Factory methods ----

    MolecularGraph createMolecularGraph();

    AtomNode createAtomNode(String symbol);

    BondEdge createBondEdge(AtomNode source, AtomNode target, BondEdge.BondOrder order);

    ReactionGraph createReactionGraph();

    // ---- Global registry ----

    static ChemToolkit get() {
        return ChemToolkitRegistry.INSTANCE;
    }

    static void register(ChemToolkit toolkit) {
        ChemToolkitRegistry.INSTANCE = toolkit;
    }
}

/**
 * Internal holder for the global toolkit singleton.
 */
class ChemToolkitRegistry {
    static volatile ChemToolkit INSTANCE;
}
