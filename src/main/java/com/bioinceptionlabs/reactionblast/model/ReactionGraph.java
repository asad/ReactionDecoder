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
 * Reaction as a graph transformation: reactant graphs → product graphs
 * with atom-atom mapping between them. Toolkit-agnostic.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface ReactionGraph {

    int getReactantCount();

    int getProductCount();

    MolecularGraph getReactant(int index);

    MolecularGraph getProduct(int index);

    Iterable<MolecularGraph> getReactants();

    Iterable<MolecularGraph> getProducts();

    void addReactant(MolecularGraph mol);

    void addProduct(MolecularGraph mol);

    String getId();

    void setId(String id);

    Map<AtomNode, AtomNode> getAtomMapping();

    void setAtomMapping(Map<AtomNode, AtomNode> mapping);

    boolean isMapped();

    ReactionGraph clone() throws CloneNotSupportedException;
}
