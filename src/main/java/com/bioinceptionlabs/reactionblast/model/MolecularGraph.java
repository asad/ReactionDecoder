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

import java.util.List;
import java.util.Map;

/**
 * Labeled molecular graph — nodes are atoms, edges are bonds.
 * Toolkit-agnostic — implementations wrap CDK IAtomContainer,
 * RDKit RWMol, OpenBabel OBMol, etc.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface MolecularGraph {

    int getNodeCount();

    int getEdgeCount();

    AtomNode getNode(int index);

    BondEdge getEdge(int index);

    int indexOf(AtomNode node);

    Iterable<AtomNode> nodes();

    Iterable<BondEdge> edges();

    List<BondEdge> getEdges(AtomNode node);

    List<AtomNode> getNeighbors(AtomNode node);

    BondEdge getEdge(AtomNode a, AtomNode b);

    String getId();

    void setId(String id);

    Object getProperty(String key);

    void setProperty(String key, Object value);

    MolecularGraph clone() throws CloneNotSupportedException;

    void addNode(AtomNode node);

    void addEdge(BondEdge edge);

    void removeNode(AtomNode node);

    void removeEdge(BondEdge edge);
}
