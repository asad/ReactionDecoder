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

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import com.bioinceptionlabs.reactionblast.model.AtomNode;
import com.bioinceptionlabs.reactionblast.model.BondEdge;
import com.bioinceptionlabs.reactionblast.model.MolecularGraph;

/**
 * CDK adapter for MolecularGraph. Wraps a CDK IAtomContainer as a labeled graph.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CDKMolecularGraph implements MolecularGraph {

    private final IAtomContainer container;

    public CDKMolecularGraph(IAtomContainer container) {
        if (container == null) throw new IllegalArgumentException("Container cannot be null");
        this.container = container;
    }

    public IAtomContainer getCDKContainer() {
        return container;
    }

    @Override public int getNodeCount() { return container.getAtomCount(); }
    @Override public int getEdgeCount() { return container.getBondCount(); }

    @Override
    public AtomNode getNode(int index) {
        return new CDKAtomNode(container.getAtom(index));
    }

    @Override
    public BondEdge getEdge(int index) {
        return new CDKBondEdge(container.getBond(index));
    }

    @Override
    public int indexOf(AtomNode node) {
        if (node instanceof CDKAtomNode) {
            return container.indexOf(((CDKAtomNode) node).getCDKAtom());
        }
        return -1;
    }

    @Override
    public Iterable<AtomNode> nodes() {
        List<AtomNode> nodes = new ArrayList<>(container.getAtomCount());
        for (IAtom atom : container.atoms()) {
            nodes.add(new CDKAtomNode(atom));
        }
        return nodes;
    }

    @Override
    public Iterable<BondEdge> edges() {
        List<BondEdge> edges = new ArrayList<>(container.getBondCount());
        for (IBond bond : container.bonds()) {
            edges.add(new CDKBondEdge(bond));
        }
        return edges;
    }

    @Override
    public List<BondEdge> getEdges(AtomNode node) {
        List<BondEdge> result = new ArrayList<>();
        if (node instanceof CDKAtomNode) {
            IAtom cdkAtom = ((CDKAtomNode) node).getCDKAtom();
            for (IBond bond : container.getConnectedBondsList(cdkAtom)) {
                result.add(new CDKBondEdge(bond));
            }
        }
        return result;
    }

    @Override
    public List<AtomNode> getNeighbors(AtomNode node) {
        List<AtomNode> result = new ArrayList<>();
        if (node instanceof CDKAtomNode) {
            IAtom cdkAtom = ((CDKAtomNode) node).getCDKAtom();
            for (IAtom neighbor : container.getConnectedAtomsList(cdkAtom)) {
                result.add(new CDKAtomNode(neighbor));
            }
        }
        return result;
    }

    @Override
    public BondEdge getEdge(AtomNode a, AtomNode b) {
        if (a instanceof CDKAtomNode && b instanceof CDKAtomNode) {
            IBond bond = container.getBond(
                    ((CDKAtomNode) a).getCDKAtom(),
                    ((CDKAtomNode) b).getCDKAtom());
            return bond != null ? new CDKBondEdge(bond) : null;
        }
        return null;
    }

    @Override public String getId() { return container.getID(); }
    @Override public void setId(String id) { container.setID(id); }
    @Override public Object getProperty(String key) { return container.getProperty(key); }
    @Override public void setProperty(String key, Object value) { container.setProperty(key, value); }

    @Override
    public MolecularGraph clone() throws CloneNotSupportedException {
        return new CDKMolecularGraph((IAtomContainer) container.clone());
    }

    @Override
    public void addNode(AtomNode node) {
        if (node instanceof CDKAtomNode) {
            container.addAtom(((CDKAtomNode) node).getCDKAtom());
        }
    }

    @Override
    public void addEdge(BondEdge edge) {
        if (edge instanceof CDKBondEdge) {
            container.addBond(((CDKBondEdge) edge).getCDKBond());
        }
    }

    @Override
    public void removeNode(AtomNode node) {
        if (node instanceof CDKAtomNode) {
            container.removeAtom(((CDKAtomNode) node).getCDKAtom());
        }
    }

    @Override
    public void removeEdge(BondEdge edge) {
        if (edge instanceof CDKBondEdge) {
            container.removeBond(((CDKBondEdge) edge).getCDKBond());
        }
    }

    @Override
    public String toString() {
        return "CDKMolecularGraph{" + getId() + ", atoms=" + getNodeCount() + ", bonds=" + getEdgeCount() + "}";
    }
}
