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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.model.AtomNode;
import com.bioinceptionlabs.reactionblast.model.MolecularGraph;
import com.bioinceptionlabs.reactionblast.model.ReactionGraph;

/**
 * CDK adapter for ReactionGraph. Wraps a CDK IReaction.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CDKReactionGraph implements ReactionGraph {

    private final IReaction cdkReaction;

    public CDKReactionGraph(IReaction cdkReaction) {
        if (cdkReaction == null) throw new IllegalArgumentException("Reaction cannot be null");
        this.cdkReaction = cdkReaction;
    }

    public IReaction getCDKReaction() {
        return cdkReaction;
    }

    @Override public int getReactantCount() { return cdkReaction.getReactantCount(); }
    @Override public int getProductCount() { return cdkReaction.getProductCount(); }

    @Override
    public MolecularGraph getReactant(int index) {
        return new CDKMolecularGraph(cdkReaction.getReactants().getAtomContainer(index));
    }

    @Override
    public MolecularGraph getProduct(int index) {
        return new CDKMolecularGraph(cdkReaction.getProducts().getAtomContainer(index));
    }

    @Override
    public Iterable<MolecularGraph> getReactants() {
        List<MolecularGraph> result = new ArrayList<>();
        for (IAtomContainer ac : cdkReaction.getReactants().atomContainers()) {
            result.add(new CDKMolecularGraph(ac));
        }
        return result;
    }

    @Override
    public Iterable<MolecularGraph> getProducts() {
        List<MolecularGraph> result = new ArrayList<>();
        for (IAtomContainer ac : cdkReaction.getProducts().atomContainers()) {
            result.add(new CDKMolecularGraph(ac));
        }
        return result;
    }

    @Override
    public void addReactant(MolecularGraph mol) {
        if (mol instanceof CDKMolecularGraph) {
            cdkReaction.addReactant(((CDKMolecularGraph) mol).getCDKContainer());
        }
    }

    @Override
    public void addProduct(MolecularGraph mol) {
        if (mol instanceof CDKMolecularGraph) {
            cdkReaction.addProduct(((CDKMolecularGraph) mol).getCDKContainer());
        }
    }

    @Override public String getId() { return cdkReaction.getID(); }
    @Override public void setId(String id) { cdkReaction.setID(id); }

    @Override
    public Map<AtomNode, AtomNode> getAtomMapping() {
        Map<AtomNode, AtomNode> mapping = new HashMap<>();
        for (IMapping m : cdkReaction.mappings()) {
            IAtom a1 = (IAtom) m.getChemObject(0);
            IAtom a2 = (IAtom) m.getChemObject(1);
            if (a1 != null && a2 != null) {
                mapping.put(new CDKAtomNode(a1), new CDKAtomNode(a2));
            }
        }
        return mapping;
    }

    @Override
    public void setAtomMapping(Map<AtomNode, AtomNode> mapping) {
        // Clear existing mappings
        int count = cdkReaction.getMappingCount();
        for (int i = count - 1; i >= 0; i--) {
            cdkReaction.removeMapping(i);
        }
        // Add new ones
        for (Map.Entry<AtomNode, AtomNode> entry : mapping.entrySet()) {
            if (entry.getKey() instanceof CDKAtomNode && entry.getValue() instanceof CDKAtomNode) {
                IAtom a1 = ((CDKAtomNode) entry.getKey()).getCDKAtom();
                IAtom a2 = ((CDKAtomNode) entry.getValue()).getCDKAtom();
                cdkReaction.addMapping(new org.openscience.cdk.Mapping(a1, a2));
            }
        }
    }

    @Override
    public boolean isMapped() {
        return cdkReaction.getMappingCount() > 0;
    }

    @Override
    public ReactionGraph clone() throws CloneNotSupportedException {
        return new CDKReactionGraph((IReaction) cdkReaction.clone());
    }

    @Override
    public String toString() {
        return "CDKReactionGraph{" + getId() + ", R=" + getReactantCount() + ", P=" + getProductCount() + "}";
    }
}
