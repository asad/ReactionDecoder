/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.graphics.direct;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 * A union of sub-spanning trees over a molecular graph, each of which has a
 * root atom.
 *
 * @author maclean
 *
 */
public class RootSystem {

    /**
     * The atom roots of the trees that make up the system.
     */
    private final List<IAtom> roots;

    /**
     * All the other atoms in the system.
     */
    private final List<IAtom> leaves;

    /**
     *
     */
    public RootSystem() {
        this.roots = new ArrayList<>();
        this.leaves = new ArrayList<>();
    }

    /**
     *
     * @param root
     */
    public void addRoot(IAtom root) {
        if (roots.contains(root)) {
        } else {
            roots.add(root);
        }
    }

    /**
     *
     * @param bond
     */
    public void addRootsFromBond(IBond bond) {
        addRoot(bond.getAtom(0));
        addRoot(bond.getAtom(1));
    }

    /**
     *
     * @param leaf
     */
    public void addLeaf(IAtom leaf) {
        if (leaves.contains(leaf)) {
        } else {
            leaves.add(leaf);
        }
    }

    /**
     *
     * @return
     */
    public List<IAtom> getRoots() {
        return roots;
    }

    /**
     *
     * @return
     */
    public List<IAtom> getLeaves() {
        return leaves;
    }

    /**
     *
     * @param otherRootSystem
     * @return
     */
    public RootSystem merge(RootSystem otherRootSystem) {
        RootSystem merged = new RootSystem();
        merged.roots.addAll(roots);
        merged.roots.addAll(otherRootSystem.roots);
        merged.leaves.addAll(leaves);
        merged.leaves.addAll(otherRootSystem.leaves);
        return merged;
    }

    private void printAtomList(List<IAtom> atoms, StringBuilder sb) {
        sb.append("{");
        for (int index = 0; index < atoms.size(); index++) {
            IAtom root = atoms.get(index);
            sb.append(root.getID());
            if (index < atoms.size() - 1) {
                sb.append(",");
            }
        }
        sb.append("}");
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Roots ");
        printAtomList(roots, sb);
        sb.append(" Leaves ");
        printAtomList(leaves, sb);
        return sb.toString();
    }

}
