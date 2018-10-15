/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.blocks;

import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import signature.AbstractGraphSignature;
import signature.AbstractVertexSignature;

/**
 * A specialized signature that covers only a (connected) part of a molecule.
 * This is different to the idea of a signature of a certain height, since the
 * subgraph may be of varying heights. It is really no different to extracting
 * the atoms and bonds from the original molecule, and making a signature from
 * that - just more convenient.
 *
 * @author maclean
 *
 */
public class SubgraphMoleculeSignature extends AbstractGraphSignature {

    private final int vertexCount;

    private final IAtomContainer atomContainer;

    private Map<Integer, int[]> subgraphAdjacencyLists;

    // XXX FIXME - this is necessary because the canonise routine will
    // run over the indices (0 to subgraphSize) without working out what
    // the actual atomIndex is...
    private List<Integer> indexLookup;

    // XXX FIXME : avoids the problem of subgraphs not working... :( 
    private final boolean useAdjLists;

    /**
     *
     * @param fullContainer
     * @param subgraphAtoms
     * @param dummy
     */
    public SubgraphMoleculeSignature(
            IAtomContainer fullContainer, List<IAtom> subgraphAtoms, int dummy) {
        this.atomContainer
                = fullContainer.getBuilder().newInstance(IAtomContainer.class);

//        for (IAtom atom : subgraphAtoms) {
//            this.atomContainer.addAtom(atom);
//        }
        // sort the atoms
        IAtom[] sortedAtoms = new IAtom[fullContainer.getAtomCount()];
        List<String> debugList = new ArrayList<>();
        subgraphAtoms.stream().forEach((atom) -> {
            int atomNumber = fullContainer.indexOf(atom);
            sortedAtoms[atomNumber] = atom;
            debugList.add(atom.getSymbol() + atomNumber);
        });
        sort(debugList);
//        System.out.println("atoms of " + fullContainer.getID() + " " + debugList);

        List<Integer> debugList2 = new ArrayList<>();
        for (int i = 0; i < sortedAtoms.length; i++) {
            if (sortedAtoms[i] != null) {
                this.atomContainer.addAtom(sortedAtoms[i]);
                debugList2.add(i);
            }
        }
//        System.out.println("atoms (2) of " + atomContainer.getID() + " " + debugList2);

        for (IBond bond : fullContainer.bonds()) {
            IAtom a0 = bond.getAtom(0);
            IAtom a1 = bond.getAtom(1);
            int a0n = fullContainer.indexOf(a0);
            int a1n = fullContainer.indexOf(a1);
            if (subgraphAtoms.contains(a0) && subgraphAtoms.contains(a1)) {
                this.atomContainer.addBond(bond);
//                System.out.println("adding bond " + a0n + " " + a1n);
            } else {
//                System.out.println("rejecting bond " + a0n + " " + a1n);
            }
        }

        this.vertexCount = subgraphAtoms.size();
        useAdjLists = false;
    }

    /**
     *
     * @param atomContainer
     * @param subgraphAtoms
     */
    public SubgraphMoleculeSignature(
            IAtomContainer atomContainer, List<IAtom> subgraphAtoms) {

        this.atomContainer = atomContainer;

        this.indexLookup = new ArrayList<>();

        // create a lookup table that contains, for each atom index, the list
        // of connected atom indices that are also in the subgraph
        subgraphAdjacencyLists = new HashMap<>();
        for (int i = 0; i < atomContainer.getAtomCount(); i++) {
            IAtom atom = atomContainer.getAtom(i);
            if (subgraphAtoms.contains(atom)) {
                List<Integer> connectedIndices = new ArrayList<>();
                atomContainer.getConnectedAtomsList(atom).stream().filter((neighbour)
                        -> (subgraphAtoms.contains(neighbour))).map((neighbour)
                        -> atomContainer.indexOf(neighbour)).forEach(connectedIndices::add);
                int[] connectedIndicesArray = new int[connectedIndices.size()];
                int x = 0;
                for (int index : connectedIndices) {
                    connectedIndicesArray[x++] = index;
                }
                subgraphAdjacencyLists.put(i, connectedIndicesArray);
                indexLookup.add(i);
            }
        }

        // this is correct, but slightly strange
//        vertexCount = subgraphAtoms.size();
        vertexCount = atomContainer.getAtomCount();
        useAdjLists = true;// XXX!!
    }

    /**
     *
     * @return
     */
    @Override
    protected int getVertexCount() {
        return vertexCount;
    }

    /**
     *
     * @param subgraphAtomIndex
     * @return
     */
    @Override
    public AbstractVertexSignature signatureForVertex(int subgraphAtomIndex) {
        if (useAdjLists) {
            int actualIndex = indexLookup.get(subgraphAtomIndex);   // XXX FIXME!

            return new SubgraphAtomSignature(
                    atomContainer, actualIndex, subgraphAdjacencyLists);
        } else {
            return new SubgraphAtomSignature(atomContainer, subgraphAtomIndex);
        }
    }

    /**
     *
     * @param subgraphAtomIndex
     * @return
     */
    @Override
    public String signatureStringForVertex(int subgraphAtomIndex) {
        return signatureForVertex(subgraphAtomIndex).toCanonicalString();
    }

    /**
     *
     * @param atomIndex
     * @param height
     * @return
     */
    @Override
    public String signatureStringForVertex(int atomIndex, int height) {
        // TODO Auto-generated method stub
        return null;
    }

}
