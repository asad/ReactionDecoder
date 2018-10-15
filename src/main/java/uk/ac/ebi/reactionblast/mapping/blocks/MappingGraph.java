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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;

/**
 * A graph of the atom-atom mappings in a reaction - each vertex of the graph is
 * a mapping between a pair of atoms, and an edge connects mappings whose
 * endpoints (atoms) are connected (bonded).
 *
 * @author maclean
 *
 */
public class MappingGraph {

    /**
     *
     */
    public List<DefinedMapping> vertices;

    /**
     *
     */
    public List<DefinedMapping>[] adjacencyTable;

    /**
     *
     * @param reaction
     */
    public MappingGraph(IReaction reaction) {
        vertices = createDefinedMappings(reaction);
//        System.out.println(vertices.size() + " vertices");
        adjacencyTable = makeAdjacencyTable();
    }

    /**
     *
     * @param mappings
     */
    public MappingGraph(List<DefinedMapping> mappings) {
        vertices = mappings;
        adjacencyTable = makeAdjacencyTable();
    }

    /**
     *
     * @return
     */
    public List<List<DefinedMapping>> calculateConnectedComponents() {
        int[] componentLabels = new int[vertices.size()];
        List<List<DefinedMapping>> components
                = new ArrayList<>();
        int currentLabel = 1;

        int i = 0;
        while (i < vertices.size()) {
            List<DefinedMapping> component
                    = new ArrayList<>();
            DefinedMapping vertex = vertices.get(i);
            search(vertex, currentLabel, componentLabels, component);
            if (component.size() > 0) {
                components.add(component);
            }
            while (i < vertices.size() && componentLabels[i] != 0) {
                i++;
            }
            currentLabel++;
        }
//        System.out.println(components.size() + " components");
        return components;
    }

    private List<DefinedMapping> createDefinedMappings(IReaction reaction) {
        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet products = reaction.getProducts();

        List<DefinedMapping> definedMappings = new ArrayList<>();
        int i = 0;
        for (IMapping mapping : reaction.mappings()) {
            String id = mapping.getChemObject(0).getID();
            AtomContainerAtomPair reactantPair = getByID(reactants, id);
            AtomContainerAtomPair productPair = getByID(products, id);

            if (reactantPair != null && productPair != null) {
                int rIndex = reactantPair.getIndex();
                int pIndex = productPair.getIndex();
                definedMappings.add(
                        new DefinedMapping(rIndex, pIndex, i,
                                reactantPair.atomContainer, productPair.atomContainer));
                i++;
            }
        }
        return definedMappings;
    }

    /**
     *
     * @param reaction
     * @return
     */
    public List<BlockPair> createBlockPairs(IReaction reaction) {
        List<BlockPair> blockPairs = new ArrayList<>();

        calculateConnectedComponents().stream().map((mappingComponent) -> {
            // all components have at least one member
            DefinedMapping aMapping = mappingComponent.get(0);
            // initialise with this member
            BlockPair blockPair
                    = new BlockPair(aMapping.getrAtomContainer(), aMapping.getpAtomContainer());
            // add the mappings
            for (int i = 0; i < mappingComponent.size(); i++) {
                DefinedMapping definedMapping = mappingComponent.get(i);
                IMapping mapping = reaction.getMapping(definedMapping.getIndex());
                blockPair.addMapping(mapping,
                        definedMapping.getRAtom(), definedMapping.getPAtom());
            }
            return blockPair;
        }).forEach((blockPair) -> {
            blockPairs.add(blockPair);
        });

        return blockPairs;
    }

    private AtomContainerAtomPair getByID(IAtomContainerSet moleculeSet, String id) {
//        System.out.println("getting id " + id);
        for (IAtomContainer ac : moleculeSet.atomContainers()) {
            for (IAtom atom : ac.atoms()) {
                String atomID = atom.getID();
                if (atomID != null && atomID.equals(id)) {
                    int index = ac.indexOf(atom);
                    return new AtomContainerAtomPair(ac, atom, index);
                }
            }
        }
        return null;
    }

    private List<DefinedMapping>[] makeAdjacencyTable() {
        Map<IAtom, DefinedMapping> lookup;
        lookup = new HashMap<>();

        vertices.stream().map((mapping) -> {
            lookup.put(mapping.getRAtom(), mapping);
            return mapping;
        }).forEach((mapping) -> {
            lookup.put(mapping.getPAtom(), mapping);
        });

        List<DefinedMapping>[] adjTable;
        adjTable = new List[vertices.size()];
        for (int i = 0; i < vertices.size(); i++) {
            DefinedMapping vertexI = vertices.get(i);
            List<IAtom> pNeighbours = vertexI.getPAtomNeighbours();
            List<IAtom> rNeighbours = vertexI.getRAtomNeighbours();
            List<DefinedMapping> neighbours;
            if (adjTable[i] == null) {
                neighbours = new ArrayList<>();
                adjTable[i] = neighbours;
            } else {
                neighbours = adjTable[i];
            }
            for (int j = i + 1; j < vertices.size(); j++) {
                DefinedMapping vertexJ = vertices.get(j);
                if (pNeighbours.contains(vertexJ.getPAtom())
                        && rNeighbours.contains(vertexJ.getRAtom())) {
                    neighbours.add(vertexJ);
                    if (adjTable[j] == null) {
                        adjTable[j] = new ArrayList<>();
                    }
                    adjTable[j].add(vertexI);
                }
            }
        }
//        System.out.println(java.util.Arrays.deepToString(adjacencyTable));
        return adjTable;
    }

    /**
     * Search for connected components of the mapping graph.
     *
     * @param vertex the current vertex
     * @param currentLabel the label of the current connected component
     * @param labels the component labels
     * @param component the current component that is being built
     */
    private void search(DefinedMapping vertex, int currentLabel,
            int[] labels, List<DefinedMapping> component) {
        if (vertex.isVisited()) {
        } else {
            vertex.setVisited(true);
            labels[vertex.getIndex()] = currentLabel;
            component.add(vertex);
            adjacencyTable[vertex.getIndex()].stream().forEach((neighbour) -> {
                if (neighbour.isVisited()) {
                } else {
                    search(neighbour, currentLabel, labels, component);
                }
            });
        }
    }

    private class AtomContainerAtomPair {

        public IAtomContainer atomContainer = null;
        public IAtom atom = null;
        public int index = -1;

        AtomContainerAtomPair(IAtomContainer atomContainer, IAtom atom, int index) {
            this.atom = atom;
            this.atomContainer = atomContainer;
            this.index = index;
        }

        public int getIndex() {
            return this.index;
        }
    }
}
