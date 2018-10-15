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

import static java.lang.System.out;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.tools.labelling.AtomContainerPrinter;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalReactionLabeller;

/**
 *
 * @author asad
 */
public class BlockReactionCanoniser implements ICanonicalReactionLabeller {

    private IChemObjectBuilder builder;

    /**
     *
     */
    public BlockReactionCanoniser() {
        builder = null;
    }

    @Override
    public IReaction getCanonicalReaction(IReaction reaction) {
        BlockMapping blockMapping = new BlockMapping(reaction);

//        printReaction(reaction);
        builder = reaction.getBuilder();

        IReaction permutedReaction = builder.newInstance(IReaction.class);

        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet permutedReactants = permuteR(reactants, blockMapping);
        permutedReaction.setReactants(permutedReactants);

        IAtomContainerSet products = reaction.getProducts();
        IAtomContainerSet permutedProducts = permuteP(products, blockMapping);
        permutedReaction.setProducts(permutedProducts);

        replaceMappings(reaction, permutedReaction);
        permutedReaction.setID(reaction.getID());

        return permutedReaction;
    }

    private void replaceMappings(IReaction reaction, IReaction permutedReaction) {
        // replace the mappings, sorting them by the position of the first
        // atom (chemObject 0) in the mapping.
        final Map<IAtom, Integer> indexMap = new HashMap<>();
        int globalIndex = fillIndexMap(indexMap, permutedReaction.getReactants(), 0);
        fillIndexMap(indexMap, permutedReaction.getProducts(), globalIndex);

        Map<Integer, IMapping> mappingMap = new HashMap<>();
        List<IMapping> orphanMappings = new ArrayList<>();
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            if (a0 == null || a1 == null) {
                orphanMappings.add(mapping);
            } else {
                Integer id = indexMap.get(a0);
                if (id == null) {
//                    System.out.println("null for " + a0.getID() + " " + a1.getID());
                    orphanMappings.add(mapping);
                    continue;
                }
                mappingMap.put(id, mapping);
            }
        }
        List<Integer> keys = new ArrayList<>(mappingMap.keySet());

//        System.out.println("presorted keys = " + keys);
        sort(keys);
//        System.out.println("sorted keys = " + keys);
        for (Integer key : keys) {
            permutedReaction.addMapping(mappingMap.get(key));
        }

        // add any 'orphan' mappings that map atoms that are not referenced in
        // the atom containers of the reaction...
        orphanMappings.stream().forEach((mapping) -> {
            permutedReaction.addMapping(mapping);
        });//        assert reaction.getMappingCount() == permutedReaction.getMappingCount();
//        System.out.println(reaction.getMappingCount() + " " + permutedReaction.getMappingCount());
    }

    private int fillIndexMap(
            Map<IAtom, Integer> indexMap, IAtomContainerSet moleculeSet, int index) {
        for (IAtomContainer atomContainer : moleculeSet.atomContainers()) {
            for (IAtom atom : atomContainer.atoms()) {
                indexMap.put(atom, index);
                index++;
            }
        }
        return index;
    }

    private IAtomContainerSet permuteR(IAtomContainerSet original, BlockMapping mapping) {
        int[] reactantPermutation = mapping.getPermutationOfReactants();
        IAtomContainerSet permutedContainers = builder.newInstance(IAtomContainerSet.class);
        // System.out.println("REACTANTS --------------------");

        List<IAtomContainer> unusedContainers = new ArrayList<>();
        for (IAtomContainer ac : original.atomContainers()) {
            unusedContainers.add(ac);
        }
        for (int i = 0; i < original.getAtomContainerCount(); i++) {
            int pi = reactantPermutation[i];
            IAtomContainer container = original.getAtomContainer(pi);
            BlockList blockList = mapping.getBlockListForReactant(container);

            if (blockList == null) {
                out.println("blocklist null for " + i
                        + new AtomContainerPrinter().toString(container)
                        + " in " + java.util.Arrays.toString(reactantPermutation));
                continue;
            }

            unusedContainers.remove(container);

            IAtomContainer permutedContainer = builder.newInstance(IAtomContainer.class);
            permutedContainer.setID(container.getID());
            // System.out.println("AC " + i + " " + blockList);

            IAtom[] atoms = getPermutedAtomsForReactant(blockList, container);
            // System.out.println("AC " + i + " " + blockList);
            permutedContainer.setAtoms(atoms);
            for (IBond bond : container.bonds()) {
                permutedContainer.addBond(bond);
            }
            permutedContainers.addAtomContainer(permutedContainer);
        }
        unusedContainers.stream().forEach((unusedContainer) -> {
            permutedContainers.addAtomContainer(unusedContainer);
        });
        return permutedContainers;
    }

    private IAtomContainerSet permuteP(IAtomContainerSet original, BlockMapping mapping) {
        IAtomContainerSet permutedContainers = builder.newInstance(IAtomContainerSet.class);
        int[] productPermutation = mapping.getPermutationOfProducts();
        // System.out.println("PRODUCTS --------------------");

        List<IAtomContainer> unusedContainers = new ArrayList<>();
        for (IAtomContainer ac : original.atomContainers()) {
            unusedContainers.add(ac);
        }
        for (int i = 0; i < original.getAtomContainerCount(); i++) {
            int pi = productPermutation[i];
            IAtomContainer container = original.getAtomContainer(pi);
            BlockList blockList = mapping.getBlockListForProduct(container);
            if (blockList == null) {
                out.println("blocklist null for " + i
                        + new AtomContainerPrinter().toString(container)
                        + " in " + java.util.Arrays.toString(productPermutation));
                continue;
            }
            unusedContainers.remove(container);

            IAtomContainer permutedContainer = builder.newInstance(IAtomContainer.class);
            permutedContainer.setID(container.getID());
            // System.out.println("AC " + i + " " + blockList);
            IAtom[] atoms = getPermutedAtomsForProduct(blockList, container);
            permutedContainer.setAtoms(atoms);
            for (IBond bond : container.bonds()) {
                permutedContainer.addBond(bond);
            }
            permutedContainers.addAtomContainer(permutedContainer);
        }
        unusedContainers.stream().forEach((unusedContainer) -> {
            permutedContainers.addAtomContainer(unusedContainer);
        });
        return permutedContainers;
    }

    private void bucketSort(List<IAtom> atoms, IAtomContainer ac) {
        final Map<IAtom, Integer> indexMap = new HashMap<>();

        atoms.stream().forEach((atom) -> {
            indexMap.put(atom, ac.indexOf(atom));
        });
        Comparator<IAtom> sorter = (IAtom o1, IAtom o2) -> indexMap.get(o1).compareTo(indexMap.get(o2));
        sort(atoms, sorter);
    }

    private IAtom[] getPermutedAtomsForReactant(BlockList blockList,
            IAtomContainer container) {
        int[] blockPermutation = blockList.getBlockPermutation();

        int m = blockPermutation.length;
        IAtom[] atoms = new IAtom[container.getAtomCount()];

        // keeps track of atoms added to the permutation
        List<IAtom> unusedAtoms = new ArrayList<>();

        for (IAtom atom : container.atoms()) {
            unusedAtoms.add(atom);
        }

        int blockStart = 0;
        for (int blockIndex = 0; blockIndex < m; blockIndex++) {
            Block block = blockList.get(blockPermutation[blockIndex]);

            int[] labels = block.getLabels();

            List<IAtom> subgraphAtoms = block.getAtoms();
            bucketSort(subgraphAtoms, container);

            int indexInBlock = 0;

            for (IAtom atom : subgraphAtoms) {
                int newIndex = labels[indexInBlock];
                atoms[blockStart + newIndex] = atom;
                unusedAtoms.remove(atom);
                indexInBlock++;
            }

            blockStart += labels.length;
        }

        int numberUnused = unusedAtoms.size();
        if (numberUnused != 0) {
            int index = atoms.length - numberUnused;
            for (IAtom unusedAtom : unusedAtoms) {
                atoms[index] = unusedAtom;
                index++;
            }
        }
        return atoms;
    }

    private IAtom[] getPermutedAtomsForProduct(BlockList blockList,
            IAtomContainer container) {
        int[] blockPermutation = blockList.getBlockPermutation();

        int m = blockPermutation.length;
        IAtom[] atoms = new IAtom[container.getAtomCount()];

        // keeps track of atoms added to the permutation
        List<IAtom> unusedAtoms = new ArrayList<>();

        for (IAtom atom : container.atoms()) {
            unusedAtoms.add(atom);
        }

        int blockStart = 0;

        for (int blockIndex = 0; blockIndex < m; blockIndex++) {
            Block block = blockList.get(blockPermutation[blockIndex]);

            int[] labels = block.getPartner().getLabels();

            List<IAtom> subgraphAtoms = block.getAtoms();
            bucketSort(subgraphAtoms, container);

            int[] mappingPermutation = block.getMappingPermutation();

            int indexInBlock = 0;
            printAtomIndices(subgraphAtoms, container);

            for (IAtom atom : subgraphAtoms) {
                int newIndex = labels[mappingPermutation[indexInBlock]];
                atoms[blockStart + newIndex] = atom;
                unusedAtoms.remove(atom);
                indexInBlock++;
            }

            blockStart += labels.length;
        }

        int numberUnused = unusedAtoms.size();

        if (numberUnused != 0) {
            int index = atoms.length - numberUnused;

            for (IAtom unusedAtom : unusedAtoms) {
                atoms[index] = unusedAtom;
                index++;
            }
        }
        return atoms;
    }

    private void printReaction(IReaction reaction) {
        for (IAtomContainer atomContainer : reaction.getReactants().atomContainers()) {
            for (IAtom atom : atomContainer.atoms()) {
                out.println("PRINTING ATOM " + atom.getID());
            }
        }
        for (IAtomContainer atomContainer : reaction.getProducts().atomContainers()) {
            for (IAtom atom : atomContainer.atoms()) {
                out.println("PRINTING ATOM " + atom.getID());
            }
        }
    }

    private void printAtomIndices(List<IAtom> atoms, IAtomContainer container) {
        List<Integer> indices = new ArrayList<>();

        for (int i = 0; i < atoms.size(); i++) {
            IAtom atom = atoms.get(i);

            if (atom == null) {
                // System.out.println("atom " + i + " is null");
            }
            int index = container.indexOf(atom);
            indices.add(index);

        } // System.out.println("list atoms " + indices);
    }

    private void printAtomIndices(IAtom[] atoms, IAtomContainer container) {
        List<Integer> indices = new ArrayList<>();

        for (int i = 0; i < container.getAtomCount(); i++) {
            IAtom atom = atoms[i];

            if (atom == null) {
                out.println("atom " + i + " is null");
            }
            int index = container.indexOf(atom);
            indices.add(index);

        }
        out.println("array atoms " + indices);

    }

    private String tmpAtomPrint(Iterable<IAtom> atoms, IAtomContainer ac) {
        String s = "";

        for (IAtom atom : atoms) {
            if (atom == null) {
                s += "!";
            } else {
                s += atom.getSymbol() + ac.indexOf(atom) + "("
                        + atom.getID() + ")";
            }
        }
        return s;

    }
}
