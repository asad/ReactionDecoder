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

import static java.lang.System.getProperty;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAtomCount;

/**
 * Converts the flat list of atom-atom mappings into blocks of mappings such
 * that two mappings [(r1, p1), (r2, p2)] are in the same block only if r1 and
 * r2 are in the same reactant atom container and p1 and p2 are in the same
 * product atom container AND there is a bond between r1 and r2 and a bond
 * between p1 and p2.
 *
 * @author maclean
 *
 */
public class BlockMapping {

    static final String NEW_LINE = getProperty("line.separator");

    private final IReaction reaction;

    private final List<Block> reactantBlocks;

    private final List<Block> productBlocks;

    private final Map<IAtomContainer, BlockList> reactantBlockMap;

    private final Map<IAtomContainer, BlockList> productBlockMap;

    private final List<BlockPair> blockPairs;

    /**
     *
     * @param reaction
     */
    public BlockMapping(IReaction reaction) {
        this.reaction = reaction;
        reactantBlocks = new ArrayList<>();
        productBlocks = new ArrayList<>();

        reactantBlockMap = new HashMap<>();
        productBlockMap = new HashMap<>();

        MappingGraph mappingGraph = new MappingGraph(reaction);
        blockPairs = mappingGraph.createBlockPairs(reaction);
        blockPairs.stream().forEach((pair) -> {
            Block reactantBlock = pair.getReactantBlock();
            Block productBlock = pair.getProductBlock();
            reactantBlocks.add(reactantBlock);
            productBlocks.add(productBlock);
            addBlockToAtomContainerMap(reactantBlock, reactantBlockMap);
            addBlockToAtomContainerMap(productBlock, productBlockMap);
        });
    }

    /**
     *
     * @return
     */
    public List<BlockPair> getBlockPairs() {
        return blockPairs;
    }

    /**
     *
     * @return
     */
    public List<Block> getReactantBlocks() {
        return reactantBlocks;
    }

    /**
     *
     * @return
     */
    public List<Block> getProductBlocks() {
        return productBlocks;
    }

    /**
     *
     * @return
     */
    public IReaction getReaction() {
        return this.reaction;
    }

    private void addBlockToAtomContainerMap(
            Block block, Map<IAtomContainer, BlockList> map) {
        BlockList blocksForContainer;
        IAtomContainer key = block.getAtomContainer();
        if (map.containsKey(key)) {
            blocksForContainer = map.get(key);
        } else {
            blocksForContainer = new BlockList();
            map.put(key, blocksForContainer);
        }
        blocksForContainer.add(block);
    }

    /**
     *
     * @return
     */
    public int[] getTotalPermutation() {
        int n = getAtomCount(reaction);
        int[] totalPermutation = new int[n];
        int totalIndex = 0;
        int offset = 0;

        int[] reactantPermutation = getPermutationOfReactants();
        IAtomContainerSet reactants = reaction.getReactants();
        for (int i = 0; i < reaction.getReactantCount(); i++) {
            int permutedReactantIndex = reactantPermutation[i];
            IAtomContainer reactant
                    = reactants.getAtomContainer(permutedReactantIndex);
            BlockList blockList = reactantBlockMap.get(reactant);
            int[] atomContainerPermutation
                    = blockList.getAtomContainerPermutation();
            for (int j = 0; j < atomContainerPermutation.length; j++) {
                int x = atomContainerPermutation[j];
                totalPermutation[totalIndex] = offset + x;
                totalIndex++;
            }
            offset += atomContainerPermutation.length;
        }

        int[] productPermutation = getPermutationOfProducts();
        IAtomContainerSet products = reaction.getProducts();
        for (int i = 0; i < reaction.getProductCount(); i++) {
            int permutedProductIndex = productPermutation[i];
            IAtomContainer product
                    = products.getAtomContainer(permutedProductIndex);
            BlockList blockList = productBlockMap.get(product);
            int[] atomContainerPermutation
                    = blockList.getAtomContainerPermutation();
            for (int j = 0; j < atomContainerPermutation.length; j++) {
                int x = atomContainerPermutation[j];
                totalPermutation[totalIndex] = offset + x;
                totalIndex++;
            }
            offset += atomContainerPermutation.length;
        }

        return totalPermutation;
    }

    /**
     * Gets a permutation of the reactant atom containers.
     *
     * @return
     */
    public int[] getPermutationOfReactants() {
        return getContainerPermutation(reaction.getReactants(), reactantBlockMap);
    }

    /**
     *
     * @param reactant
     * @return
     */
    public int[] getContainerPermutationForReactant(IAtomContainer reactant) {
        return reactantBlockMap.get(reactant).getAtomContainerPermutation();
    }

    /**
     * Gets a permutation of the product atom containers.
     *
     * @return
     */
    public int[] getPermutationOfProducts() {
        return getContainerPermutation(reaction.getProducts(), productBlockMap);
    }

    /**
     *
     * @param product
     * @return
     */
    public int[] getContainerPermutationForProduct(IAtomContainer product) {
        return productBlockMap.get(product).getAtomContainerPermutation();
    }

    private int[] getContainerPermutation(
            IAtomContainerSet containers,
            Map<IAtomContainer, BlockList> map) {
        int n = containers.getAtomContainerCount();

        // map the block lists for each atom container to the original index
        Map<BlockList, Integer> blockListToOriginalOrderMap
                = new HashMap<>();
        for (int i = 0; i < n; i++) {
            IAtomContainer container = containers.getAtomContainer(i);
            BlockList blockList = map.get(container);
            if (blockList != null) {
                blockListToOriginalOrderMap.put(blockList, i);
            }
        }

        // bucket-sort the original indices of the containers
        List<BlockList> blockListKeys
                = new ArrayList<>(blockListToOriginalOrderMap.keySet());
//        System.out.println(blockListKeys);
        sort(blockListKeys);

        // construct a permutation of the atom containers
        int[] containerPermutation = new int[n];
        int index = 0;
        for (BlockList blockList : blockListKeys) {
            int originalOrder = blockListToOriginalOrderMap.get(blockList);
            containerPermutation[index] = originalOrder;
            index++;
        }
        return containerPermutation;
    }

    /**
     *
     * @param container
     * @return
     */
    public BlockList getBlockListForReactant(IAtomContainer container) {
        return reactantBlockMap.get(container);
    }

    /**
     *
     * @param container
     * @return
     */
    public BlockList getBlockListForProduct(IAtomContainer container) {
        return productBlockMap.get(container);
    }

    @Override
    public String toString() {
        String rbm = "{";
        rbm = reactantBlockMap.keySet().stream().map((r) -> r.getID() + ":" + reactantBlockMap.get(r) + NEW_LINE).reduce(rbm, String::concat);
        rbm += "}" + NEW_LINE;

        String pbm = "{";
        pbm = productBlockMap.keySet().stream().map((p) -> p.getID() + ":" + productBlockMap.get(p) + NEW_LINE).reduce(pbm, String::concat);
        pbm += "}" + NEW_LINE;

        return reactantBlocks.toString() + NEW_LINE
                + productBlocks.toString() + NEW_LINE
                + rbm
                + pbm;
    }
}
