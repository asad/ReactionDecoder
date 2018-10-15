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
import java.util.List;

/**
 * The list of blocks in an atom container.
 *
 * @author maclean
 *
 */
public class BlockList implements Comparable<BlockList> {

    private final List<Block> blocks;
    private String signatureString;
    private int totalSize;

    /**
     *
     */
    public BlockList() {
        blocks = new ArrayList<>();
        totalSize = 0;
    }

    /**
     *
     * @param block
     */
    public void add(Block block) {
        blocks.add(block);
    }

    /**
     *
     * @return
     */
    public int[] getBlockPermutation() {
        totalSize = 0;
        List<StringIntPair> signatureToBlockIndexMap;
        signatureToBlockIndexMap = new ArrayList<>();
        for (int i = 0; i < blocks.size(); i++) {
            Block block = blocks.get(i);
            String signature = block.getSignatureString();
            totalSize += block.getAtomCount();
            signatureToBlockIndexMap.add(new StringIntPair(signature, i));
        }

        // bucket sort
        sort(signatureToBlockIndexMap);

        int index = 0;
        int[] blockPermutation = new int[blocks.size()];
        for (StringIntPair stringIntPair : signatureToBlockIndexMap) {
            blockPermutation[index] = stringIntPair.intValue;
            index++;
        }
        return blockPermutation;
    }

    /**
     *
     * @return
     */
    public int[] getAtomContainerPermutation() {
        int[] blockPermutation = getBlockPermutation();

        int totalPermutationIndex = 0;
        int[] totalPermutation = new int[totalSize];
        int blockStart = 0;
        for (int blockIndex = 0; blockIndex < blocks.size(); blockIndex++) {
            int permutedBlockIndex = blockPermutation[blockIndex];
            Block block = blocks.get(permutedBlockIndex);
            SubgraphMoleculeSignature subgraphSignature
                    = block.getSubgraphSignature();
//            System.out.println("sig = " + subgraphSignature.toCanonicalString());
            int[] labels = subgraphSignature.getCanonicalLabels();
            for (int labelIndex = 0; labelIndex < labels.length; labelIndex++) {
                int x = blockStart + labels[labelIndex];
//                System.out.println("blockIndex = " + blockIndex
//                        + " blockStart = " + blockStart
//                        + " labelIndex = " + labelIndex
//                        + " labels = " + Arrays.toString(labels)
//                        + " totalPermutation = " + Arrays.toString(totalPermutation));
                totalPermutation[totalPermutationIndex] = x;
                totalPermutationIndex++;
            }
            blockStart += labels.length;

        }
        return totalPermutation;
    }

    private String calculateSignatureString() {
        List<String> signatures = new ArrayList<>();
        blocks.stream().forEach((block) -> {
            signatures.add(block.getSignatureString());
        });
        sort(signatures);
        StringBuilder sb = new StringBuilder();
        signatures.stream().forEach((s) -> {
            sb.append(s).append("|");
        });
        return sb.toString();
    }

    /**
     *
     * @return
     */
    public String getSignatureString() {
        if (signatureString == null) {
            signatureString = calculateSignatureString();
        }
        return signatureString;
    }

    /**
     *
     * @return
     */
    public int numberOfAtoms() {
        int total = 0;
        total = blocks.stream().map((block) -> block.getAtomCount()).reduce(total, Integer::sum);
        return total;
    }

    @Override
    public int compareTo(BlockList o) {
        int countThis = numberOfAtoms();
        int countOther = o.numberOfAtoms();
        if (countThis > countOther) {
            return -1;
        } else if (countThis < countOther) {
            return 1;
        } else {
            return getSignatureString().compareTo(o.getSignatureString());
        }
    }

    @Override
    public int hashCode() {
        return getSignatureString().hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof BlockList) {
            boolean signaturesEqual
                    = getSignatureString().equals(((BlockList) o).getSignatureString());
            if (signaturesEqual) {
                return ((BlockList) o).blocks.get(0).getAtomContainer()
                        == blocks.get(0).getAtomContainer();
            }
        }
        return false;
    }

    /**
     *
     * @param m
     * @return
     */
    public Block get(int m) {
        return blocks.get(m);
    }

    @Override
    public String toString() {
        return blocks.toString();
    }

    private class StringIntPair implements Comparable<StringIntPair> {

        public String stringValue;
        public int intValue;

        StringIntPair(String stringValue, int intValue) {
            this.stringValue = stringValue;
            this.intValue = intValue;
        }

        @Override
        public int compareTo(StringIntPair o) {
            return stringValue.compareTo(o.stringValue);
        }
    }
}
