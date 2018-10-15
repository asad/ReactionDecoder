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
package uk.ac.ebi.reactionblast.fingerprints.tools;

import static java.lang.Math.sqrt;
import java.util.BitSet;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 *
 * <PRE>
 *
 *   The Measures program takes as input any fixed length bit strings,
 *   these can be from the Mesa Fingerprint programs or user supplied fingerprints.
 *   User supplied fingerprints must take the form of ASCII 1's and 0's, (e.g. 011100001111000....),
 *   ASCII CDK fingerprints inside the FP<> data type are also valid input to  Measures .
 *   The Measures program  produces a similarity or dissimilarity matrix (user's choice)
 *   using one of the following user selected measures:
 *   <B> Tversky, Tanimoto, Euclidean, Hamman, or Ochia (1-Cosine).
 *
 * In similarity form:
 *
 *                            Tanimoto(bitset1,bitset2)  = c / [a + b - c]  (symmetric)
 *
 *                            Euclidean(bitset1,bitset2) = 1 - {[(a + b)] / n}(1/2)   (symmetric)
 *
 *                            Hamman(bitset1,bitset2)  = [c + d] /n  (symmetric)
 *
 *                            Ochia(bitset1,bitset2) = 1 - Cosine(bitset1,bitset2) = c / [(c + a) * (c + b)](1/2)  (symmetric)
 *
 *                            Tversky(bitset1,bitset2) = c / [(alpha) * a + (beta) * b + c]  (asymmetric)
 *
 *                            a : Unique bits turned on in molecule "bitset1"
 *                            b:  Unique bits turned on in molecule "bitset2"
 *                            c:  Common bits turned on in both molecule "bitset1" and molecule "bitset2"
 *                            d:  Common bits turned off in both molecule "bitset1" and molecule "bitset2"
 *                            n:  The total number of bits in the fingerprint
 * </B>
 *
 *  <B> Note:</B>The Tanimoto, Euclidean, Hamman, and Ochai are all symmetric measures.
 *
 *  <U> This means that the comparison of bitset1 to bitset2 yields the same number as the comparison of compound bitset2 to compound bitset1.</U>
 *  <B> Note:</B> The dissimilarity is just 1 - similarity.
 *
 *
 *
 * </PRE>
 *
 * @ref <B>Willett et.al., Chemical Similarity Searching,</B> <I>J.Chem. Inf.
 * Comput. Sci.</I>, Vol. 38, No. 6, 1998
 *
 *
 */
public class Similarity {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Similarity.class);

    /**
     *
     * @param Molecule1 BitSet
     * @param Molecule2 BitSet
     * @return <B>Similarity <U>Tanimoto, Jaccard</U> </B>
     * <B>c/(a+b-c)></B>
     * @throws java.lang.Exception
     */
    public static synchronized float getTanimotoSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
        BitSet bitset1 = (BitSet) Molecule1.clone();
        BitSet bitset2 = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = bitset1.cardinality();
        float _bitset2_cardinality = bitset2.cardinality();

//        System.out.println("bitset1: "+ bitset1.size() + " " + " bitset2" + bitset2.size());
        if (bitset1.size() != bitset2.size()) {
            throw new Exception("BitSets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) bitset1.clone();
        one_and_two.and(bitset2);
        float _common_bit_count = one_and_two.cardinality();
        return _common_bit_count / (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     * @return <B>Similarity <U>Cosine,Ochiai,Carbo</U></B>
     * <B>c/sqrt(a*b)</B>
     * @throws Exception
     */
    public static synchronized double getCosineSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
        BitSet bitset1 = (BitSet) Molecule1.clone();
        BitSet bitset2 = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = bitset1.cardinality();
        float _bitset2_cardinality = bitset2.cardinality();

        if (bitset1.size() != bitset2.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) bitset1.clone();
        one_and_two.and(bitset2);
        float _common_bit_count = one_and_two.cardinality();

        return _common_bit_count / (sqrt(_bitset1_cardinality * _bitset2_cardinality));
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     * @return <B>Similarity <U>Dice, Sorensen, Czekanowski,
     * Hodgkin-Richards</U></B>
     * <B>2c/(a+b)</B>
     * @throws Exception
     *
     */
    public static synchronized double getDiceSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
        BitSet bitset1 = (BitSet) Molecule1.clone();
        BitSet bitset2 = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = bitset1.cardinality();
        float _bitset2_cardinality = bitset2.cardinality();

        if (bitset1.size() != bitset2.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) bitset1.clone();
        one_and_two.and(bitset2);
        float _common_bit_count = one_and_two.cardinality();

        return 2 * _common_bit_count / (_bitset1_cardinality + _bitset2_cardinality);
    }

    private Similarity() {
    }
}
