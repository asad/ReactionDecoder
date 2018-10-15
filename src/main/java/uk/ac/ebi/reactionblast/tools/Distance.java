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
package uk.ac.ebi.reactionblast.tools;

import static java.lang.Math.sqrt;
import java.util.BitSet;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
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
 *                            Tanimoto(A,B)  = c / [a + b - c]  (symmetric)
 *
 *                            Euclidean(A,B) = 1 - {[(a + b)] / n}(1/2)   (symmetric)
 *
 *                            Hamman(A,B)  = [c + d] /n  (symmetric)
 *
 *                            Ochia(A,B) = 1 - Cosine(A,B) = c / [(c + a) * (c + b)](1/2)  (symmetric)
 *
 *                            Tversky(A,B) = c / [(alpha) * a + (beta) * b + c]  (asymmetric)
 *
 *                            a : Unique bits turned on in molecule "A"
 *                            b:  Unique bits turned on in molecule "B"
 *                            c:  Common bits turned on in both molecule "A" and molecule "B"
 *                            d:  Common bits turned off in both molecule "A" and molecule "B"
 *                            n:  The total number of bits in the fingerprint
 * </B>
 *
 *  <B> Note:</B>The Tanimoto, Euclidean, Hamman, and Ochai are all symmetric measures.
 *
 *  <U> This means that the comparison of A to B yields the same number as the comparison of compound B to compound A.</U>
 *  <B> Note:</B> The dissimilarity is just 1 - similarity.
 *
 *
 *
 * </PRE>
 *
 * @ref <B>Willett et.al., Chemical Similarity Searching,</B> <I>J.Chem. Inf.
 * Comput. Sci.</I>, Vol. 38, No. 6, 1998
 *
 */
public class Distance {

    private static BitSet A;
    private static BitSet B;

    /**
     *
     * @param Molecule1 BitSet
     * @param Molecule2 BitSet
     * <B>Note:</B> Can be Normalized to the range 1 to 0 if all attributes are
     * normalized to this range and the result is divided in N
     * @return <B>Distance (Dissimilarity) <U>Euclidean</U></B>
     * <B>sqrt(a+b-2c)</B>
     * <B>sqrt(OnA+OnB-2*c)</B>
     * @throws java.lang.Exception
     */
    public static double getEuclideanDistance(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        double Val = sqrt((_bitset1_cardinality + _bitset2_cardinality - _common_bit_count) - _common_bit_count);

        return Val;
    }

    /**
     *
     * @param Molecule1 BitSet
     * @param Molecule2 BitSet
     * <B>Note:</B> Can be Normalized to the range 1 to 0 if all attributes are
     * normalized to this range and the result is divided in N
     * @return <B>Distance (Dissimilarity) <U>Euclidean</U></B>
     * <B>sqrt(a+b-2c)</B>
     * <B>sqrt(OnA+OnB-2*c)</B>
     * @throws java.lang.Exception
     */
    public static double getNormalizedEuclideanDistance(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        double Val = sqrt((_bitset1_cardinality + _bitset2_cardinality - _common_bit_count) - _common_bit_count);

        float N1 = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
        float N2 = (1 - _bitset1_cardinality - _bitset2_cardinality + _common_bit_count);

        float N = N1 + N2;
        return Val / N;
    }

    /*
     * <B>Note:</B> Can be Normalized to the range 1 to 0 if all 
     * attributes are normalized to this range and the result is divided in N
     * @param Molecule1 
     * @param Molecule2 
     * @return <B>Distance (Dissimilarity)<U> Hamming</U> </B>
     * <B>a+b-2c</B>
     * @throws Exception 
     */
    /**
     *
     * @param Molecule1
     * @param Molecule2
     * @return
     * @throws Exception
     */
    public static double getHammingDistance(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        double Val = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count) - _common_bit_count;

//        float N1 = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
//        float N2 = (1 - _bitset1_cardinality - _bitset2_cardinality + _common_bit_count);
//
//        float N = N1 + N2;
        return Val;
    }

    /*
     * <B>Note:</B> Can be Normalized to the range 1 to 0 if all 
     * attributes are normalized to this range and the result is divided in N
     * @param Molecule1 
     * @param Molecule2 
     * @return <B>Distance (Dissimilarity)<U> Hamming</U> </B>
     * <B>a+b-2c</B>
     * @throws Exception 
     */
    /**
     *
     * @param Molecule1
     * @param Molecule2
     * @return
     * @throws Exception
     */
    public static double getNormalizedHammingDistance(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        double Val = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count) - _common_bit_count;

        float N1 = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
        float N2 = (1 - _bitset1_cardinality - _bitset2_cardinality + _common_bit_count);

        float N = N1 + N2;
        return Val / N;
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     * @return <B>Distance (Dissimilarity)<U> Soergel</U></B>
     * <B>1-c/[a+b-c]=[a+b-2*c]/[a+b-c]</B>
     * @throws Exception
     */
    public static double getSoergelDistance(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        double Val = (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count) / _common_bit_count;

        return Val;

    }

    private Distance() {
    }
}
