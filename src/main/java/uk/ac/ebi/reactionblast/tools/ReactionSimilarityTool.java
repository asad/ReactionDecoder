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

import static java.lang.Double.parseDouble;
import static java.lang.Math.log;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import static uk.ac.ebi.reactionblast.tools.utility.EBIDoubleUtility.append;

/**
 * This tool finds reaction similarity and distance based on our in-house
 * scoring functions The similarity score can be further weight by alpha and
 * beta factor coding for bond change and structure weight respectively
 *
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionSimilarityTool {

    /**
     *
     * @param alpha
     * @param beta
     * @param bondSimilarity
     * @param structuralSimilarity
     * @return
     * @throws CDKException
     */
    public static double getSimilarityScore(double alpha,
            double beta,
            double bondSimilarity,
            double structuralSimilarity) throws CDKException {

        double score = 0;
        // In maths S=e pow(x) * e pow(y) = e pow(x+y)
        double BC = bondSimilarity;
        double SC = structuralSimilarity;

        if (alpha == 0 && beta == 0) {
            throw new CDKException("Both alpha & beta can't be zero at the same time");
        } else if (BC == 0 && SC == 0) {
            score = 0;
        } else {
            double structureScore = (beta / (alpha + beta)) * SC;
            double bondScore = (alpha / (alpha + beta)) * BC;
            score = bondScore + structureScore;
        }
//        System.out.println("alpha: " + alpha + "\tbeta: " + beta + "\tBondSimilarity: " + BC + "\tStructuralSimilarity: " + SC + "\tScore: " + score);

        DecimalFormat df = new DecimalFormat("0.00");
        df.setMaximumFractionDigits(2);
        String a = df.format(score);
        score = parseDouble(a);
        return score;
    }

    /**
     *
     * @param alpha
     * @param beta
     * @param gamma
     * @param bondSimilarity
     * @param structuralSimilarity
     * @param stereoSimilarity
     * @return
     * @throws CDKException
     */
    public static double getSimilarityScore(
            double alpha,
            double beta,
            double gamma,
            double bondSimilarity,
            double structuralSimilarity,
            double stereoSimilarity) throws CDKException {

        double score = 0;
        // In maths S=e pow(x) * e pow(y) = e pow(x+y)
        double BC = bondSimilarity;
        double SC = structuralSimilarity;
        double ST = stereoSimilarity;

        if (alpha == 0 && beta == 0 && gamma == 0) {
            throw new CDKException("Both alpha & beta can't be zero at the same time");
        } else if (BC == 0 && SC == 0) {
            score = 0;
        } else {
            double structureScore = (beta / (alpha + beta + gamma)) * SC;
            double bondScore = (alpha / (alpha + beta + gamma)) * BC;
            double stereoScore = (gamma / (alpha + beta + gamma)) * ST;
            score = bondScore + structureScore + stereoScore;
        }
//        System.out.println("alpha: " + alpha + "\tbeta: " + beta + "\tBondSimilarity: " + BC + "\tStructuralSimilarity: " + SC + "\tScore: " + score);

        DecimalFormat df = new DecimalFormat("0.00");
        df.setMaximumFractionDigits(2);
        String a = df.format(score);
        score = parseDouble(a);
        return score;
    }

    /**
     *
     * @param alpha weight for BondFP
     * @param beta weight for StructuralFP
     * @param StructFP1
     * @param bondFeatures1
     * @param StructFP2
     * @param bondFeatures2
     * @param oxidationPenalty Penalty of 25% on the bond changes if they don't
     * have same similar oxidation type (oxidation and non-oxidative reactions
     * will face penalty) @Equation (alpha/(alpha+beta)) * BondChange +
     * (beta/alpha+beta)* StructuralChange
     * @return Gives overall similarity between two reactions (Higher the score,
     * farthest the reaction is)
     * @throws Exception
     */
    public static double getReactionSimilarity(
            double alpha,
            double beta,
            IPatternFingerprinter StructFP1,
            double[] bondFeatures1,
            IPatternFingerprinter StructFP2,
            double[] bondFeatures2,
            boolean oxidationPenalty)
            throws Exception {

        double score = 0.0;

        int size1 = StructFP1.getFingerprintSize();
        int size2 = StructFP2.getFingerprintSize();

        if (size1 != size2) {
            throw new CDKException("Features vectors must be of the same length");

        } else if (alpha < 0. && beta < 0.) {
            throw new CDKException("both alpha and beta can't =< be zero");
        } else {

            double similarityOfBondChanges = getSimilarity(bondFeatures1, bondFeatures2);

            double[] structFeatures1 = StructFP1.getWeightedHashedFingerPrint();
            double[] structFeatures2 = StructFP2.getWeightedHashedFingerPrint();
            double similarityStruct = getSimilarity(structFeatures1, structFeatures2);
//
//            System.out.println("Alpha " + alpha + " Beta " + beta
//                    + " mechanism similarity Score: " + similarityOfBondChanges
//                    + " structure similarity Score: " + similarityStruct);

            score = getSimilarityScore(alpha, beta, similarityOfBondChanges, similarityStruct);
        }

        return score;
    }

    /**
     *
     * @param alpha weight for BondFP
     * @param beta weight for StructuralFP
     * @param StructFP1
     * @param BondsCF1 Reactant Bond Cleaved/Formed FP for Query Reaction
     * @param BondsCF2 Product Bond Cleaved/Formed FP for Target Reaction
     * @param BondsOC1 Reactant Bond Order/Change FP for Query Reaction
     * @param StructFP2
     * @param BondsOC2 Product Bond Order/Change FP for Target Reaction
     * @param BondsST1 Reactant Bond Stereo FP for Query Reaction
     * @param BondsST2 Product Bond Stereo FP for Query Reaction
     *
     * @Equation (alpha/(alpha+beta)) * BondChange + (beta/alpha+beta)*
     * StructuralChange
     * @return Gives overall similarity between two reactions (Higher the score,
     * farthest the reaction is)
     * @throws Exception
     */
    public static double getReactionSimilarity(
            double alpha,
            double beta,
            IPatternFingerprinter StructFP1,
            IPatternFingerprinter BondsCF1,
            IPatternFingerprinter BondsOC1,
            IPatternFingerprinter BondsST1,
            IPatternFingerprinter StructFP2,
            IPatternFingerprinter BondsCF2,
            IPatternFingerprinter BondsOC2,
            IPatternFingerprinter BondsST2)
            throws Exception {
        double score = 0.0;
        if (alpha < 0. && beta < 0.) {
            throw new CDKException("both alpha and beta can't =< be zero");
        } else {

            double[] bondFeatures1 = append(BondsCF1.getWeightedHashedFingerPrint(),
                    BondsOC1.getWeightedHashedFingerPrint(), BondsST1.getWeightedHashedFingerPrint());
            double[] bondFeatures2 = append(BondsCF2.getWeightedHashedFingerPrint(),
                    BondsOC2.getWeightedHashedFingerPrint(), BondsST2.getWeightedHashedFingerPrint());
            double similarityOfBondChanges = getSimilarity(bondFeatures1, bondFeatures2);

            double[] structFeatures1 = StructFP1.getWeightedHashedFingerPrint();
            double[] structFeatures2 = StructFP2.getWeightedHashedFingerPrint();
            double similarityStruct = getSimilarity(structFeatures1, structFeatures2);
//
//            System.out.println("Alpha " + alpha + " Beta " + beta
//                    + " mechanism similarity Score: " + similarityOfBondChanges
//                    + " structure similarity Score: " + similarityStruct);

            score = getSimilarityScore(alpha, beta, similarityOfBondChanges, similarityStruct);
        }

        return score;
    }

    /**
     *
     * @param alpha weight for BondFP
     * @param beta weight for StructuralFP
     * @param gamma weight for Bond Stereo FP
     * @param StructFP1
     * @param BondsCF1 Reactant Bond Cleaved/Formed FP for Query Reaction
     * @param BondsCF2 Product Bond Cleaved/Formed FP for Target Reaction
     * @param BondsOC1 Reactant Bond Order/Change FP for Query Reaction
     * @param StructFP2
     * @param BondsOC2 Product Bond Order/Change FP for Target Reaction
     * @param BondsST1 Reactant Bond Stereo FP for Query Reaction
     * @param BondsST2 Product Bond Stereo FP for Query Reaction
     *
     * @Equation (alpha/(alpha+beta)) * BondChange + (beta/alpha+beta)*
     * StructuralChange
     * @return Gives overall similarity between two reactions (Higher the score,
     * farthest the reaction is)
     * @throws Exception
     */
    public static double getReactionSimilarity(
            double alpha,
            double beta,
            double gamma,
            IPatternFingerprinter StructFP1,
            IPatternFingerprinter BondsCF1,
            IPatternFingerprinter BondsOC1,
            IPatternFingerprinter BondsST1,
            IPatternFingerprinter StructFP2,
            IPatternFingerprinter BondsCF2,
            IPatternFingerprinter BondsOC2,
            IPatternFingerprinter BondsST2)
            throws Exception {

        double score = 0.0;
        if (alpha < 0. && beta < 0.) {
            throw new CDKException("both alpha and beta can't =< be zero");
        } else {

            double[] bondFeatures1 = append(BondsCF1.getWeightedHashedFingerPrint(),
                    BondsOC1.getWeightedHashedFingerPrint());
            double[] bondFeatures2 = append(BondsCF2.getWeightedHashedFingerPrint(),
                    BondsOC2.getWeightedHashedFingerPrint());
            double similarityOfBondChanges = getSimilarity(bondFeatures1, bondFeatures2);
            double similarityOfStereoChanges = getSimilarity(BondsST1.getWeightedHashedFingerPrint(), BondsST2.getWeightedHashedFingerPrint());

            double[] structFeatures1 = StructFP1.getWeightedHashedFingerPrint();
            double[] structFeatures2 = StructFP2.getWeightedHashedFingerPrint();
            double similarityStruct = getSimilarity(structFeatures1, structFeatures2);
//
//            System.out.println("Alpha " + alpha + " Beta " + beta
//                    + " mechanism similarity Score: " + similarityOfBondChanges
//                    + " structure similarity Score: " + similarityStruct);

            score = getSimilarityScore(alpha, beta, gamma, similarityOfBondChanges, similarityStruct, similarityOfStereoChanges);
        }

        return score;
    }

    /**
     *
     * @param alpha weight for BondFP
     * @param beta weight for StructuralFP
     * @param StructFP1
     * @param BondsCF1 Reactant Bond Cleaved/Formed FP for Query Reaction
     * @param BondsCF2 Product Bond Cleaved/Formed FP for Target Reaction
     * @param BondsOC1 Reactant Bond Order/Change FP for Query Reaction
     * @param StructFP2
     * @param BondsOC2 Product Bond Order/Change FP for Target Reaction
     * @param BondsST1 Reactant Bond Stereo FP for Query Reaction
     * @param BondsST2 Product Bond Stereo FP for Query Reaction
     *
     * @Equation (alpha/(alpha+beta)) * BondChange + (beta/alpha+beta)*
     * StructuralChange
     * @return Gives overall similarity between two reactions (Higher the score,
     * farthest the reaction is)
     * @throws CDKException
     * @throws Exception
     */
    public static double getPointWiseMutualInformation(double alpha, double beta, IPatternFingerprinter StructFP1,
            IPatternFingerprinter BondsCF1, IPatternFingerprinter BondsOC1, IPatternFingerprinter BondsST1, IPatternFingerprinter StructFP2,
            IPatternFingerprinter BondsCF2, IPatternFingerprinter BondsOC2, IPatternFingerprinter BondsST2)
            throws CDKException, Exception {

//    System.out.println("Calling weighted score");
        double score = 0.0;

        int size1 = StructFP1.getFingerprintSize();
        int size2 = StructFP2.getFingerprintSize();

        if (size1 != size2) {
            throw new CDKException("Features vectors must be of the same length");

        } else if (alpha < 0. && beta < 0.) {
            throw new CDKException("both alpha and beta can't =< be zero");
        } else {

//            double similarityOfBondChanges = new Float(getBondChangeDistance(BondsCF1, BondsCF2, BondsOC1, BondsOC2)).floatValue();
//
            double[] bondFeatures1 = append(BondsCF1.getWeightedHashedFingerPrint(),
                    BondsOC1.getWeightedHashedFingerPrint(), BondsST1.getWeightedHashedFingerPrint());
            double[] bondFeatures2 = append(BondsCF2.getWeightedHashedFingerPrint(),
                    BondsOC2.getWeightedHashedFingerPrint(), BondsST2.getWeightedHashedFingerPrint());
            double similarityOfBondChanges = getPointWiseMutualInformation(bondFeatures1, bondFeatures2);

            double[] structFeatures1 = StructFP1.getWeightedHashedFingerPrint();
            double[] structFeatures2 = StructFP2.getWeightedHashedFingerPrint();
            double similarityStruct = getPointWiseMutualInformation(structFeatures1, structFeatures2);

            score = getSimilarityScore(alpha, beta, similarityOfBondChanges, similarityStruct);

        }

        return score;
    }

    /**
     *
     * @param alpha weight for BondFP
     * @param beta weight for StructuralFP
     * @param StructFP1
     * @param BondsCF1 Reactant Bond Cleaved/Formed FP for Query Reaction
     * @param BondsCF2 Product Bond Cleaved/Formed FP for Target Reaction
     * @param BondsOC1 Reactant Bond Order/Change FP for Query Reaction
     * @param StructFP2
     * @param BondsOC2 Product Bond Order/Change FP for Target Reaction
     * @param BondsST1 Reactant Bond Stereo FP for Query Reaction
     * @param BondsST2 Product Bond Stereo FP for Query Reaction
     *
     * @Equation Unweighted comparison of patterns
     * @return returns true if first fingerprint is a subset of second
     * @throws CDKException
     * @throws Exception
     */
    public static boolean isSubset(double alpha, double beta, double[] StructFP1, double[] BondsCF1,
            double[] BondsOC1, double[] BondsST1, double[] StructFP2, double[] BondsCF2,
            double[] BondsOC2, double[] BondsST2) throws CDKException, Exception {

//    public static double getWeightedScore(double alpha, double beta, BitSet reactants1, BitSet reactants2, double[] BondsCF1, double[] BondsCF2, BitSet products1, BitSet products2, double[] BondsOC1, double[] BondsOC2) throws CDKException, Exception {
//    System.out.println("Calling weighted score");
        boolean score = false;

        int size1 = StructFP1.length;
        int size2 = StructFP2.length;

        if (size1 != size2) {
            throw new CDKException("Features vectors must be of the same length");

        } else if (alpha < 0. && beta < 0.) {
            throw new CDKException("both alpha and beta can't =< be zero");
        } else {
            boolean similarityOfBondChanges = true;
            boolean similarityStruct = true;
            if (alpha > 0.) {
                similarityOfBondChanges = isSubset(BondsCF1, BondsOC1, BondsST1, BondsCF2, BondsOC2, BondsST2);
            }
            if (beta > 0.) {
                double[] structFeatures1 = StructFP1;
                double[] structFeatures2 = StructFP2;
                similarityStruct = isSubset(structFeatures1, structFeatures2);
            }

            if (similarityOfBondChanges && similarityStruct) {
                score = true;
            }
        }

        return score;
    }

    /**
     *
     * @param BondsCF1 Reactant Bond Cleaved/Formed FP for Query Reaction
     * @param BondsCF2 Product Bond Cleaved/Formed FP for Target Reaction
     * @param BondsOC1 Reactant Bond Order/Change FP for Query Reaction
     * @param BondsOC2 Product Bond Order/Change FP for Target Reaction
     * @param BondsST1 Reactant Bond Stereo FP for Query Reaction
     * @param BondsST2 Product Bond Stereo FP for Query Reaction
     *
     * @Equation Unweighted comparison of patterns
     * @return returns true if first fingerprint is a subset of second
     * @throws CDKException
     * @throws Exception
     */
    public static boolean isSubset(double[] BondsCF1, double[] BondsOC1, double[] BondsST1, double[] BondsCF2,
            double[] BondsOC2, double[] BondsST2) throws CDKException, Exception {
        double[] bondFeatures1 = append(BondsCF1, BondsOC1, BondsST1);
        double[] bondFeatures2 = append(BondsCF2, BondsOC2, BondsST2);
        return isSubset(bondFeatures1, bondFeatures2);
    }

    /**
     *
     * @param query
     * @param target @Equation Unweighted comparison of patterns
     * @return returns true if query fingerprint is a subset of target
     * @throws CDKException
     */
    public static boolean isSubset(double[] query, double target[]) throws CDKException {

        if (query.length != target.length) {
            throw new CDKException("Unequal Fingerprint size can't be processed");
        }

        for (int i = 0; i < query.length; i++) {
            if (query[i] > target[i]) {
                return false;
            }
        }
        return true;

    }

    /**
     *
     * @param query
     * @param target @Equation Unweighted comparison of patterns
     * @return returns true if query fingerprint is a subset of target
     * @throws CDKException
     */
    public static boolean isSubset(IPatternFingerprinter query, IPatternFingerprinter target) throws CDKException {

        if (query.getFingerprintSize() != target.getFingerprintSize()) {
            throw new CDKException("Unequal Fingerprint size can't be processed");
        }
        BitSet q = query.getHashedFingerPrint();
        BitSet t = target.getHashedFingerPrint();

        for (int i = 0; i < q.length(); i++) {
            if (!(t.get(i) & q.get(i))) {
                return false;
            }
        }
        return true;
    }

    /**
     *
     * @param fp1
     * @param fp2
     * @return reaction rawScore between two reaction based on the reactant
     * product structure
     * @throws CDKException
     * @throws java.lang.Exception
     */
    public static double getSimilarity(IPatternFingerprinter fp1, IPatternFingerprinter fp2) throws CDKException, Exception {

//        System.out.println("Calling StructuralSimilarity");
        double score = 0.0;
        int size1 = fp1.getFingerprintSize();
        int size2 = fp2.getFingerprintSize();

        if (size1 != size2) {
            throw new CDKException("Features vectors must be of the same length");

        } else {

            double[] structFeatures1 = fp1.getWeightedHashedFingerPrint();
            double[] structFeatures2 = fp2.getWeightedHashedFingerPrint();
            score = getSimilarity(structFeatures1, structFeatures2);
        }
        return score;
    }

    /**
     *
     * @param bondFeatures1
     * @param bondFeatures2
     * @return
     * @throws CDKException
     */
    private static double getSimilarity(double[] bondFeatures1, double[] bondFeatures2) throws CDKException {
        double similarity = 0.0;

        if (bondFeatures1.length != bondFeatures2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }

        int n = bondFeatures1.length;
        double ab = 0.0;
        double a2 = 0.0;
        double b2 = 0.0;

        for (int i = 0; i < n; i++) {
            ab += bondFeatures1[i] * bondFeatures2[i];
            a2 += bondFeatures1[i] * bondFeatures1[i];
            b2 += bondFeatures2[i] * bondFeatures2[i];
        }

        if (a2 > 0.0 && b2 > 0.0) {
            similarity = ab / (a2 + b2 - ab);
        }
        return similarity;
    }

    // log2:  Logarithm base 2
    private static double log2(double d) {
        return log(d) / log(2.0);
    }

    // log10: Logarithm base 10
    private static double log10(double d) {
        return log(d) / log(10.0);
    }

    // logx: Logarithm base 10
    private static double logX(double value, double base) {
        return log(value) / log(base);
    }

    private static double getPointWiseMutualInformation(double[] bondFeatures1, double[] bondFeatures2) {

        BitSet q = new BitSet(bondFeatures1.length);
        BitSet t = new BitSet(bondFeatures2.length);

        for (int i = 0; i < bondFeatures2.length; i++) {
            if (bondFeatures1[i] > 0.) {
                q.set(i, true);
            } else {
                q.set(i, false);
            }
            if (bondFeatures2[i] > 0.) {
                t.set(i, true);
            } else {
                t.set(i, false);
            }
        }

        List<Double> Pxy = new ArrayList<>(4);

        for (int i = 0; i < 4; i++) {
            Pxy.add(i, 0.);
        }

        for (int i = 0; i < q.size(); i++) {
            if (q.get(i) && t.get(i)) {
                Pxy.add(3, Pxy.get(3) + 1);
            }
            if (q.get(i) && !t.get(i)) {
                Pxy.add(2, Pxy.get(2) + 1);
            }
            if (!q.get(i) && t.get(i)) {
                Pxy.add(1, Pxy.get(1) + 1);
            }
            if (!q.get(i) && !t.get(i)) {
                Pxy.add(0, Pxy.get(0) + 1);
            }
        }

        for (int i = 0; i < 4; i++) {
            Pxy.set(i, Pxy.get(i) / q.size());
        }

        double Px0 = Pxy.get(0) + Pxy.get(1);
        double Px1 = Pxy.get(2) + Pxy.get(3);
        double Py0 = Pxy.get(0) + Pxy.get(2);
        double Py1 = Pxy.get(1) + Pxy.get(3);

        double SI_0_0 = log2(Pxy.get(0) / (Px0 * Py0));
        double SI_0_1 = log2(Pxy.get(1) / (Px0 * Py1));
        double SI_1_0 = log2(Pxy.get(2) / (Px1 * Py0));
        double SI_1_1 = log2(Pxy.get(3) / (Px1 * Py1));

        double MI = Pxy.get(0) * SI_0_0
                + Pxy.get(1) * SI_0_1
                + Pxy.get(2) * SI_1_0
                + Pxy.get(3) * SI_1_1;

        return MI;
    }

    private ReactionSimilarityTool() {
    }
}
