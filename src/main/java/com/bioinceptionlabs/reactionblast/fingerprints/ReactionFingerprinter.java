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
package com.bioinceptionlabs.reactionblast.fingerprints;

import java.io.Serializable;
import java.util.BitSet;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static java.lang.Long.toHexString;
import static java.lang.Math.sqrt;
import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.fingerprint.CircularFingerprinter.CLASS_ECFP4;
import static org.openscience.cdk.geometry.GeometryUtil.has2DCoordinates;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import static org.openscience.smsd.MoleculeInitializer.initializeMolecule;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.Feature;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.IFeature;


/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ReactionFingerprinter implements Serializable {

    private static final long serialVersionUID = 7867867834118778L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(ReactionFingerprinter.class);

    /**
     *
     * @param molSet
     * @throws CDKException
     */
    private static IPatternFingerprinter getSumOfFingerprints(IAtomContainerSet molSet) throws CDKException, Exception {
        FingerprintGenerator molFingerprint = new FingerprintGenerator();
        IPatternFingerprinter fp = new PatternFingerprinter(FingerprintGenerator.getFingerprinterSize());
        for (IAtomContainer mol : molSet.atomContainers()) {
            BitSet booleanArray = molFingerprint.getFingerprint(mol);
            for (int i = 0; i < booleanArray.size(); i++) {
                if (booleanArray.get(i)) {
                    fp.add(new Feature(valueOf(i), 1.0));
                }
            }
        }
        return fp;
    }

    /**
     *
     * @param bondFeatures1
     * @param bondFeatures2
     * @return
     * @throws CDKException
     */
    private static IPatternFingerprinter summationPatterns(IPatternFingerprinter pattern1, IPatternFingerprinter pattern2) throws CDKException {

        PatternFingerprinter patternFingerprinter = null;
        if (pattern1 != null && pattern2 != null
                && pattern1.getFingerprintSize()
                == pattern2.getFingerprintSize()) {
            patternFingerprinter = new PatternFingerprinter(pattern1.getFingerprintSize());

            patternFingerprinter.add(pattern1);
            patternFingerprinter.add(pattern2);
        } else {
            throw new CDKException("Index < 0: ");
        }
        return patternFingerprinter;
    }

    /*
     * @param reaction
     * @return
     */
    /**
     *
     * @param reaction
     * @return
     */
    public static IReaction expandReactionAndRemoveHydrogens(IReaction reaction) {
        IReaction r = new Reaction();
        /*
        * imp. to set reactin ID
         */
        String rid = reaction.getID() == null ? toHexString(currentTimeMillis()).toUpperCase() : reaction.getID();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer ac1 = ac.getBuilder().newInstance(IAtomContainer.class, ac);
            String id = ac.getID() == null ? toHexString(currentTimeMillis()).toUpperCase() : ac.getID();
            Double reactantCoefficient = reaction.getReactantCoefficient(ac);
            try {
                try {
                    ac1 = removeHydrogensExceptSingleAndPreserveAtomID(ac1);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                initializeMolecule(ac1);
            } catch (CDKException ex) {
                LOGGER.debug("ERROR: while configuring the reaction");
            }
            ac1.setID(id);
            for (int i = 0; i < reactantCoefficient; i++) {
                r.addReactant(ac1, 1.0);
            }
        }
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer ac1 = ac.getBuilder().newInstance(IAtomContainer.class, ac);
            String id = ac.getID() == null ? toHexString(currentTimeMillis()).toUpperCase() : ac.getID();
            Double productCoefficient = reaction.getProductCoefficient(ac);

            try {
                try {
                    ac1 = removeHydrogensExceptSingleAndPreserveAtomID(ac1);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                initializeMolecule(ac1);
            } catch (CDKException ex) {
                LOGGER.debug("ERROR: while configuring the reaction");
            }
            ac1.setID(id);
            for (int i = 0; i < productCoefficient; i++) {
                r.addProduct(ac1, 1.0);
            }
        }
        r.setID(rid);
        return r;
    }
    private final IPatternFingerprinter reactionFingerprint;

    /**
     *
     * @param reaction
     * @throws CDKException
     */
    public ReactionFingerprinter(IReaction reaction) throws CDKException {
        IReaction r = expandReactionAndRemoveHydrogens(reaction);
        IPatternFingerprinter fpr = null;
        try {
            fpr = getSumOfFingerprints(r.getReactants());
        } catch (Exception ex) {
            LOGGER.debug("ERROR: while get SumOfFingerprints for Reactants");
        }
        IPatternFingerprinter fpp = null;
        try {
            fpp = getSumOfFingerprints(r.getProducts());
        } catch (Exception ex) {
            LOGGER.debug("ERROR: while get SumOfFingerprints for Products");
        }
        this.reactionFingerprint = summationPatterns(fpr, fpp);
        reactionFingerprint.setFingerprintID(r.getID());
    }

    /**
     *
     * @return
     */
    public IPatternFingerprinter getReactionStruturalFingerprint() {
        return this.reactionFingerprint;
    }

    // === Inner classes merged from separate files ===


    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static interface IFingerprintGenerator {

        /**
         *
         * @param mol
         * @return
         * @throws CDKException
         */
        BitSet getFingerprint(IAtomContainer mol) throws CDKException;
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class FingerprintGenerator implements IFingerprintGenerator {

        private final static ILoggingTool LOGGER
                = createLoggingTool(FingerprintGenerator.class);

        /**
         * Size of the fingerprint
         *
         * @return
         */
        public static int getFingerprinterSize() {
            return new CircularFingerprinter(CLASS_ECFP4).getSize();
        }

        //define the FINGERPRINT_SIZE of the fingerprint
        //NOTE: this should be a multiple of 64 and preferably not 1024 or 2048
        //as for these values we often get the random numbers for one-atom or
        //two-atom paths the same!
        final CircularFingerprinter fingerprinter;

        /**
         *
         */
        public FingerprintGenerator() {
            fingerprinter = new CircularFingerprinter(CLASS_ECFP4);
        }

        /**
         *
         * @param mol
         * @return
         * @throws CDKException
         */
        @Override
        public BitSet getFingerprint(IAtomContainer mol) throws CDKException {
            if (!has2DCoordinates(mol)) {
                StructureDiagramGenerator structureDiagramGenerator = new StructureDiagramGenerator();
                structureDiagramGenerator.setMolecule(mol, true);
                if (isConnected(mol)) {
                    structureDiagramGenerator.generateCoordinates();
                    mol = structureDiagramGenerator.getMolecule();
                } else {
                    LOGGER.debug("Disconnected components needs to be layout separately");
                }
            }
            return fingerprinter.getBitFingerprint(mol).asBitSet();
        }

    }



    /**
     *
     * @author Syed Asad Rahman, BioInception
     * @contact asad.rahman@bioinceptionlabs.com
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
    public static class Similarity {

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
        public static float getTanimotoSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
            BitSet bitset1 = (BitSet) Molecule1.clone();
            BitSet bitset2 = (BitSet) Molecule2.clone();

            float _bitset1_cardinality = bitset1.cardinality();
            float _bitset2_cardinality = bitset2.cardinality();

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
        public static double getCosineSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
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
        public static double getDiceSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {
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


}
