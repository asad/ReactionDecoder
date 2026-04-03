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
package com.bioinceptionlabs.reactionblast.mechanism;

import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.Feature;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mapping.ReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.mapping.Reactor;
import com.bioinceptionlabs.reactionblast.mapping.SmsdReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.signature.RBlastMoleculeSignature;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.silent.RingSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.MoleculeInitializer;
import org.openscience.smsd.BaseMapping;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.CASE_INSENSITIVE_ORDER;
import static java.util.Arrays.sort;
import static java.util.Collections.sort;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import static org.openscience.cdk.smiles.smarts.parser.SMARTSParser.parse;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import static org.openscience.smsd.ExtAtomContainerManipulator.aromatizeDayLight;
import static org.openscience.smsd.ExtAtomContainerManipulator.cloneWithIDs;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 * Consolidated mechanism helper classes.
 * Merges: MatrixPrinter, Utility, CountSubstructures, AtomAtomMappingContainer,
 *         AtomStereoChangeInformation, BondChange, MoleculeMoleculePair,
 *         ReactantProductPair, ReactionCenterFragment
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class MechanismHelpers {

    private static final ReactionMappingEngine MAPPING_ENGINE
            = SmsdReactionMappingEngine.getInstance();

    private MechanismHelpers() { /* utility class */ }


    // ========== MatrixPrinter ==========

    public static class MatrixPrinter extends Object {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(MatrixPrinter.class);


        /**
         * This method prints the matrix to the standard output
         *
         * @param rMatrix R-Matrix to be Printed
         */
        public static void printReactionMatrix(RMatrix rMatrix) {
            try {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                    sb.append("\t\t").append(i);
                }
                sb.append(System.lineSeparator());
                for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                    sb.append("\t\t").append(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol())
                            .append(rMatrix.getReactantBEMatrix().getAtom(i).getID());
                }
                sb.append(System.lineSeparator());
                for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                    sb.append("\t\t").append(rMatrix.getProductBEMatrix().getAtom(i).getSymbol())
                            .append(rMatrix.getProductBEMatrix().getAtom(i).getID());
                }
                sb.append(System.lineSeparator());
                for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                    if (i == rMatrix.getRowDimension() - 1) {
                        sb.append("\t");
                    } else {
                        sb.append(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol())
                                .append(rMatrix.getReactantBEMatrix().getAtom(i).getID())
                                .append("\t").append(rMatrix.getProductBEMatrix().getAtom(i).getSymbol())
                                .append(rMatrix.getProductBEMatrix().getAtom(i).getID());
                    }
                    for (int j = 0; j < rMatrix.getColumnDimension(); j++) {
                        sb.append("\t").append(rMatrix.getValue(i, j));
                    }
                    sb.append(System.lineSeparator());
                }
                LOGGER.debug(sb.toString());
            } catch (CDKException ex) {
                LOGGER.debug("A CDKException has been arisen while printing the RMatrix");
            }
        }

        /**
         *
         * @param outputFile
         * @param rMatrix R-Matrix
         * @throws IOException
         */
        public static void writeReactionMatrix(File outputFile, RMatrix rMatrix) throws IOException {
            try (BufferedWriter matrixFileWriter = new BufferedWriter(new FileWriter(outputFile))) {
                matrixFileWriter.newLine();
                try {
                    for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                        matrixFileWriter.write("\t" + i);
                    }
                    matrixFileWriter.newLine();
                    for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                        matrixFileWriter.write("\t" + rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                                + rMatrix.getReactantBEMatrix().getAtom(i).getID());
                    }
                    matrixFileWriter.newLine();
                    for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                        matrixFileWriter.write("\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                                + rMatrix.getProductBEMatrix().getAtom(i).getID());
                    }
                    matrixFileWriter.newLine();
                    for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                        if (i == rMatrix.getRowDimension() - 1) {
                            matrixFileWriter.write("\t");
                        } else {
                            matrixFileWriter.write(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                                    + rMatrix.getReactantBEMatrix().getAtom(i).getID()
                                    + "\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                                    + rMatrix.getProductBEMatrix().getAtom(i).getID());
                        }
                        for (int j = 0; j < rMatrix.getColumnDimension(); j++) {
                            matrixFileWriter.write("\t" + rMatrix.getValue(i, j));
                        }
                        matrixFileWriter.newLine();
                    }

                } catch (CDKException ex) {
                    LOGGER.debug("A CDKException has been arisen while printing the RMatrix");
                }
            }
        }

        /**
         * This method prints the matrix to the standard output
         *
         * @param beMatrix
         *
         */
        public static void printBEMatrix(BEMatrix beMatrix) {
            List<IAtom> atomArray = beMatrix.getAtoms();
            StringBuilder sb = new StringBuilder();
            sb.append(atomArray.size()).append(System.lineSeparator());
            for (int i = 0; i < atomArray.size(); i++) {
                sb.append(atomArray.get(i).getSymbol()).append(atomArray.get(i).getID()).append("\t");
            }
            sb.append(System.lineSeparator());
            for (int i = 0; i < beMatrix.getRowDimension(); i++) {
                for (int j = 0; j < beMatrix.getColumnDimension(); j++) {
                    sb.append(beMatrix.getValue(i, j)).append("\t");
                }
                sb.append(System.lineSeparator());
            }
            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param outputFile
         * @param beMatrix BE-Matrix
         * @throws IOException
         */
        public static void writeBEMatrix(File outputFile, BEMatrix beMatrix) throws IOException {
            List<IAtom> atomArray = beMatrix.getAtoms();
            try (BufferedWriter matrixFileWriter = new BufferedWriter(new FileWriter(outputFile))) {
                matrixFileWriter.newLine();
                matrixFileWriter.write(atomArray.size());
                matrixFileWriter.newLine();
                for (int i = 0; i < atomArray.size(); i++) {
                    matrixFileWriter.write(atomArray.get(i).getSymbol() + atomArray.get(i).getID() + "\t");
                }
                matrixFileWriter.newLine();
                for (int i = 0; i < beMatrix.getRowDimension(); i++) {
                    for (int j = 0; j < beMatrix.getColumnDimension(); j++) {
                        matrixFileWriter.write(beMatrix.getValue(i, j) + "\t");
                    }
                    matrixFileWriter.newLine();
                }
            }
        }

        MatrixPrinter() {
        }
    }


    // ========== Utility ==========

    public static abstract class Utility extends MatrixPrinter implements Serializable {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(Utility.class);

        /**
         * Used CDK to generate moiety
         *
         * @param reactant
         * @param product
         * @param remove_AAM
         * @return
         */
        public static String getMoietyAsSMILES(IAtomContainer reactant, IAtomContainer product, boolean remove_AAM) {
            AtomAtomMapping atomAtomMapping = new AtomAtomMapping(reactant, product);
            for (IAtom a : reactant.atoms()) {
                for (IAtom b : product.atoms()) {
                    if (a.getID() == null ? b.getID() == null : a.getID().equals(b.getID())) {
                        atomAtomMapping.put(a, b);//store mapping if they share IDs
                    }
                }
            }
            StringBuilder sb = new StringBuilder("");
            try {
                sb.append(getSMILES(atomAtomMapping.getCommonFragment(), remove_AAM));
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            return sb.toString();
        }

        /**
         * Used CDK to generate smikrs
         *
         * @param reaction
         * @param remove_AAM
         * @return
         */
        public static String getSMILES(IReaction reaction, boolean remove_AAM) {
            StringBuilder sb = new StringBuilder("");
            try {
                for (IAtomContainer mol : reaction.getReactants().atomContainers()) {
                    sb.append(getSMILES(mol, remove_AAM));
                }

                sb.append(">>");

                for (IAtomContainer mol : reaction.getProducts().atomContainers()) {
                    sb.append(getSMILES(mol, remove_AAM));
                }

            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            return sb.toString();
        }

        /**
         * Used CDK to generate smiles
         *
         * @param mol
         * @param atom
         * @param level
         * @param remove_AAM
         * @return
         * @throws Exception
         */
        public static String getCircularSMILES(
                IAtomContainer mol, IAtom atom, int level, boolean remove_AAM) throws Exception {
            int refAtom = getAtomIndexByID(mol, atom);
            IAtomContainer fragment = getCircularFragment(mol, refAtom, level);
            String smiles = getSMILES(fragment, remove_AAM);
            return smiles;
        }

        /**
         * Used CDK to generate smiles
         *
         * @param mol
         * @param remove_AAM
         * @return
         */
        public static String getSMILES(
                IAtomContainer mol, boolean remove_AAM) {
            String smiles = "";
            try {
                return new CDKSMILES(mol, true, remove_AAM).getCanonicalSMILES();
            } catch (CloneNotSupportedException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            return smiles;
        }

        /**
         *
         * @param mol
         * @return
         */
        protected static List<IAtom> getAtoms(
                IAtomContainer mol) {
            List<IAtom> atoms = new ArrayList<>(mol.getAtomCount());
            for (IAtom atom : mol.atoms()) {
                atoms.add(atom);
            }
            return atoms;
        }

        /**
         *
         * @param bond
         * @param molset
         * @return
         */
        protected static String getMoleculeID(IBond bond, IAtomContainerSet molset) {
            for (IAtomContainer mol : molset.atomContainers()) {
                if (mol.contains(bond)) {
                    return mol.getID();
                }
            }
            return null;
        }

        /**
         *
         * @param atom
         * @param molset
         * @return
         */
        protected static String getMoleculeID(IAtom atom, IAtomContainerSet molset) {
            for (IAtomContainer mol : molset.atomContainers()) {
                if (mol.contains(atom)) {
                    return mol.getID();
                }
            }
            return null;
        }

        /**
         *
         * @param bond
         * @param molset
         * @return
         */
        protected static IAtomContainer getAtomContainer(IBond bond, IAtomContainerSet molset) {
            for (IAtomContainer mol : molset.atomContainers()) {
                if (mol.contains(bond)) {
                    return mol;
                }
            }
            return null;
        }

        /**
         *
         * @param atom
         * @param molset
         * @return
         */
        protected static IAtomContainer getAtomContainer(IAtom atom, IAtomContainerSet molset) {
            for (IAtomContainer mol : molset.atomContainers()) {
                if (mol.contains(atom)) {
                    return mol;
                }
            }
            return null;
        }

        /**
         * Return atom by ID match
         *
         * @param molWithoutH
         * @param refAtom
         * @return
         */
        protected static int getAtomIndexByID(IAtomContainer molWithoutH, IAtom refAtom) {
            for (IAtom atom : molWithoutH.atoms()) {
                if (atom.getID().equalsIgnoreCase(refAtom.getID())) {
                    return molWithoutH.indexOf(atom);
                }
            }
            return -1;
        }

        /**
         * Return Signature of height h
         *
         * @param mol
         * @param atom
         * @param height
         * @return
         * @throws CloneNotSupportedException
         */
        protected static String getSignature(IAtomContainer mol, IAtom atom, int height) throws CloneNotSupportedException {
            IAtomContainer molWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(mol);
            int atomIndex = getAtomIndexByID(molWithoutH, atom);
            RBlastMoleculeSignature moleculeSignature = new RBlastMoleculeSignature(molWithoutH);
            moleculeSignature.setUseCharge(true);
            moleculeSignature.setBondSensitive(true);
            moleculeSignature.setUseAromatics(true);
            if (atomIndex >= 0) {
                return moleculeSignature.getAtomSignature(atomIndex, height).toCanonicalString();
            } else {
                return "";
            }
        }

        /**
         *
         * @param atomRCChangesMap
         * @param fragments
         * @throws CloneNotSupportedException
         */
        protected static void setFragmentMatches(SortedMap<String, Integer> atomRCChangesMap, List<IAtomContainer> fragments) throws CloneNotSupportedException {
            for (IAtomContainer fragment : fragments) {
                CountSubstructures countSubstructures = new CountSubstructures(fragment);
                atomRCChangesMap.keySet().stream().forEach((pattern) -> {
                    int hit = 0;
                    try {
                        hit = countSubstructures.substructureSize(pattern);
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                    int val = hit == 0 ? 0 : atomRCChangesMap.get(pattern) + 1;
                    atomRCChangesMap.put(pattern, val);
                });
            }
        }

        /**
         *
         * @param atomRCChangesMap
         * @param signatures
         */
        protected static void setReactionCenterMatches(IPatternFingerprinter atomRCChangesMap, List<String> signatures) {
            signatures.stream().forEach((fragment) -> {
                try {
                    atomRCChangesMap.add(new Feature(fragment, 1.0));
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            });
        }

        /**
         *
         * @param rid
         * @param mol
         * @param atom
         * @param patternFP
         * @throws CDKException
         * @throws CloneNotSupportedException
         */
        protected static void setCircularSignatureFingerprints(String rid, IAtomContainer mol, IAtom atom, Map<Integer, IPatternFingerprinter> patternFP) throws CDKException, CloneNotSupportedException {

            for (int i = 1; i < 5; i++) {
                if (!patternFP.containsKey(i)) {
                    IPatternFingerprinter fp = new PatternFingerprinter();
                    fp.setFingerprintID(rid + ":" + "Signature: " + i);
                    patternFP.put(i, fp);
                }
                String signature = getSignature(mol, atom, i);
                patternFP.get(i).add(new Feature(signature, 1.0));
            }

            if (!patternFP.containsKey(-1)) {
                IPatternFingerprinter fp = new PatternFingerprinter();
                fp.setFingerprintID(rid + ":" + "Signature: " + -1);
                patternFP.put(-1, fp);
            }
            String signature = getSignature(mol, atom, -1);
            patternFP.get(-1).add(new Feature(signature, 1.0));
        }

        /**
         *
         * @param reactBond
         * @param prodBond
         * @return
         */
        protected static String getCanonisedBondChangePattern(IBond reactBond, IBond prodBond) {

            String concatE = getCanonicalisedBondChangePattern(reactBond);
            String concatP = getCanonicalisedBondChangePattern(prodBond);

            List<String> pattern = new ArrayList<>(2);
            pattern.add(0, concatE);
            pattern.add(1, concatP);

            sort(pattern, CASE_INSENSITIVE_ORDER);
            return pattern.get(0).concat("*").concat(pattern.get(1));

        }

        /**
         *
         * @param bond
         * @return
         */
        protected static String getCanonicalisedBondChangePattern(IBond bond) {
            String symbol = getBondOrderSign(bond);
            List<String> atoms = new ArrayList<>(2);
            atoms.add(0, bond.getAtom(0).getSymbol());
            atoms.add(1, bond.getAtom(1).getSymbol());
            sort(atoms, CASE_INSENSITIVE_ORDER);
            String concatenatedSymbols = atoms.get(0).concat(symbol).concat(atoms.get(1));
            return concatenatedSymbols.trim();
        }

        /**
         *
         * @param bond
         * @return
         */
        public static String getBondOrderSign(IBond bond) {
            String bondSymbol = "";
            if (bond.isAromatic()) {
                bondSymbol += "@";
            } else if (bond.isInRing()) {
                bondSymbol += "%";
            } else if (bond.getOrder() == SINGLE) {
                bondSymbol += "-";
            } else if (bond.getOrder() == DOUBLE) {
                bondSymbol += "=";
            } else if (bond.getOrder() == TRIPLE) {
                bondSymbol += "#";
            } else if (bond.getOrder() == QUADRUPLE) {
                return "$";
            }

            return bondSymbol;
        }

        /**
         *
         * @param ringBond
         * @param singleRings
         * @return
         */
        public static int getNeighbourBondOrderCountFromRing(IBond ringBond, IRingSet singleRings) {
            int minValue = 9999;
            for (IAtomContainer ring : singleRings.atomContainers()) {
                int value = 0;
                if (ring.contains(ringBond.getAtom(0)) && ring.contains(ringBond.getAtom(1))) {
                    for (IBond bond : ring.bonds()) {
                        if (bond.contains(ringBond.getAtom(0)) || bond.contains(ringBond.getAtom(1))) {
                            value += bond.getOrder().numeric();
                        }
                    }
                }
                if (value < minValue) {
                    minValue = value;
                }
            }
            return minValue;
        }

        /**
         *
         * @param ringBond
         * @param singleRings
         * @return
         */
        public static IRingSet getSmallestRingSet(IBond ringBond, IRingSet singleRings) {
            IRingSet rs = new RingSet();
            for (IAtomContainer ring : singleRings.atomContainers()) {
                if (ring.contains(ringBond.getAtom(0)) && ring.contains(ringBond.getAtom(1))) {
                    if (rs.getAtomContainerCount() == 0) {
                        rs.addAtomContainer(ring);
                        continue;
                    }
                    for (IAtomContainer smallestRing : rs.atomContainers()) {
                        if (ring.getAtomCount() == smallestRing.getAtomCount()) {
                            if (!rs.contains(ring)) {
                                rs.addAtomContainer(ring);
                            }
                        } else if (ring.getAtomCount() < smallestRing.getAtomCount()) {
                            rs.removeAllAtomContainers();
                            rs.addAtomContainer(ring);
                        }
                    }
                }
            }
            return rs;
        }

        /**
         *
         * @param atomContainer Atom container where rings are to be marked
         * @throws CDKException if there is a problem in ring perception or
         * aromaticity detection, which is usually related to a timeout in the ring
         * finding code.
         */
        protected static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
            MoleculeInitializer.initializeMolecule(atomContainer);
        }

        /**
         *
         * @param mol
         * @param startAtomIndex
         * @param radius
         * @return
         * @throws Exception
         */
        public static IAtomContainer getCircularFragment(IAtomContainer mol, int startAtomIndex, int radius) throws Exception {
            IAtomContainer fragment = cloneWithIDs(mol);
            Set<IAtom> removeList = new HashSet<>();
            Collection<IAtom> solutionSphereList = circularFragment(fragment, startAtomIndex, radius);

            for (IAtom atom : fragment.atoms()) {
                if (!solutionSphereList.contains(atom)) {
                    removeList.add(atom);
                }
            }

            for (Iterator<IAtom> it = removeList.iterator(); it.hasNext();) {
                fragment.removeAtom(it.next());
            }

            IAtomContainer canonicalise = canonicalise(fragment);
            aromatizeDayLight(canonicalise);

            return fragment;
        }

        /**
         *
         * @param org_mol
         * @return cloned canonicalized molecule
         * @throws CloneNotSupportedException
         * @throws CDKException
         */
        public static IAtomContainer canonicalise(IAtomContainer org_mol) throws CloneNotSupportedException, CDKException {

            IAtomContainer cloneMolecule = cloneWithIDs(org_mol);

            int[] p = new int[cloneMolecule.getAtomCount()];
    //        /*
    //         Signature based canonical Permutations
    //         */
    //        p = new com.bioinceptionlabs.reactionblast.tools.
    //        SignatureMoleculeLabeller().getCanonicalPermutation(cloneMolecule);

            /*
            Use the Canonical labelling from the SMILES
            IMP: Suggested by John May
             */
            try {
    //            unique().create(cloneMolecule, p);
                SmilesGenerator smiles = new SmilesGenerator(
                        SmiFlavor.AtomAtomMap
                        | SmiFlavor.Unique
                        //                    | SmiFlavor.UseAromaticSymbols
                        | SmiFlavor.Stereo);

                String sm = smiles.create(cloneMolecule, p);
            } catch (Exception e) {
                LOGGER.debug("Fragment not fit to canonicalise: " + e.getMessage());
            }

            permuteWithoutClone(p, cloneMolecule);

            /*
            Set the IDs to container
             */
            if (org_mol.getID() != null) {
                cloneMolecule.setID(org_mol.getID());
            }

            return cloneMolecule;
        }

        /*
        This is a very imp code modified by John May
        The idea is to canonicalise the atoms and bonds
         */
        private static void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
            int n = atomContainer.getAtomCount();
            int[] permutation = normalizePermutation(p, n);
            IAtom[] permutedAtoms = new IAtom[n];

            for (int i = 0; i < n; i++) {
                IAtom atom = atomContainer.getAtom(i);
                permutedAtoms[permutation[i]] = atom;
                atom.setProperty("label", permutation[i]);
            }
            atomContainer.setAtoms(permutedAtoms);

            IBond[] bonds = getBondArray(atomContainer);
            sort(bonds, (IBond o1, IBond o2) -> {
                int u = o1.getAtom(0).getProperty("label");
                int v = o1.getAtom(1).getProperty("label");
                int x = o2.getAtom(0).getProperty("label");
                int y = o2.getAtom(1).getProperty("label");
                int min1 = min(u, v);
                int min2 = min(x, y);
                int max1 = max(u, v);
                int max2 = max(x, y);

                int minCmp = Integer.compare(min1, min2);
                if (minCmp != 0) {
                    return minCmp;
                }
                int maxCmp = Integer.compare(max1, max2);
                if (maxCmp != 0) {
                    return maxCmp;
                }
                LOGGER.debug("pokemon!");
                throw new InternalError();
            });
            atomContainer.setBonds(bonds);
        }

        private static int[] normalizePermutation(int[] permutation, int size) {
            if (permutation == null || permutation.length != size) {
                return identityPermutation(size);
            }

            boolean[] seen = new boolean[size];
            for (int value : permutation) {
                if (value < 0 || value >= size || seen[value]) {
                    return identityPermutation(size);
                }
                seen[value] = true;
            }
            return permutation;
        }

        private static int[] identityPermutation(int size) {
            int[] identity = new int[size];
            for (int i = 0; i < size; i++) {
                identity[i] = i;
            }
            return identity;
        }

        /**
         * Performs a breadthFirstSearch in an AtomContainer starting with a
         * particular sphere, which usually consists of one start atom. While
         * searching the graph, the method marks each visited atom. It then puts all
         * the atoms connected to the atoms in the given sphere into a new vector
         * which forms the sphere to search for the next recursive method call. All
         * atoms that have been visited are put into a molecule container. This
         * breadthFirstSearch does thus find the connected graph for a given start
         * atom.
         *
         * @param atomContainer The AtomContainer to be searched
         * @param rootAtom
         * @param max
         * @return
         * @throws CDKException
         */
        public static Collection<IAtom> circularFragment(IAtomContainer atomContainer, int rootAtom, int max) throws CDKException {
            IAtom root = atomContainer.getAtom(rootAtom);
            Set<IAtom> paths = new HashSet<>();
            // list of visited nodes
            LinkedList<IAtom> closedList = new LinkedList<>();
            // list of nodes to visit (sorted)
            LinkedList<IAtom> openList = new LinkedList<>();
            openList.add(root);

            // list of nodes to visit (sorted)
            LinkedList<IAtom> neighbours = new LinkedList<>();

            int level = 0;
            while (!openList.isEmpty()) {
                IAtom currentPath = openList.removeFirst();

                // path found!
                paths.add(currentPath);
                closedList.add(currentPath);

                // addBinary neighbors to the open list
                neighbours.addAll(atomContainer.getConnectedAtomsList(currentPath));

                if (openList.isEmpty() && !neighbours.isEmpty() && (max > level || max == -1)) {
                    neighbours.stream().filter((a) -> (!closedList.contains(a))).forEach((a) -> {
                        openList.add(a);
                    });
                    level += 1;
                    neighbours.clear();
                }
            }
            return paths;
        }

        /**
         *
         * @param rid
         * @param molOrignal
         * @param atom
         * @param patternFP
         * @throws Exception
         * @throws CloneNotSupportedException
         */
        protected static void setCircularFingerprints(String rid,
                IAtomContainer molOrignal,
                IAtom atom, Map<Integer, IPatternFingerprinter> patternFP)
                throws Exception, CloneNotSupportedException {
            IAtomContainer clone = molOrignal.clone();
            for (int i = 0; i < 3; i++) {
                if (!patternFP.containsKey(i)) {
                    IPatternFingerprinter fp = new PatternFingerprinter();
                    fp.setFingerprintID(rid + ":" + "Signature: " + i);
                    patternFP.put(i, fp);
                }
                String circularSMILES = getCircularSMILES(clone, atom, i, true);
                patternFP.get(i).add(new Feature(circularSMILES, 1.0));
            }
            if (!patternFP.containsKey(-1)) {
                IPatternFingerprinter fp = new PatternFingerprinter();
                fp.setFingerprintID(rid + ":" + "Signature: " + -1);
                patternFP.put(-1, fp);
            }

            String circularSMILES = getCircularSMILES(clone, atom, -1, true);
            patternFP.get(-1).add(new Feature(circularSMILES, 1.0));

        }

        /**
         *
         * @param molOrignal
         * @param type
         * @param atom
         * @return
         * @throws Exception
         * @throws CloneNotSupportedException
         */
        protected static List<ReactionCenterFragment> getCircularReactionPatternFingerprints(IAtomContainer molOrignal,
                IAtom atom,
                BondChangeCalculator.EnumSubstrateProduct type)
                throws Exception, CloneNotSupportedException {
            List<ReactionCenterFragment> fragmentsRC = new ArrayList<>();
            IAtomContainer clone = molOrignal.clone();
            for (int i = 0; i < 3; i++) {
                String smiles = getCircularSMILES(clone, atom, i, true);
                ReactionCenterFragment reactionCenterFragment = new ReactionCenterFragment(smiles, i, type);
                fragmentsRC.add(reactionCenterFragment);
            }
            String smiles = getCircularSMILES(clone, atom, -1, true);
            ReactionCenterFragment reactionCenterFragment = new ReactionCenterFragment(smiles, -1, type);
            fragmentsRC.add(reactionCenterFragment);
            return fragmentsRC;
        }
    }


    // ========== CountSubstructures ==========

    static class CountSubstructures extends MoleculeInitializer implements Serializable {

        private static final ILoggingTool LOGGER
                = createLoggingTool(CountSubstructures.class);
        private static final long serialVersionUID = 12343289751445148L;
        private final SmilesParser sp;
        private IAtomContainer mol;

        CountSubstructures(IAtomContainer atomContainer) throws CloneNotSupportedException {
            sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
            try {
                this.mol = null;
                mol = removeHydrogensExceptSingleAndPreserveAtomID(atomContainer);
                initializeMolecule(mol);
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        public int substructureSize(String smiles) throws CDKException {
            AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
            BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

            try {
                IAtomContainer parseSmiles = sp.parseSmiles(smiles);
                BaseMapping sub = MAPPING_ENGINE.findSubstructure(
                        parseSmiles, mol, atomMatcher, bondMatcher, false);
                return sub.isSubgraph() ? sub.getFirstAtomMapping().getCount() : 0;
            } catch (InvalidSmilesException ex) {
                BaseMapping sub = MAPPING_ENGINE.findSubstructure(
                        parse(smiles, mol.getBuilder()), mol, false);
                return sub.isSubgraph() ? sub.getFirstAtomMapping().getCount() : 0;
            }
        }
    }


    // ========== AtomAtomMappingContainer ==========

    public static class AtomAtomMappingContainer extends Object implements Serializable {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(AtomAtomMappingContainer.class);

        private static final String NEW_LINE = System.lineSeparator();
        private static final long serialVersionUID = 17879096958755L;

        private List<IAtom> reactantAtomArray = new ArrayList<>();
        private List<IAtom> productAtomArray = new ArrayList<>();
    //    private Reactor myReaction = null;

        /**
         * Class constructor. Creates the mapping of a given reaction.
         *
         * @param reactor Reactor for which the AtomAtomMappingContainer is required
         * @param withoutH Store Mapping without H
         * @throws Exception
         */
        public AtomAtomMappingContainer(Reactor reactor, boolean withoutH) throws Exception {
            this(reactor.getReactionWithAtomAtomMapping(), withoutH);
        }

        /**
         *
         * @param reactants
         * @param products
         * @param withoutH
         */
        public AtomAtomMappingContainer(IAtomContainerSet reactants, IAtomContainerSet products, boolean withoutH) {
            int atomNo = 0;
            int mappedAtomsR = 0;
            int mappedAtomsP = 0;
            IAtom[] atVect = null;

            //REACTANTS
            if (withoutH) {
                for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                    IAtomContainer M = reactants.getAtomContainer(i);
                    for (IAtom a : M.atoms()) {
                        if (!a.getSymbol().equalsIgnoreCase("H")) {
                            atomNo++;
                        }
                    }
                }
            } else {
                for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                    atomNo += reactants.getAtomContainer(i).getAtomCount();
                }
            }
            atVect = new IAtom[atomNo];
            for (int i = 0; i < atVect.length; i++) {
                atVect[i] = null;
            }
            for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                for (int j = 0; j < reactants.getAtomContainer(i).getAtomCount(); j++) {
                    IAtom at = reactants.getAtomContainer(i).getAtom(j);
                    if (withoutH && at.getSymbol().equalsIgnoreCase("H")) {
                        continue;
                    }
                    if (at.getID() == null) {
                        continue;
                    }
                    int atomID;
                    try {
                        atomID = Integer.valueOf(at.getID());
                    } catch (NumberFormatException e) {
                        continue;
                    }
                    if (atomID <= 0) {
                        continue;
                    }
                    if (atomID - 1 >= atVect.length) {
                        continue;
                    }
                    atVect[atomID - 1] = at;
                    mappedAtomsR++;
                }
            }
            for (int i = 0; i < mappedAtomsR; i++) {
                reactantAtomArray.add(atVect[i]);
            }
            //Checking for holes in the vector. 
            boolean findNull = false;
            boolean error = false;
            for (IAtom atVect1 : atVect) {
                if (findNull && (atVect1 != null)) {
                    error = true;
                }
                if (atVect1 == null) {
                    findNull = true;
                }
            }
            if (error) {
                StringBuilder sb = new StringBuilder();
                sb.append("ERROR in AtomAtomMapping-found hole in the mapping (reactants atomIDs)");
                for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                    sb.append(NEW_LINE).append("Mol:").append(reactants.getAtomContainer(i).getID());
                    for (int j = 0; j < reactants.getAtomContainer(i).getAtomCount(); j++) {
                        IAtom at = reactants.getAtomContainer(i).getAtom(j);
                        sb.append(NEW_LINE).append(at.getSymbol()).append(at.getID());
                    }
                }
                LOGGER.error(sb.toString());
            }
            //end of checking statements

            //PRODUCTS
            atomNo = 0;
            if (withoutH) {
                for (int i = 0; i < products.getAtomContainerCount(); i++) {
                    IAtomContainer M = products.getAtomContainer(i);
                    for (IAtom a : M.atoms()) {
                        if (!a.getSymbol().equalsIgnoreCase("H")) {
                            atomNo++;
                        }
                    }
                }
            } else {
                for (int i = 0; i < products.getAtomContainerCount(); i++) {
                    atomNo += products.getAtomContainer(i).getAtomCount();
                }
            }
            atVect = new IAtom[atomNo];
            for (int i = 0; i < atVect.length; i++) {
                atVect[i] = null;
            }
            for (int i = 0; i < products.getAtomContainerCount(); i++) {
                for (int j = 0; j < products.getAtomContainer(i).getAtomCount(); j++) {
                    IAtom at = products.getAtomContainer(i).getAtom(j);
                    if (withoutH && at.getSymbol().equalsIgnoreCase("H")) {
                        continue;
                    }

                    if (at.getID() == null) {
                        continue;
                    }
                    int atomID;
                    try {
                        atomID = Integer.valueOf(at.getID());
                    } catch (NumberFormatException e) {
                        continue;
                    }
                    if (atomID <= 0) {
                        continue;
                    }
                    if (atomID - 1 >= atVect.length) {
                        continue;
                    }
                    atVect[atomID - 1] = at;
                    mappedAtomsP++;
                }
            }
            for (int i = 0; i < mappedAtomsP; i++) {
                productAtomArray.add(atVect[i]);
            }
            //Checking for holes in the vector. 
            findNull = false;
            error = false;
            for (IAtom atVect1 : atVect) {
                if (findNull && (atVect1 != null)) {
                    error = true;
                }
                if (atVect1 == null) {
                    findNull = true;
                }
            }
            if (mappedAtomsP != mappedAtomsR) {
                error = true;
            }
            if (error) {
                StringBuilder sb = new StringBuilder();
                sb.append("ERROR in AtomAtomMapping-found hole in the mapping (products atomIDs)");
                sb.append("mapped reactants atoms: ").append(mappedAtomsR).append(", mapped products atoms: ").append(mappedAtomsP);
                for (int i = 0; i < products.getAtomContainerCount(); i++) {
                    sb.append(NEW_LINE).append("Mol:").append(products.getAtomContainer(i).getID());
                    for (int j = 0; j < products.getAtomContainer(i).getAtomCount(); j++) {
                        IAtom at = products.getAtomContainer(i).getAtom(j);
                        sb.append(NEW_LINE).append(at.getSymbol()).append(at.getID());
                    }
                }
                LOGGER.error(sb.toString());
            }
            //end of checking statements
        }

        /**
         *
         * @param reaction
         * @param withoutH
         */
        public AtomAtomMappingContainer(IReaction reaction, boolean withoutH) {
            for (IMapping m : reaction.mappings()) {
                IAtom rAtom = (IAtom) m.getChemObject(0);
                IAtom pAtom = (IAtom) m.getChemObject(1);
                if (withoutH && rAtom != null && pAtom != null
                        && (rAtom.getSymbol().equalsIgnoreCase("H")
                        || pAtom.getSymbol().equalsIgnoreCase("H"))) {
                } else {
                    reactantAtomArray.add(rAtom);
                    productAtomArray.add(pAtom);
                }
            }
        }

        /**
         * This method prints the matrix to the standard output
         *
         * @return
         */
        @Override
        public String toString() {
            StringBuilder result = new StringBuilder();
            result.append(reactantAtomArray.size()).append(NEW_LINE);
            for (int i = 0; i < reactantAtomArray.size(); i++) {
                result.append(i).append("\t");
            }
            result.append(NEW_LINE);
            for (int i = 0; i < reactantAtomArray.size(); i++) {
                result.append((reactantAtomArray.get(i)).getSymbol()).append((reactantAtomArray.get(i)).getID()).append("\t");
            }
            result.append(NEW_LINE);
            for (int i = 0; i < productAtomArray.size(); i++) {
                result.append((productAtomArray.get(i)).getSymbol()).append((productAtomArray.get(i)).getID()).append("\t");
            }
            result.append(NEW_LINE);
            return result.toString();
        }

        /**
         * The method returns the product atom mapped to the reactant atom passed as
         * parameter. Returns null if the reactantAtom is not mapped to any product
         * atom.
         *
         * @param reactantAtom The IAtom for which the product atom is required.
         * @return The product atom mapped to the given reactant atom.
         */
        public IAtom getMappedProductAtom(IAtom reactantAtom) {
            IAtom a = null;
            int reactantIdx = -1;
            for (int i = 0; i < reactantAtomArray.size(); i++) {
                if (reactantAtomArray.get(i).getID().equals(reactantAtom.getID())) {
                    reactantIdx = i;
                }
            }
            if (reactantIdx != -1) {
                a = productAtomArray.get(reactantIdx);
            }
            return a;
        }

        /**
         * The method returns the idx-th reactant atom which has been mapped.
         *
         * @param idx The index of the reactant atom which is required.
         * @return The idx-th a reactant atom mapped.
         */
        public IAtom getReactantAtom(int idx) {
            IAtom ret = null;
            if ((idx < reactantAtomArray.size()) && (idx > -1)) {
                ret = reactantAtomArray.get(idx);
            }
            return ret;
        }

        /**
         *
         * @param idx
         * @return
         */
        public IAtom getProductAtom(int idx) {
            IAtom ret = null;
            if ((idx < productAtomArray.size()) && (idx > -1)) {
                ret = productAtomArray.get(idx);
            }
            return ret;
        }

        /**
         * Returns the number of mappings which the AtomAtomMappingContainer
         * contains.
         *
         * @return the number of mappings which the AtomAtomMappingContainer
         * contains.
         */
        public int getSize() {
            return reactantAtomArray.size();
        }

        /**
         *
         * @return
         */
        public int getSizeNoHydrogens() {
            int count = 0;
            count = reactantAtomArray.stream().filter((a)
                    -> (!a.getSymbol().equals("H"))).map((_item) -> 1)
                    .reduce(count, Integer::sum);
            return count;
        }

        /**
         * Returns true if the reactant atom is present
         *
         * @param atom
         * @return
         */
        public boolean isReactantAtomPresent(IAtom atom) {
            return reactantAtomArray.contains(atom) == true;
        }

        /**
         * Return true if the product atom is present
         *
         * @param atom
         * @return
         */
        public boolean isProductAtomPresent(IAtom atom) {
            return productAtomArray.contains(atom) == true;
        }
    }


    // ========== AtomStereoChangeInformation ==========

    public static class AtomStereoChangeInformation implements Serializable {

        private static final long serialVersionUID = 1896986585959789L;
        private final IAtom reactantAtom;
        private final IAtom productAtom;
        private boolean stereoChange = false;
        private BondChangeCalculator.IStereoAndConformation atomStereoR;
        private BondChangeCalculator.IStereoAndConformation atomStereoP;

        /**
         *
         * @param rAtom
         * @param pAtom
         */
        public AtomStereoChangeInformation(IAtom rAtom, IAtom pAtom) {
            this.reactantAtom = rAtom;
            this.productAtom = pAtom;
            setStereoChange(true);
        }

        /**
         *
         * @param atomE
         * @param atomP
         * @param aStereoR
         * @param aStereoP
         */
        public AtomStereoChangeInformation(IAtom atomE, IAtom atomP, BondChangeCalculator.IStereoAndConformation aStereoR, BondChangeCalculator.IStereoAndConformation aStereoP) {
            this(atomE, atomP);
            this.atomStereoR = aStereoR;
            this.atomStereoP = aStereoP;
        }

        /**
         * @return the reactantAtom
         */
        public IAtom getReactantAtom() {
            return reactantAtom;
        }

        /**
         * @return the productAtom
         */
        public IAtom getProductAtom() {
            return productAtom;
        }

        /**
         * @return the stereoChange
         */
        public boolean isStereoChange() {
            return stereoChange;
        }

        /**
         * @param stereoChange the stereoChange to set
         */
        private void setStereoChange(boolean stereoChange) {
            this.stereoChange = stereoChange;
        }

        /**
         * @return the atomStereo
         */
        public BondChangeCalculator.IStereoAndConformation getReactantAtomStereo() {
            return atomStereoR;
        }

        /**
         * @return the atomStereo
         */
        public BondChangeCalculator.IStereoAndConformation getProductAtomStereo() {
            return atomStereoP;
        }

        /**
         * @param atomStereoR
         * @param atomStereoP
         */
        public void setAtomStereo(BondChangeCalculator.IStereoAndConformation atomStereoR, BondChangeCalculator.IStereoAndConformation atomStereoP) {
            this.atomStereoR = atomStereoR;
            this.atomStereoP = atomStereoP;
        }
    }


    // ========== BondChange ==========

    public static class BondChange implements Serializable {

        private static final String NEW_LINE = System.lineSeparator();
        private static final long serialVersionUID = 9890766688070991L;

        /**
         *
         * @param bond
         * @return
         */
        public static int convertBondOrder(IBond bond) {
            if (bond.getOrder() == null) {
                return bond.isAromatic() ? 2 : 1;
            }
            switch (bond.getOrder()) {
                case QUADRUPLE:
                    return 4;
                case TRIPLE:
                    return 3;
                case DOUBLE:
                    return 2;
                case SINGLE:
                    return 1;
                default:
                    // Handle UNSET or other cases — check aromaticity flag
                    return bond.isAromatic() ? 2 : 1;
            }
        }

        /**
         *
         * @param bond
         * @return
         */
        @SuppressWarnings("deprecation")
        public static int convertBondStereo(IBond bond) {
            int value;
            switch (bond.getStereo()) {
                case UP:
                    value = 1;
                    break;
                case UP_INVERTED:
                    value = 1;
                    break;
                case DOWN:
                    value = 6;
                    break;
                case DOWN_INVERTED:
                    value = 6;
                    break;
                case UP_OR_DOWN:
                    value = 4;
                    break;
                case UP_OR_DOWN_INVERTED:
                    value = 4;
                    break;
                case E_OR_Z:
                    value = 3;
                    break;
                default:
                    value = 0;
            }
            return value;
        }

        private final IBond reactantBond;
        private final IBond productBond;
        private final float bondChangeDelta;

        /**
         *
         * @param reactantBond
         * @param productBond
         */
        public BondChange(IBond reactantBond, IBond productBond) {
            this.reactantBond = reactantBond;
            this.productBond = productBond;
            if (this.reactantBond != null && this.productBond != null) {
                this.bondChangeDelta = convertBondOrder(this.productBond) - convertBondOrder(this.reactantBond);
            } else if (this.reactantBond == null && this.productBond != null) {
                this.bondChangeDelta = convertBondOrder(this.productBond);
            } else if (this.reactantBond != null && this.productBond == null) {
                this.bondChangeDelta = convertBondOrder(this.reactantBond);
            } else {
                this.bondChangeDelta = 0;
            }
        }

        /**
         * @return the reactantBond
         */
        public IBond getReactantBond() {
            return reactantBond;
        }

        /**
         * @return the productBond
         */
        public IBond getProductBond() {
            return productBond;
        }

        /**
         * @return the bondChangeDelta
         */
        public float getBondChangeDelta() {
            return bondChangeDelta;
        }

        @Override
        public String toString() {
            StringBuilder result = new StringBuilder();
            result.append("\t");
            result.append(NEW_LINE);
            if (reactantBond != null) {
                result.append("R: ").append(reactantBond.getAtom(0).getSymbol());
                result.append("(").append(reactantBond.getAtom(0).getID()).append(")");
                result.append("[").append(convertBondOrder(reactantBond)).append("]");
                result.append(reactantBond.getAtom(1).getSymbol());
                result.append("(").append(reactantBond.getAtom(1).getID()).append(")");

            } else {
                result.append("NA");
            }

            if (productBond != null) {
                result.append(", P: ").append(productBond.getAtom(0).getSymbol());
                result.append("(").append(productBond.getAtom(0).getID()).append(")");
                result.append("[").append(convertBondOrder(productBond)).append("]");
                result.append(productBond.getAtom(1).getSymbol());
                result.append("(").append(productBond.getAtom(1).getID()).append(")");
                result.append(NEW_LINE);
            } else {
                result.append(", NA");
                result.append(NEW_LINE);
            }
            return result.toString();
        }

    }


    // ========== MoleculeMoleculePair ==========

    public static class MoleculeMoleculePair implements Serializable, Comparable<MoleculeMoleculePair>, Comparator<MoleculeMoleculePair> {

        private static final long serialVersionUID = 107097779868968L;
        private final ReactantProductPair name;
        private final ReactantProductPair smarts;
        private final ReactantProductPair signature;
        private final String smirks;
        private final String moiety;
        private ReactantProductPair smarts1;
        private ReactantProductPair signature1;
        private String smirks1;
        private ReactantProductPair smarts2;
        private ReactantProductPair signature2;
        private String smirks2;
        private ReactantProductPair smarts3;
        private ReactantProductPair signature3;
        private String smirks3;

        /**
         *
         * @param name
         * @param smarts
         * @param signature
         * @param smirks
         * @param moiety
         */
        public MoleculeMoleculePair(
                ReactantProductPair name,
                ReactantProductPair smarts,
                ReactantProductPair signature,
                String smirks,
                String moiety) {
            this.name = name;
            this.smarts = smarts;
            this.signature = signature;
            this.smirks = smirks;
            this.moiety = moiety;
        }

        /**
         * @return the name
         */
        public ReactantProductPair getName() {
            return name;
        }

        /**
         * @return the smarts
         */
        public ReactantProductPair getSmarts() {
            return smarts;
        }

        /**
         * @return the signature
         */
        public ReactantProductPair getSignature() {
            return signature;
        }

        /**
         * @return the smirks
         */
        public String getSmirks() {
            return smirks;
        }

        /**
         * @return the smarts at level 1
         */
        public ReactantProductPair getSmarts1() {
            return smarts1;
        }

        /**
         * @param smarts1 the smarts at level 1 to set
         */
        public void setSmarts1(ReactantProductPair smarts1) {
            this.smarts1 = smarts1;
        }

        /**
         * @return the signature at level 1
         */
        public ReactantProductPair getSignature1() {
            return signature1;
        }

        /**
         * @param signature1 the signature1 to set
         */
        public void setSignature1(ReactantProductPair signature1) {
            this.signature1 = signature1;
        }

        /**
         * @return the smirks at level 1
         */
        public String getSmirks1() {
            return smirks1;
        }

        /**
         * @param smirks1 the smirks at level 1 to set
         */
        public void setSmirks1(String smirks1) {
            this.smirks1 = smirks1;
        }

        /**
         * @return the smarts at level 2
         */
        public ReactantProductPair getSmarts2() {
            return smarts2;
        }

        /**
         * @param smarts2 the smarts at level 2 to set
         */
        public void setSmarts2(ReactantProductPair smarts2) {
            this.smarts2 = smarts2;
        }

        /**
         * @return the signature at level 2
         */
        public ReactantProductPair getSignature2() {
            return signature2;
        }

        /**
         * @param signature2 the signature at level 2 to set
         */
        public void setSignature2(ReactantProductPair signature2) {
            this.signature2 = signature2;
        }

        /**
         * @return the smirks at level 2
         */
        public String getSmirks2() {
            return smirks2;
        }

        /**
         * @param smirks2 the smirks at level 2 to set
         */
        public void setSmirks2(String smirks2) {
            this.smirks2 = smirks2;
        }

        /**
         * @return the smarts at level 3
         */
        public ReactantProductPair getSmarts3() {
            return smarts3;
        }

        /**
         * @param smarts3 the smarts3 to set
         */
        public void setSmarts3(ReactantProductPair smarts3) {
            this.smarts3 = smarts3;
        }

        /**
         * @return the signature at level 3
         */
        public ReactantProductPair getSignature3() {
            return signature3;
        }

        /**
         * @param signature3 the signature at level 3 to set
         */
        public void setSignature3(ReactantProductPair signature3) {
            this.signature3 = signature3;
        }

        /**
         * @return the smirks at level 3
         */
        public String getSmirks3() {
            return smirks3;
        }

        /**
         * @param smirks3 the smirks at level 3 to set
         */
        public void setSmirks3(String smirks3) {
            this.smirks3 = smirks3;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Name:").append(this.name);
            sb.append("signature:").append(this.signature);
            sb.append("smarts:").append(this.smarts);
            sb.append("smirks:").append(this.smirks);
            return sb.toString();
        }

        @Override
        public int compareTo(MoleculeMoleculePair o) {
            String local = this.name + this.smirks;
            String other = o.getName() + o.getSmirks();
            return local.compareTo(other);
        }

        @Override
        public int compare(MoleculeMoleculePair o1, MoleculeMoleculePair o2) {
            return o1.compareTo(o2);
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 37 * hash + Objects.hashCode(this.name);
            hash = 37 * hash + Objects.hashCode(this.smarts);
            hash = 37 * hash + Objects.hashCode(this.signature);
            hash = 37 * hash + Objects.hashCode(this.smirks);
            return hash;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final MoleculeMoleculePair other = (MoleculeMoleculePair) obj;
            if (!Objects.equals(this.name, other.name)) {
                return false;
            }
            if (!Objects.equals(this.smarts, other.smarts)) {
                return false;
            }
            if (!Objects.equals(this.signature, other.signature)) {
                return false;
            }
            return Objects.equals(this.smirks, other.smirks);
        }

        /**
         * @return the moiety
         */
        public String getMoiety() {
            return moiety;
        }
    }


    // ========== ReactantProductPair ==========

    public static class ReactantProductPair implements Serializable,
            Comparable<ReactantProductPair>, Comparator<ReactantProductPair> {

        private static final long serialVersionUID = 19876565735478L;
        private final String query;
        private final String target;

        /**
         *
         * @param query
         * @param target
         */
        public ReactantProductPair(String query, String target) {
            this.query = query;
            this.target = target;
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 97 * hash + Objects.hashCode(this.query);
            hash = 97 * hash + Objects.hashCode(this.target);
            return hash;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final ReactantProductPair other = (ReactantProductPair) obj;
            if (!Objects.equals(this.query, other.query)) {
                return false;
            }
            return Objects.equals(this.target, other.target);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("R:").append(this.getQuery());
            sb.append(", P:").append(this.getTarget());
            return sb.toString();
        }

        /**
         * @return the query
         */
        public String getQuery() {
            return query;
        }

        /**
         * @return the target
         */
        public String getTarget() {
            return target;
        }

        @Override
        public int compareTo(ReactantProductPair o) {
            String local = this.query + this.target;
            String object = o.getQuery() + o.getTarget();
            return local.compareTo(object);
        }

        @Override
        public int compare(ReactantProductPair o1, ReactantProductPair o2) {
            return o1.compareTo(o2);
        }
    }


    // ========== ReactionCenterFragment ==========

    public static class ReactionCenterFragment implements Serializable {

        private static final long serialVersionUID = 9879878799977781L;

        private final String signature;
        private final int level;
        private final BondChangeCalculator.EnumSubstrateProduct rpf;

        /**
         *
         * @param signature
         * @param level
         * @param rpf
         */
        public ReactionCenterFragment(String signature, int level, BondChangeCalculator.EnumSubstrateProduct rpf) {
            this.signature = signature;
            this.level = level;
            this.rpf = rpf;
        }

        @Override
        public String toString() {
            return "ReactionCenterFragment{" + "signature=" + signature + ", level=" + level + ", rpf=" + rpf + '}';
        }

        /**
         *
         * @return
         */
        public int getLevel() {
            return level;
        }

        /**
         *
         * @return
         */
        public BondChangeCalculator.EnumSubstrateProduct getReactantProductInfo() {
            return rpf;
        }

        /**
         *
         * @return
         */
        public String getSignature() {
            return signature;
        }
    }

}
