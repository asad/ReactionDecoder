/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping.algorithm;

import java.io.IOException;
import java.io.Serializable;
import static java.lang.String.valueOf;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Calendar;
import static java.util.Calendar.DATE;
import static java.util.Calendar.HOUR;
import static java.util.Calendar.MILLISECOND;
import static java.util.Calendar.MINUTE;
import static java.util.Calendar.MONTH;
import static java.util.Calendar.YEAR;
import java.util.Collection;
import static java.util.Collections.sort;
import static java.util.Collections.unmodifiableList;
import java.util.Comparator;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.MoleculeInitializer;
import org.openscience.smsd.BaseMapping.Algorithm;
import org.openscience.smsd.ExtAtomContainerManipulator;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import static com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.Similarity.getTanimotoSimilarity;
import com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.FingerprintGenerator;
import static com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.FingerprintGenerator.getFingerprinterSize;
import com.bioinceptionlabs.reactionblast.mapping.ThreadSafeCache;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.ReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.BestMatchContainer;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.HydrogenFreeFingerPrintContainer;
import com.bioinceptionlabs.reactionblast.mapping.MappingDiagnostics;
import com.bioinceptionlabs.reactionblast.mapping.MappingKeyUtil;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.MoleculeMoleculeMapping;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.MolMapping;
import com.bioinceptionlabs.reactionblast.mapping.SmsdReactionMappingEngine;
import static com.bioinceptionlabs.reactionblast.mapping.GraphMatcher.matcher;
import com.bioinceptionlabs.reactionblast.mapping.GraphMatcher.MCSSolution;
import com.bioinceptionlabs.reactionblast.mapping.GraphMatcher.GraphMatching;
import com.bioinceptionlabs.reactionblast.mapping.Reactor.Debugger;
import com.bioinceptionlabs.reactionblast.mapping.GraphMatcher.AbstractGraphMatching;
import com.bioinceptionlabs.reactionblast.mapping.BestMatch;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import com.bioinceptionlabs.reactionblast.tools.ICanonicalMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.tools.SmilesMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.tools.MoleculeTools.AtomContainerSetComparator;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
interface IGameTheory {

    int getDelta();

    MoleculeMoleculeMapping getReactionMolMapping();

    String getSuffix() throws IOException;

    void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping);

    void UpdateMatrix(Holder mh, boolean removeHydrogen) throws Exception;

    void UpdateMatrix(Collection<MCSSolution> mcsSolutions, Holder mh, boolean removeHydrogen) throws Exception;
}

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
interface IGraphTheoryMatrix {

    void Clear() throws IOException;

    int getDelta();

    List<String> getEductCounter();

    Holder getMatrixHolder();

    List<String> getProductCounter();

    MoleculeMoleculeMapping getReactionMolMapping();

    void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping);
}

/**
 * Consolidated game-theory engine for atom-atom mapping.
 * Merges: BaseGameTheory, GameTheoryFactory, GameTheoryMax, GameTheoryMin,
 *         GameTheoryMixture, GameTheoryRings, GameTheoryMatrix
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public abstract class GameTheoryEngine extends Debugger implements IGameTheory, Serializable {

    private final static ILoggingTool LOGGER
            = createLoggingTool(GameTheoryEngine.class);
    private static final long serialVersionUID = 1698688633678282L;
    private static final ReactionMappingEngine MAPPING_ENGINE
            = SmsdReactionMappingEngine.getInstance();

    // ---- BaseGameTheory methods inlined into outer class ----

    protected static boolean isPseudoAtoms(IAtomContainer atomContainer) {
        for (IAtom atoms : atomContainer.atoms()) {
            if (atoms instanceof IPseudoAtom || atoms instanceof PseudoAtom) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String getSuffix() throws IOException {
        Calendar cal = new GregorianCalendar();
        int ms = cal.get(YEAR);
        String suffix = valueOf(ms);
        ms = cal.get(MONTH);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(DATE);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(HOUR);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(MINUTE);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(MILLISECOND);
        suffix = suffix.concat(valueOf(ms));
        return suffix;
    }

    @Override
    public void UpdateMatrix(Holder mh, boolean removeHydrogen) throws InterruptedException {
        try {
            LOGGER.debug("**********Updated Matrix And Calculate Similarity**************");
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Collection<MCSSolution> mcsSolutions = null;
            try {
                mcsSolutions = matcher(mh);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw e;
            }
            Map<ReactionContainer.Key, MCSSolution> indexedSolutions = indexSolutions(mcsSolutions);
            for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
                for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
                    try {
                        IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                        IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                        LOGGER.debug("===================Mapped==================");
                        LOGGER.debug(" educt " + educt.getID() + ", product " + product.getID());
                        LOGGER.debug(" educt Index " + substrateIndex + ", product Index " + productIndex);
                        LOGGER.debug(" educt count " + educt.getAtomCount() + ", product count " + product.getAtomCount());
                        LOGGER.debug("mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) "
                                + mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex));
                        LOGGER.debug("reactionStructureInformation.isEductModified(substrateIndex) "
                                + reactionStructureInformation.isEductModified(substrateIndex));
                        LOGGER.debug("reactionStructureInformation.isProductModified(productIndex) "
                                + reactionStructureInformation.isProductModified(productIndex));
                        if ((educt != null && product != null)
                                && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                                && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                                || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
                            if (reactionStructureInformation.isEductModified(substrateIndex)
                                    || reactionStructureInformation.isProductModified(productIndex)) {
                                refillMatrixWithNewData(mh, substrateIndex, productIndex, indexedSolutions);
                            } else {
                                refillMatrixWithOldData(mh, substrateIndex, productIndex);
                            }
                        } else {
                            mh.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getStereoMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getCliqueMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                        }
                    } catch (IOException | CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                }
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw e;
        } catch (Exception e) {
            LOGGER.error("Error in matching molecules, check Graph Matcher module! ", e.getMessage());
        }
        try {
            resetFLAGS(mh);
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    @Override
    public void UpdateMatrix(Collection<MCSSolution> mcsSolutions, Holder mh, boolean removeHydrogen) throws Exception {
        try {
            LOGGER.debug("**********Updated Matrix And Calculate Similarity**************");
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Map<ReactionContainer.Key, MCSSolution> indexedSolutions = indexSolutions(mcsSolutions);
            for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
                for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
                    IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                    IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                    if ((educt != null && product != null)
                            && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                            && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                            || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
                        if (reactionStructureInformation.isEductModified(substrateIndex)
                                || reactionStructureInformation.isProductModified(productIndex)) {
                            refillMatrixWithNewData(mh, substrateIndex, productIndex, indexedSolutions);
                        } else {
                            refillMatrixWithOldData(mh, substrateIndex, productIndex);
                        }
                    } else {
                        mh.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getStereoMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getCliqueMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                    }
                }
            }
            resetFLAGS(mh);
        } catch (Exception e) {
            LOGGER.error("Error in matching molecules, check Graph Matcher module! ", e.getMessage());
        }
    }

    private void refillMatrixWithNewData(
            Holder holder, int substrateIndex, int productIndex,
            Map<ReactionContainer.Key, MCSSolution> solutionIndex) {
        LOGGER.debug("**********Generate NEW MCS And Calculate Similarity**************");
        try {
            ReactionContainer reactionContainer = holder.getReactionContainer();
            BestMatch initMcsAtom = holder.getBestMatchContainer();
            double stereoVal = 0.0;
            int fragmentVal = 0;
            double energyVal = 0.0;
            double graphSimilarity = 0.0;
            double mappingSize = 0.0;
            double fpSim = 0.0;
            double carbonCount = 0.0;
            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);
            LOGGER.debug("--Get matches--");
            MCSSolution atomatomMapping = getMappings(
                    holder.getReactionID(),
                    holder.getTheory() == null ? "UNKNOWN" : holder.getTheory().name(),
                    substrateIndex, productIndex, educt, product, solutionIndex);
            if (atomatomMapping == null) {
                clearScores(holder, substrateIndex, productIndex);
                return;
            }
            carbonCount = countMappedCarbons(atomatomMapping.getAtomAtomMapping());
            if (atomatomMapping.getStereoScore() != null) {
                stereoVal = atomatomMapping.getStereoScore();
            }
            if (atomatomMapping.getFragmentSize() != null) {
                fragmentVal = atomatomMapping.getFragmentSize();
            }
            if (atomatomMapping.getEnergy() != null) {
                energyVal = atomatomMapping.getEnergy();
            }
            AtomAtomMapping fam = atomatomMapping.getAtomAtomMapping();
            initMcsAtom.putBestMapping(substrateIndex, productIndex, fam);
            double ACount = educt.getAtomCount();
            double BCount = product.getAtomCount();
            mappingSize = atomatomMapping.getAtomAtomMapping().getCount();
            graphSimilarity = mappingSize / (ACount + BCount - mappingSize);
            initMcsAtom.setTotalFragmentCount(substrateIndex, productIndex, fragmentVal);
            initMcsAtom.setBondEnergy(substrateIndex, productIndex, energyVal);
            initMcsAtom.setStereoScore(substrateIndex, productIndex, stereoVal);
            initMcsAtom.setGraphSimilarity(substrateIndex, productIndex, graphSimilarity);
            BitSet a = reactionContainer.getFingerPrintofEduct(substrateIndex);
            BitSet b = reactionContainer.getFingerPrintofProduct(productIndex);
            if (a != null && b != null) {
                try {
                    fpSim = getTanimotoSimilarity(a, b);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, " error in calculating fingerprint ", ex.getMessage());
                }
            }
            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, carbonCount);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);
        } catch (IOException | CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
            clearScores(holder, substrateIndex, productIndex);
        }
    }

    private MCSSolution getMappings(
            String reactionId, String algorithmName,
            int queryPosition, int targetPosition,
            IAtomContainer educt, IAtomContainer product,
            Map<ReactionContainer.Key, MCSSolution> solutionIndex) throws CDKException {
        if (solutionIndex == null || solutionIndex.isEmpty()) {
            return quickMapping(reactionId, algorithmName, educt, product, queryPosition, targetPosition);
        }
        MCSSolution solution = solutionIndex.get(new ReactionContainer.Key(queryPosition, targetPosition));
        if (solution == null) {
            return null;
        }
        if (solution.getAtomAtomMapping().isEmpty()) {
            Set<String> atomMaps = new HashSet<>();
            for (IAtom a : educt.atoms()) {
                atomMaps.add(a.getSymbol());
            }
            boolean mappingPossible = false;
            for (IAtom a : product.atoms()) {
                if (atomMaps.contains(a.getSymbol())) {
                    mappingPossible = true;
                }
            }
            atomMaps.clear();
            if (mappingPossible) {
                return quickMapping(reactionId, algorithmName, educt, product, queryPosition, targetPosition);
            }
        }
        return solution;
    }

    private Map<ReactionContainer.Key, MCSSolution> indexSolutions(Collection<MCSSolution> mcsSolutions) {
        int initialCapacity = mcsSolutions == null ? 0 : Math.max(16, mcsSolutions.size() * 2);
        Map<ReactionContainer.Key, MCSSolution> indexedSolutions = new HashMap<>(initialCapacity);
        if (mcsSolutions == null) {
            return indexedSolutions;
        }
        for (MCSSolution solution : mcsSolutions) {
            if (solution == null) {
                continue;
            }
            indexedSolutions.put(
                    new ReactionContainer.Key(solution.getQueryPosition(), solution.getTargetPosition()),
                    solution);
        }
        return indexedSolutions;
    }

    private MCSSolution quickMapping(String reactionId, String algorithmName,
            IAtomContainer educt, IAtomContainer product,
            int queryPosition, int targetPosition) {
        ThreadSafeCache<String, MCSSolution> mappingcache = ThreadSafeCache.getInstance();
        LOGGER.debug("====Quick Mapping====");
        MappingDiagnostics.recordQuickMappingCall(reactionId, algorithmName);
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(educt);
            MoleculeInitializer.initializeMolecule(educt);
        } catch (CDKException ex) {
            LOGGER.error(Level.SEVERE, "Error in config. mol ", ex.getMessage());
        }
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(product);
            MoleculeInitializer.initializeMolecule(product);
        } catch (CDKException ex) {
            LOGGER.error(Level.SEVERE, "Error in config. mol ", ex.getMessage());
        }
        try {
            CycleFinder cycles = Cycles.vertexShort();
            Cycles rings = cycles.find(educt);
            int numberOfCyclesEduct = rings.numberOfCycles();
            rings = cycles.find(product);
            int numberOfCyclesProduct = rings.numberOfCycles();
            String key = generateUniqueKey(
                    educt, product,
                    educt.getID(), product.getID(),
                    educt.getAtomCount(), product.getAtomCount(),
                    educt.getBondCount(), product.getBondCount(),
                    false, false, false, false,
                    numberOfCyclesEduct, numberOfCyclesProduct);
            MCSSolution cached = mappingcache.get(key);
            if (cached != null) {
                MappingDiagnostics.recordQuickMappingCacheHit(reactionId, algorithmName);
                MCSSolution mcs = copyOldSolutionToNew(
                        queryPosition, targetPosition,
                        educt, product, cached);
                return mcs;
            } else {
                AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
                BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);
                MappingDiagnostics.recordQuickMappingSearch(reactionId, algorithmName);
                BaseMapping isomorphism = MAPPING_ENGINE.findMcs(
                        educt, product, Algorithm.DEFAULT, atomMatcher, bondMatcher);
                MCSSolution mcs = addMCSSolution(
                        queryPosition, targetPosition,
                        educt, product,
                        key, mappingcache, isomorphism);
                return mcs;
            }
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return null;
    }

    private void resetFLAGS(Holder mh) throws Exception {
        ReactionContainer reactionStructureInformation = mh.getReactionContainer();
        for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
            reactionStructureInformation.setEductModified(substrateIndex, false);
        }
        for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
            reactionStructureInformation.setProductModified(productIndex, false);
        }
    }

    protected final String canonicalMatchedSmiles(
            ICanonicalMoleculeLabeller canonLabeler,
            IAtomContainer matchedPart) throws Exception {
        IAtomContainer canonical = canonLabeler.getCanonicalMolecule(matchedPart);
        CDKSMILES cdkSmiles = new CDKSMILES(canonical, true, false);
        return cdkSmiles.getCanonicalSMILES();
    }

    private void refillMatrixWithOldData(Holder holder, int substrateIndex, int productIndex) {
        LOGGER.debug("**********REFILL MCS And Calculate Similarity**************");
        try {
            ReactionContainer reactionContainer = holder.getReactionContainer();
            BestMatch initMcsAtom = holder.getBestMatchContainer();
            double stereoVal = 0.0;
            int fragmentVal = 0;
            double energyVal = 0.0;
            double graphSimilarity = 0.0;
            double mappingSize = 0.0;
            double fpSim = 0.0;
            double carbonCount = 0.0;
            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);
            if (initMcsAtom.containsKey(substrateIndex, productIndex)) {
                AtomAtomMapping bestAtomAtomMapping = initMcsAtom.getAtomMatch(substrateIndex, productIndex);
                if (bestAtomAtomMapping == null) {
                    clearScores(holder, substrateIndex, productIndex);
                    return;
                }
                carbonCount = countMappedCarbons(bestAtomAtomMapping);
                stereoVal = initMcsAtom.getStereoScore(substrateIndex, productIndex);
                fragmentVal = initMcsAtom.getTotalFragmentCount(substrateIndex, productIndex);
                energyVal = initMcsAtom.getBondEnergy(substrateIndex, productIndex);
                double ACount = educt.getAtomCount();
                double BCount = product.getAtomCount();
                mappingSize = bestAtomAtomMapping.getCount();
                graphSimilarity = mappingSize / (ACount + BCount - mappingSize);
                initMcsAtom.setTotalFragmentCount(substrateIndex, productIndex, fragmentVal);
                initMcsAtom.setBondEnergy(substrateIndex, productIndex, energyVal);
                initMcsAtom.setStereoScore(substrateIndex, productIndex, stereoVal);
                initMcsAtom.setGraphSimilarity(substrateIndex, productIndex, graphSimilarity);
                BitSet a = reactionContainer.getFingerPrintofEduct(substrateIndex);
                BitSet b = reactionContainer.getFingerPrintofProduct(productIndex);
                if (a != null && b != null) {
                    try {
                        fpSim = getTanimotoSimilarity(a, b);
                    } catch (Exception ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                }
            }
            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, carbonCount);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);
        } catch (CDKException ex) {
            LOGGER.debug(SEVERE, null, ex);
            clearScores(holder, substrateIndex, productIndex);
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
            clearScores(holder, substrateIndex, productIndex);
        }
    }

    private void clearScores(Holder holder, int substrateIndex, int productIndex) {
        holder.getCliqueMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getStereoMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
        holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
    }

    private double countMappedCarbons(AtomAtomMapping mapping) {
        if (mapping == null) {
            return 0.0;
        }
        double carbonCount = 0.0;
        for (IAtom atom : mapping.getMappingsByAtoms().keySet()) {
            if ("C".equalsIgnoreCase(atom.getSymbol())) {
                carbonCount++;
            }
        }
        return carbonCount;
    }

    String generateUniqueKey(
            IAtomContainer compound1, IAtomContainer compound2,
            String id1, String id2,
            int atomCount1, int atomCount2,
            int bondCount1, int bondCount2,
            boolean atomtypeMatcher, boolean bondMatcher,
            boolean ringMatcher, boolean hasPerfectRings,
            int numberOfCyclesEduct, int numberOfCyclesProduct) {
        return MappingKeyUtil.buildPairKey(
                compound1,
                compound2,
                "QUICK",
                atomtypeMatcher,
                bondMatcher,
                ringMatcher,
                hasPerfectRings);
    }

    MCSSolution copyOldSolutionToNew(int queryPosition, int targetPosition,
            IAtomContainer compound1, IAtomContainer compound2, MCSSolution oldSolution) {
        AtomAtomMapping atomAtomMapping = oldSolution.getAtomAtomMapping();
        Map<Integer, Integer> mappingsByIndex = atomAtomMapping.getMappingsByIndex();
        AtomAtomMapping atomAtomMappingNew = new AtomAtomMapping(compound1, compound2);
        mappingsByIndex.entrySet().forEach((m) -> {
            atomAtomMappingNew.put(compound1.getAtom(m.getKey()), compound2.getAtom(m.getValue()));
        });
        MCSSolution mcsSolution = new MCSSolution(queryPosition, targetPosition, compound1, compound2, atomAtomMappingNew);
        mcsSolution.setEnergy(oldSolution.getEnergy());
        mcsSolution.setFragmentSize(oldSolution.getFragmentSize());
        mcsSolution.setStereoScore(oldSolution.getStereoScore());
        return mcsSolution;
    }

    MCSSolution addMCSSolution(int queryPosition, int targetPosition,
            IAtomContainer educt, IAtomContainer product,
            String key, ThreadSafeCache<String, MCSSolution> mappingcache, BaseMapping isomorphism) {
        MAPPING_ENGINE.applyDefaultFilters(isomorphism);
        MCSSolution mcs = new MCSSolution(queryPosition, targetPosition,
                isomorphism.getQuery(), isomorphism.getTarget(), isomorphism.getFirstAtomMapping());
        mcs.setEnergy(isomorphism.getEnergyScore(0));
        mcs.setFragmentSize(isomorphism.getFragmentSize(0));
        mcs.setStereoScore(isomorphism.getStereoScore(0));
        MCSSolution cached = mappingcache.putIfAbsent(key, mcs);
        if (cached == mcs) {
            return mcs;
        }
        return copyOldSolutionToNew(queryPosition, targetPosition, educt, product, cached);
    }

    // ========== Inner class: GameTheoryFactory ==========

    public static class GameTheoryFactory implements Serializable {

        private static final long serialVersionUID = 01567272317571L;

        public static IGameTheory make(IMappingAlgorithm theory, IReaction reaction, boolean removeHydrogen,
                Map<Integer, IAtomContainer> educts, Map<Integer, IAtomContainer> products,
                GameTheoryMatrix rpsh) throws Exception {
            switch (theory) {
                case MIXTURE:
                    return new GameTheoryMixture(reaction, removeHydrogen, educts, products, rpsh);
                case MIN:
                    return new GameTheoryMin(reaction, removeHydrogen, educts, products, rpsh);
                case MAX:
                    return new GameTheoryMax(reaction, removeHydrogen, educts, products, rpsh);
                case RINGS:
                    return new GameTheoryRings(reaction, removeHydrogen, educts, products, rpsh);
                default:
                    return null;
            }
        }

        private GameTheoryFactory() {
        }
    }

    // ========== Inner class: GameTheoryMax ==========

    static final class GameTheoryMax extends GameTheoryEngine {

        private static final int MAX_MAPPING_ITERATIONS = 100;
        private static final long serialVersionUID = 1887868678797L;
        private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryMax.class);
        private final List<String> eductList;
        private final List<String> productList;
        private final MappingChecks.ChooseWinner winner;
        private final IReaction reaction;
        private final String rid;
        private final String dirSuffix;
        private final boolean removeHydrogen;
        private MoleculeMoleculeMapping reactionMolMapping = null;
        private Map<Integer, IAtomContainer> educts = null;
        private Map<Integer, IAtomContainer> products = null;
        private Holder mh;
        private int delta = 0;
        private Integer stepIndex = 0;
        private final ICanonicalMoleculeLabeller canonLabeler;
        private int counter = 0;

        GameTheoryMax(IReaction reaction, boolean removeHydrogen,
                Map<Integer, IAtomContainer> _educts, Map<Integer, IAtomContainer> _products,
                GameTheoryMatrix rpsh) throws Exception {
            LOGGER.debug("I am MAX");
            this.counter = 1;
            this.canonLabeler = new SmilesMoleculeLabeller();
            this.removeHydrogen = removeHydrogen;
            this.reaction = reaction;
            this.educts = _educts;
            this.products = _products;
            this.rid = reaction.getID();
            this.eductList = new ArrayList<>(rpsh.getEductCounter());
            this.productList = new ArrayList<>(rpsh.getProductCounter());
            this.mh = rpsh.getMatrixHolder();
            setReactionMolMapping(rpsh.getReactionMolMapping());
            this.winner = new MappingChecks.ChooseWinner(eductList, productList);
            this.dirSuffix = super.getSuffix();
            GenerateMapping(false);
        }

        private void GenerateMapping(boolean flag) throws Exception {
            boolean ruleMatchingFlag = flag;
            int iteration = 0;
            boolean continueMapping = true;
            while (continueMapping && iteration < MAX_MAPPING_ITERATIONS) {
                if (Thread.interrupted()) {
                    throw new InterruptedException("MAX mapping interrupted at iteration " + iteration);
                }
                this.counter++;
                boolean conditionmet = false;
                if (!ruleMatchingFlag) {
                    MappingChecks.RuleBasedMappingHandler ruleBasedMappingHandler = new MappingChecks.RuleBasedMappingHandler(mh, eductList, productList);
                    if (ruleBasedMappingHandler.isMatchFound()) {
                        mh = MappingChecks.Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                        conditionmet = true;
                    }
                    ruleMatchingFlag = true;
                }
                if (!conditionmet && counter <= 5) {
                    MappingChecks.MaxSelection select = new MappingChecks.MaxSelection(mh, eductList, productList);
                    if (select.isSubAndCompleteMatchFlag()) {
                        mh = select.getUpdatedHolder();
                    }
                }
                winner.searchWinners(educts, products, mh);
                if (winner.getFlag()) {
                    UpdateMapping();
                    ReactionContainer rc = mh.getReactionContainer();
                    boolean allMapped = true;
                    for (int i = 0; i < rc.getEductCount() && allMapped; i++) {
                        if (rc.getEduct(i).getAtomCount() > 0) allMapped = false;
                    }
                    for (int j = 0; j < rc.getProductCount() && allMapped; j++) {
                        if (rc.getProduct(j).getAtomCount() > 0) allMapped = false;
                    }
                    if (allMapped) { break; }
                    boolean hasRemainingPairs = false;
                    for (int i = 0; i < mh.getGraphSimilarityMatrix().getRowDimension(); i++) {
                        for (int j = 0; j < mh.getGraphSimilarityMatrix().getColumnDimension(); j++) {
                            if (mh.getGraphSimilarityMatrix().getValue(i, j) > 0) {
                                hasRemainingPairs = true;
                                break;
                            }
                        }
                        if (hasRemainingPairs) break;
                    }
                    if (!hasRemainingPairs) { break; }
                    UpdateMatrix(mh, removeHydrogen);
                    iteration++;
                } else {
                    continueMapping = false;
                }
            }
        }

        private void UpdateMapping() throws Exception {
            boolean[][] FlagMatrix = winner.getFlagMatrix();
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            for (int iIndex = 0; iIndex < reactionStructureInformation.getEductCount(); iIndex++) {
                for (int jIndex = 0; jIndex < reactionStructureInformation.getProductCount(); jIndex++) {
                    int substrateIndex = iIndex;
                    int productIndex = jIndex;
                    IAtomContainer ac1 = reactionStructureInformation.getEduct(substrateIndex);
                    IAtomContainer ac2 = reactionStructureInformation.getProduct(productIndex);
                    if (FlagMatrix[substrateIndex][productIndex]) {
                        BitSet A = reactionStructureInformation.getFingerPrintofEduct(substrateIndex);
                        BitSet B = reactionStructureInformation.getFingerPrintofProduct(productIndex);
                        ac1.setID(this.eductList.get(substrateIndex));
                        ac2.setID(this.productList.get(productIndex));
                        AbstractGraphMatching GM = new GraphMatching(rid, ac1, ac2, dirSuffix, removeHydrogen);
                        boolean mcsMatch = GM.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, A, B);
                        if (mcsMatch) {
                            delta += GM.removeMatchedAtomsAndUpdateAAM(reaction);
                            List<MolMapping> rMap = getReactionMolMapping().
                                    getMapping(rid, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                            String matchedSmiles = canonicalMatchedSmiles(canonLabeler, GM.getMatchedPart());
                            for (MolMapping map : rMap) {
                                map.setReactionMapping(true);
                                map.setMatchedSMILES(matchedSmiles, ++stepIndex);
                            }
                        }
                        IAtomContainer remainingEduct = GM.getRemainingEduct();
                        IAtomContainer remainingProduct = GM.getRemainingProduct();
                        reactionStructureInformation.putEduct(substrateIndex, remainingEduct);
                        reactionStructureInformation.putProduct(productIndex, remainingProduct);
                        reactionStructureInformation.setEductModified(substrateIndex, true);
                        reactionStructureInformation.setProductModified(productIndex, true);
                    }
                }
            }
        }

        @Override
        public MoleculeMoleculeMapping getReactionMolMapping() { return reactionMolMapping; }
        @Override
        public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) { this.reactionMolMapping = reactionMolMapping; }
        @Override
        public int getDelta() { return delta; }
    }

    // ========== Inner class: GameTheoryMin ==========

    static final class GameTheoryMin extends GameTheoryEngine {

        private static final int MAX_MAPPING_ITERATIONS = 100;
        private static final long serialVersionUID = 1808979786969868698L;
        private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryMin.class);
        private final List<String> eductList;
        private final List<String> productList;
        private Holder mh;
        private final MappingChecks.ChooseWinner winner;
        private final IReaction reaction;
        private final String reactionName;
        private final String _dirSuffix;
        private final boolean removeHydrogen;
        private MoleculeMoleculeMapping reactionMolMapping = null;
        private Map<Integer, IAtomContainer> educts = null;
        private Map<Integer, IAtomContainer> products = null;
        private int delta = 0;
        private Integer stepIndex = 0;
        private final ICanonicalMoleculeLabeller canonLabeler;
        private int counter = 0;

        GameTheoryMin(IReaction reaction, boolean removeHydrogen,
                Map<Integer, IAtomContainer> _educts, Map<Integer, IAtomContainer> _products,
                GameTheoryMatrix rpsh) throws Exception {
            LOGGER.debug("I am MIN MIX");
            this.counter = 1;
            this.canonLabeler = new SmilesMoleculeLabeller();
            this.removeHydrogen = removeHydrogen;
            this.reaction = reaction;
            this.educts = _educts;
            this.products = _products;
            this.reactionName = reaction.getID();
            this.eductList = new ArrayList<>(rpsh.getEductCounter());
            this.productList = new ArrayList<>(rpsh.getProductCounter());
            this.mh = rpsh.getMatrixHolder();
            setReactionMolMapping(rpsh.getReactionMolMapping());
            winner = new MappingChecks.ChooseWinner(eductList, productList);
            this._dirSuffix = super.getSuffix();
            MappingChecks.ReactionIsomorphismHandler RIH = new MappingChecks.ReactionIsomorphismHandler(mh, eductList, productList);
            if (RIH.getIsomorphismFlag()) {
                mh = RIH.getMatrixHolder();
                GenerateIsoMorphismMapping();
            } else {
                GenerateMapping(false);
            }
        }

        private void GenerateIsoMorphismMapping() throws Exception {
            winner.searchWinners(educts, products, mh);
            if (winner.getFlag()) {
                UpdateMapping();
                UpdateMatrix(mh, removeHydrogen);
                GenerateMapping(false);
            }
        }

        private void GenerateMapping(boolean flag) throws Exception {
            boolean ruleMatchingFlag = flag;
            int iteration = 0;
            boolean continueMapping = true;
            while (continueMapping && iteration < MAX_MAPPING_ITERATIONS) {
                if (Thread.interrupted()) {
                    throw new InterruptedException("MIN mapping interrupted at iteration " + iteration);
                }
                this.counter++;
                boolean conditionmet = false;
                if (!ruleMatchingFlag) {
                    MappingChecks.RuleBasedMappingHandler ruleBasedMappingHandler = new MappingChecks.RuleBasedMappingHandler(mh, eductList, productList);
                    if (ruleBasedMappingHandler.isMatchFound()) {
                        mh = MappingChecks.Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                        conditionmet = true;
                    }
                    ruleMatchingFlag = true;
                }
                if (!conditionmet && counter <= 5) {
                    MappingChecks.MinSelection select = new MappingChecks.MinSelection(mh, eductList, productList);
                    if (select.isSubAndCompleteMatchFlag()) {
                        mh = select.getUpdatedHolder();
                    }
                }
                winner.searchWinners(educts, products, mh);
                if (winner.getFlag()) {
                    UpdateMapping();
                    ReactionContainer rc = mh.getReactionContainer();
                    boolean allMapped = true;
                    for (int i = 0; i < rc.getEductCount() && allMapped; i++) {
                        if (rc.getEduct(i).getAtomCount() > 0) allMapped = false;
                    }
                    for (int j = 0; j < rc.getProductCount() && allMapped; j++) {
                        if (rc.getProduct(j).getAtomCount() > 0) allMapped = false;
                    }
                    if (allMapped) { break; }
                    boolean hasRemainingPairs = false;
                    for (int i = 0; i < mh.getGraphSimilarityMatrix().getRowDimension(); i++) {
                        for (int j = 0; j < mh.getGraphSimilarityMatrix().getColumnDimension(); j++) {
                            if (mh.getGraphSimilarityMatrix().getValue(i, j) > 0) {
                                hasRemainingPairs = true;
                                break;
                            }
                        }
                        if (hasRemainingPairs) break;
                    }
                    if (!hasRemainingPairs) { break; }
                    UpdateMatrix(mh, removeHydrogen);
                    iteration++;
                } else {
                    continueMapping = false;
                }
            }
        }

        private void UpdateMapping() throws Exception {
            boolean[][] FlagMatrix = winner.getFlagMatrix();
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            for (int iIndex = 0; iIndex < reactionStructureInformationContainer.getEductCount(); iIndex++) {
                for (int jIndex = 0; jIndex < reactionStructureInformationContainer.getProductCount(); jIndex++) {
                    int substrateIndex = iIndex;
                    int productIndex = jIndex;
                    IAtomContainer ac1 = reactionStructureInformationContainer.getEduct(substrateIndex);
                    IAtomContainer ac2 = reactionStructureInformationContainer.getProduct(productIndex);
                    if (FlagMatrix[substrateIndex][productIndex]) {
                        BitSet a_BitSet = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);
                        BitSet b_BitSet = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);
                        ac1.setID(this.eductList.get(substrateIndex));
                        ac2.setID(this.productList.get(productIndex));
                        AbstractGraphMatching graphMatching = new GraphMatching(reactionName, ac1, ac2, _dirSuffix, removeHydrogen);
                        boolean mcsMatch = graphMatching.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, a_BitSet, b_BitSet);
                        if (mcsMatch) {
                            delta += graphMatching.removeMatchedAtomsAndUpdateAAM(reaction);
                            List<MolMapping> rMap = getReactionMolMapping().
                                    getMapping(reactionName, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                            String matchedSmiles = canonicalMatchedSmiles(canonLabeler, graphMatching.getMatchedPart());
                            for (MolMapping map : rMap) {
                                map.setReactionMapping(true);
                                map.setMatchedSMILES(matchedSmiles, ++stepIndex);
                            }
                        }
                        IAtomContainer remainingEduct = graphMatching.getRemainingEduct();
                        IAtomContainer remainingProduct = graphMatching.getRemainingProduct();
                        reactionStructureInformationContainer.putEduct(substrateIndex, remainingEduct);
                        reactionStructureInformationContainer.putProduct(productIndex, remainingProduct);
                        reactionStructureInformationContainer.setEductModified(substrateIndex, true);
                        reactionStructureInformationContainer.setProductModified(productIndex, true);
                    }
                }
            }
        }

        @Override
        public MoleculeMoleculeMapping getReactionMolMapping() { return reactionMolMapping; }
        @Override
        public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) { this.reactionMolMapping = reactionMolMapping; }
        @Override
        public int getDelta() { return delta; }
    }

    // ========== Inner class: GameTheoryMixture ==========

    static final class GameTheoryMixture extends GameTheoryEngine {

        private static final int MAX_MAPPING_ITERATIONS = 100;
        private static final long serialVersionUID = 1808979786969868698L;
        private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryMixture.class);
        private final List<String> eductList;
        private final List<String> productList;
        private Holder mh;
        private final MappingChecks.ChooseWinner winner;
        private final IReaction reaction;
        private final String RID;
        private final String _dirSuffix;
        private final boolean removeHydrogen;
        private MoleculeMoleculeMapping reactionMolMapping = null;
        private Map<Integer, IAtomContainer> educts = null;
        private Map<Integer, IAtomContainer> products = null;
        private int delta = 0;
        private Integer stepIndex = 0;
        private final ICanonicalMoleculeLabeller canonLabeler;

        GameTheoryMixture(IReaction reaction, boolean removeHydrogen,
                Map<Integer, IAtomContainer> _educts, Map<Integer, IAtomContainer> _products,
                GameTheoryMatrix rpsh) throws Exception {
            LOGGER.debug("I am MIXTURE");
            this.canonLabeler = new SmilesMoleculeLabeller();
            this.removeHydrogen = removeHydrogen;
            this.reaction = reaction;
            this.educts = _educts;
            this.products = _products;
            this.RID = reaction.getID();
            this.eductList = new ArrayList<>(rpsh.getEductCounter());
            this.productList = new ArrayList<>(rpsh.getProductCounter());
            this.mh = rpsh.getMatrixHolder();
            setReactionMolMapping(rpsh.getReactionMolMapping());
            winner = new MappingChecks.ChooseWinner(eductList, productList);
            this._dirSuffix = super.getSuffix();
            MappingChecks.ReactionIsomorphismHandler RIH = new MappingChecks.ReactionIsomorphismHandler(mh, eductList, productList);
            if (RIH.getIsomorphismFlag()) {
                mh = RIH.getMatrixHolder();
                GenerateIsoMorphismMapping();
            } else {
                GenerateMapping(false);
            }
        }

        private void GenerateIsoMorphismMapping() throws Exception {
            winner.searchWinners(educts, products, mh);
            if (winner.getFlag()) {
                UpdateMapping();
                UpdateMatrix(mh, removeHydrogen);
                GenerateMapping(false);
            }
        }

        private void GenerateMapping(boolean flag) throws Exception {
            boolean ruleMatchingFlag = flag;
            int iteration = 0;
            boolean continueMapping = true;
            while (continueMapping && iteration < MAX_MAPPING_ITERATIONS) {
                if (Thread.interrupted()) {
                    throw new InterruptedException("RINGS mapping interrupted at iteration " + iteration);
                }
                if (!ruleMatchingFlag) {
                    MappingChecks.RuleBasedMappingHandler ruleBasedMappingHandler
                            = new MappingChecks.RuleBasedMappingHandler(mh, eductList, productList);
                    if (ruleBasedMappingHandler.isMatchFound()) {
                        mh = MappingChecks.Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                    }
                    ruleMatchingFlag = true;
                }
                winner.searchWinners(educts, products, mh);
                if (winner.getFlag()) {
                    UpdateMapping();
                    ReactionContainer rc = mh.getReactionContainer();
                    boolean allMapped = true;
                    for (int i = 0; i < rc.getEductCount() && allMapped; i++) {
                        if (rc.getEduct(i).getAtomCount() > 0) allMapped = false;
                    }
                    for (int j = 0; j < rc.getProductCount() && allMapped; j++) {
                        if (rc.getProduct(j).getAtomCount() > 0) allMapped = false;
                    }
                    if (allMapped) { break; }
                    boolean hasRemainingPairs = false;
                    for (int i = 0; i < mh.getGraphSimilarityMatrix().getRowDimension(); i++) {
                        for (int j = 0; j < mh.getGraphSimilarityMatrix().getColumnDimension(); j++) {
                            if (mh.getGraphSimilarityMatrix().getValue(i, j) > 0) {
                                hasRemainingPairs = true;
                                break;
                            }
                        }
                        if (hasRemainingPairs) break;
                    }
                    if (!hasRemainingPairs) { break; }
                    UpdateMatrix(mh, removeHydrogen);
                    iteration++;
                } else {
                    continueMapping = false;
                }
            }
        }

        private void UpdateMapping() throws Exception {
            boolean[][] FlagMatrix = winner.getFlagMatrix();
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            for (int iIndex = 0; iIndex < reactionStructureInformationContainer.getEductCount(); iIndex++) {
                for (int jIndex = 0; jIndex < reactionStructureInformationContainer.getProductCount(); jIndex++) {
                    int substrateIndex = iIndex;
                    int productIndex = jIndex;
                    IAtomContainer ac1 = reactionStructureInformationContainer.getEduct(substrateIndex);
                    IAtomContainer ac2 = reactionStructureInformationContainer.getProduct(productIndex);
                    if (FlagMatrix[substrateIndex][productIndex]) {
                        BitSet A = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);
                        BitSet B = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);
                        ac1.setID(this.eductList.get(substrateIndex));
                        ac2.setID(this.productList.get(productIndex));
                        AbstractGraphMatching GM = new GraphMatching(RID, ac1, ac2, _dirSuffix, removeHydrogen);
                        boolean mcsMatch = GM.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, A, B);
                        if (mcsMatch) {
                            delta += GM.removeMatchedAtomsAndUpdateAAM(reaction);
                            List<MolMapping> rMap = getReactionMolMapping().
                                    getMapping(RID, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                            String matchedSmiles = canonicalMatchedSmiles(canonLabeler, GM.getMatchedPart());
                            for (MolMapping map : rMap) {
                                map.setReactionMapping(true);
                                map.setMatchedSMILES(matchedSmiles, ++stepIndex);
                            }
                        }
                        IAtomContainer RemainingEduct = GM.getRemainingEduct();
                        IAtomContainer RemainingProduct = GM.getRemainingProduct();
                        reactionStructureInformationContainer.putEduct(substrateIndex, RemainingEduct);
                        reactionStructureInformationContainer.putProduct(productIndex, RemainingProduct);
                        reactionStructureInformationContainer.setEductModified(substrateIndex, true);
                        reactionStructureInformationContainer.setProductModified(productIndex, true);
                    }
                }
            }
        }

        @Override
        public MoleculeMoleculeMapping getReactionMolMapping() { return reactionMolMapping; }
        @Override
        public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) { this.reactionMolMapping = reactionMolMapping; }
        @Override
        public int getDelta() { return delta; }
    }

    // ========== Inner class: GameTheoryRings ==========

    static final class GameTheoryRings extends GameTheoryEngine {

        private static final int MAX_MAPPING_ITERATIONS = 100;
        private static final long serialVersionUID = 0x152ec264bc2L;
        private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryRings.class);
        private final List<String> eductList;
        private final List<String> productList;
        private Holder mh;
        private final MappingChecks.ChooseWinner winner;
        private final IReaction reaction;
        private final String RID;
        private final String _dirSuffix;
        private final boolean removeHydrogen;
        private MoleculeMoleculeMapping reactionMolMapping = null;
        private Map<Integer, IAtomContainer> educts = null;
        private Map<Integer, IAtomContainer> products = null;
        private int delta = 0;
        private Integer stepIndex = 0;
        private final ICanonicalMoleculeLabeller canonLabeler;

        GameTheoryRings(IReaction reaction, boolean removeHydrogen,
                Map<Integer, IAtomContainer> _educts, Map<Integer, IAtomContainer> _products,
                GameTheoryMatrix rpsh) throws Exception {
            LOGGER.debug("I am Ring");
            this.canonLabeler = new SmilesMoleculeLabeller();
            this.removeHydrogen = removeHydrogen;
            this.reaction = reaction;
            this.educts = _educts;
            this.products = _products;
            this.RID = reaction.getID();
            this.eductList = new ArrayList<>(rpsh.getEductCounter());
            this.productList = new ArrayList<>(rpsh.getProductCounter());
            this.mh = rpsh.getMatrixHolder();
            setReactionMolMapping(rpsh.getReactionMolMapping());
            winner = new MappingChecks.ChooseWinner(eductList, productList);
            this._dirSuffix = super.getSuffix();
            MappingChecks.ReactionIsomorphismHandler RIH = new MappingChecks.ReactionIsomorphismHandler(mh, eductList, productList);
            if (RIH.getIsomorphismFlag()) {
                mh = RIH.getMatrixHolder();
                GenerateIsoMorphismMapping();
            } else {
                GenerateMapping();
            }
        }

        private void GenerateIsoMorphismMapping() throws Exception {
            MappingChecks.RuleBasedMappingHandler ph = new MappingChecks.RuleBasedMappingHandler(mh, eductList, productList);
            if (ph.isMatchFound()) {
                mh = ph.getMatrixHolder();
            }
            winner.searchWinners(educts, products, mh);
            if (winner.getFlag()) {
                UpdateMapping();
                UpdateMatrix(mh, removeHydrogen);
                GenerateMapping();
            }
        }

        private void GenerateMapping() throws Exception {
            int iteration = 0;
            boolean continueMapping = true;
            while (continueMapping && iteration < MAX_MAPPING_ITERATIONS) {
                if (Thread.interrupted()) {
                    throw new InterruptedException("MIXTURE mapping interrupted at iteration " + iteration);
                }
                MappingChecks.RuleBasedMappingHandler ruleBasedMappingHandler = new MappingChecks.RuleBasedMappingHandler(mh, eductList, productList);
                if (ruleBasedMappingHandler.isMatchFound()) {
                    mh = MappingChecks.Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                }
                winner.searchWinners(educts, products, mh);
                if (winner.getFlag()) {
                    UpdateMapping();
                    ReactionContainer rc = mh.getReactionContainer();
                    boolean allMapped = true;
                    for (int i = 0; i < rc.getEductCount() && allMapped; i++) {
                        if (rc.getEduct(i).getAtomCount() > 0) allMapped = false;
                    }
                    for (int j = 0; j < rc.getProductCount() && allMapped; j++) {
                        if (rc.getProduct(j).getAtomCount() > 0) allMapped = false;
                    }
                    if (allMapped) { break; }
                    boolean hasRemainingPairs = false;
                    for (int i = 0; i < mh.getGraphSimilarityMatrix().getRowDimension(); i++) {
                        for (int j = 0; j < mh.getGraphSimilarityMatrix().getColumnDimension(); j++) {
                            if (mh.getGraphSimilarityMatrix().getValue(i, j) > 0) {
                                hasRemainingPairs = true;
                                break;
                            }
                        }
                        if (hasRemainingPairs) break;
                    }
                    if (!hasRemainingPairs) { break; }
                    UpdateMatrix(mh, removeHydrogen);
                    iteration++;
                } else {
                    continueMapping = false;
                }
            }
        }

        private void UpdateMapping() throws Exception {
            boolean[][] FlagMatrix = winner.getFlagMatrix();
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            for (int iIndex = 0; iIndex < reactionStructureInformationContainer.getEductCount(); iIndex++) {
                for (int jIndex = 0; jIndex < reactionStructureInformationContainer.getProductCount(); jIndex++) {
                    int substrateIndex = iIndex;
                    int productIndex = jIndex;
                    IAtomContainer ac1 = reactionStructureInformationContainer.getEduct(substrateIndex);
                    IAtomContainer ac2 = reactionStructureInformationContainer.getProduct(productIndex);
                    if (FlagMatrix[substrateIndex][productIndex]) {
                        BitSet A = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);
                        BitSet B = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);
                        ac1.setID(this.eductList.get(substrateIndex));
                        ac2.setID(this.productList.get(productIndex));
                        AbstractGraphMatching GM = new GraphMatching(RID, ac1, ac2, _dirSuffix, removeHydrogen);
                        boolean mcsMatch = GM.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, A, B);
                        if (mcsMatch) {
                            delta += GM.removeMatchedAtomsAndUpdateAAM(reaction);
                            List<MolMapping> rMap = getReactionMolMapping().
                                    getMapping(RID, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                            String matchedSmiles = canonicalMatchedSmiles(canonLabeler, GM.getMatchedPart());
                            for (MolMapping map : rMap) {
                                map.setReactionMapping(true);
                                map.setMatchedSMILES(matchedSmiles, ++stepIndex);
                            }
                        }
                        IAtomContainer remainingEduct = GM.getRemainingEduct();
                        IAtomContainer remainingProduct = GM.getRemainingProduct();
                        reactionStructureInformationContainer.putEduct(substrateIndex, remainingEduct);
                        reactionStructureInformationContainer.putProduct(productIndex, remainingProduct);
                        reactionStructureInformationContainer.setEductModified(substrateIndex, true);
                        reactionStructureInformationContainer.setProductModified(productIndex, true);
                    }
                }
            }
        }

        @Override
        public MoleculeMoleculeMapping getReactionMolMapping() { return reactionMolMapping; }
        @Override
        public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) { this.reactionMolMapping = reactionMolMapping; }
        @Override
        public int getDelta() { return delta; }
    }

    // ========== Inner class: GameTheoryMatrix ==========

    public static class GameTheoryMatrix extends GameTheoryEngine implements IGraphTheoryMatrix {

        private static final long serialVersionUID = 0x2c36427fd2L;
        private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryMatrix.class);
        private Holder matrixHolder;
        private MoleculeMoleculeMapping reactionBlastMolMapping;
        private final List<String> eductCounter;
        private final List<String> productCounter;
        private final IReaction reaction;
        private final Map<Integer, BitSet> substrateductFPMap;
        private final Map<Integer, BitSet> productFPMap;
        private final FingerprintGenerator fpr;
        private final HydrogenFreeFingerPrintContainer hydFreeFPContainer;
        private final boolean removeHydrogen;
        private final String reactionID;
        private final ReactionContainer structureMapObj;
        private final BestMatch bestMatchContainer;
        private final IMappingAlgorithm theory;

        public GameTheoryMatrix(IMappingAlgorithm theory, IReaction reaction, boolean removeHydrogen) throws Exception {
            this.theory = theory;
            this.removeHydrogen = removeHydrogen;
            this.reaction = reaction;
            this.reactionID = reaction.getID();
            this.substrateductFPMap = new TreeMap<>();
            this.productFPMap = new TreeMap<>();
            this.fpr = new FingerprintGenerator();
            this.eductCounter = new LinkedList<>();
            this.productCounter = new LinkedList<>();
            this.structureMapObj = new ReactionContainer();
            this.bestMatchContainer = new BestMatchContainer();
            this.hydFreeFPContainer = new HydrogenFreeFingerPrintContainer();
            this.reactionBlastMolMapping = new MoleculeMoleculeMapping();
            try {
                StoichiometricCoefficientReplicator_Structure_FingerPrint_MapGenerator();
                BuildScoringMatrix();
            } catch (Exception e) {
                LOGGER.error(e);
            }
        }

        private void BuildScoringMatrix() throws Exception {
            try {
                matrixHolder = new Holder(theory, reactionID, eductCounter, productCounter,
                        structureMapObj, bestMatchContainer, hydFreeFPContainer);
                this.reactionBlastMolMapping.setMolMappings(reactionID, matrixHolder.getMappingMolPair());
                UpdateMatrix(matrixHolder, removeHydrogen);
            } catch (Exception e) {
                LOGGER.error(SEVERE, null, e);
            }
        }

        @Override
        public void Clear() throws IOException {
            structureMapObj.Clear();
            hydFreeFPContainer.Clear();
            bestMatchContainer.Clear();
            substrateductFPMap.clear();
            productFPMap.clear();
            eductCounter.clear();
            productCounter.clear();
            matrixHolder = null;
            reactionBlastMolMapping = null;
        }

        @Override
        public List<String> getEductCounter() { return unmodifiableList(eductCounter); }

        @Override
        public List<String> getProductCounter() { return unmodifiableList(productCounter); }

        private void StoichiometricCoefficientReplicator_Structure_FingerPrint_MapGenerator() {
            List<IAtomContainer> ac = new LinkedList<>();
            List<IAtomContainer> pd = new LinkedList<>();
            sortAtomContainer(ac, pd);
            for (int key = 0; key < ac.size(); key++) {
                try {
                    IAtomContainer mol = ac.get(key).clone();
                    String eductID = ac.get(key).getID() != null ? ac.get(key).getID().trim() : String.valueOf(key);
                    mol.setID(eductID);
                    BitSet FP;
                    if (hydFreeFPContainer.isKeyPresent(eductID)) {
                        FP = hydFreeFPContainer.getFingerPrint(eductID);
                    } else if (mol.getAtomCount() > 0) {
                        IAtomContainer tempMol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
                        FP = fpr.getFingerprint(tempMol);
                    } else {
                        FP = new BitSet(getFingerprinterSize());
                    }
                    hydFreeFPContainer.setValue(eductID, FP);
                    eductCounter.add(key, eductID);
                    structureMapObj.putEduct(key, mol);
                    structureMapObj.setEductModified(key, true);
                    substrateductFPMap.put(key, FP);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
            for (int key = 0; key < pd.size(); key++) {
                try {
                    IAtomContainer mol = pd.get(key).clone();
                    String productID = pd.get(key).getID() != null ? pd.get(key).getID().trim() : String.valueOf(key);
                    mol.setID(productID);
                    BitSet fingerPrint;
                    if (hydFreeFPContainer.isKeyPresent(productID)) {
                        fingerPrint = hydFreeFPContainer.getFingerPrint(productID);
                    } else if (mol.getAtomCount() > 0) {
                        IAtomContainer tempMol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
                        fingerPrint = fpr.getFingerprint(tempMol);
                    } else {
                        fingerPrint = new BitSet(getFingerprinterSize());
                    }
                    hydFreeFPContainer.setValue(productID, fingerPrint);
                    productCounter.add(key, productID);
                    structureMapObj.putProduct(key, mol);
                    structureMapObj.setProductModified(key, true);
                    productFPMap.put(key, fingerPrint);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }

        @Override
        public Holder getMatrixHolder() { return matrixHolder; }

        @Override
        public MoleculeMoleculeMapping getReactionMolMapping() { return reactionBlastMolMapping; }

        private void sortAtomContainer(List<IAtomContainer> ac, List<IAtomContainer> pd) {
            for (IAtomContainer e : reaction.getReactants().atomContainers()) { ac.add(e); }
            for (IAtomContainer p : reaction.getProducts().atomContainers()) { pd.add(p); }
            try {
                Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
                sort(ac, comparator);
                sort(pd, comparator);
            } catch (Exception e) {
                LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            }
        }

        @Override
        public int getDelta() { throw new UnsupportedOperationException("Not supported yet."); }

        @Override
        public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) {
            this.reactionBlastMolMapping = reactionMolMapping;
        }
    }
}
