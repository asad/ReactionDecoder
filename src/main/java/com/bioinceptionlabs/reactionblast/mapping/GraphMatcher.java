/*
 * GraphMatcher - consolidated graph matching classes.
 * Merged: Combination, GraphMatching, MCSSolution, MCSThread into GraphMatcher
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinception.smsd.core.SearchEngine;
import com.bioinceptionlabs.reactionblast.mapping.BestMatch;
import com.bioinceptionlabs.reactionblast.mapping.Reactor.Debugger;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.Key;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.ThreadSafeCache;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.BaseMapping.Algorithm;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.ExtAtomContainerManipulator;
import static java.lang.Runtime.getRuntime;
import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.getProperty;
import static java.lang.System.nanoTime;
import static java.util.Collections.unmodifiableCollection;
import static java.util.Collections.unmodifiableMap;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;
import static org.openscience.cdk.CDKConstants.UNSET;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.smsd.ExtAtomContainerManipulator.cloneWithIDs;


/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class GraphMatcher extends Debugger {

    private static final int SINGLE_SUBGRAPH_MATCH = 1;
    private static final long SUBGRAPH_TIMEOUT_MS = 5_000L;
    private static final long MCS_TIMEOUT_MS = 10_000L;
    /** Hard timeout per poll() call waiting for the next completed MCS pair. */
    private static final long MCS_POLL_TIMEOUT_MS = 15_000L;
    /** Overall wall-clock budget for the entire matcher() call. */
    private static final long MATCHER_BUDGET_MS = 60_000L;
    /** Hard timeout for the executor shutdown after collection. */
    private static final long MATCHER_SHUTDOWN_TIMEOUT_MS = 2_000L;

    /**
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static abstract class AbstractGraphMatching {

        public static void setMCSUpdationFlags(Holder holder, int substrateIndex, int productIndex) throws Exception {
            ReactionContainer reactionStructureInformation = holder.getReactionContainer();
            reactionStructureInformation.setEductModified(substrateIndex, true);
            reactionStructureInformation.setProductModified(productIndex, true);
        }

        public abstract IAtomContainer getMatchedPart();

        public abstract IAtomContainer getRemainingEduct();

        public abstract IAtomContainer getRemainingProduct();

        public abstract boolean mcsMatch(Holder holder, boolean removeHydrogen, Integer I, Integer J, BitSet eductFP, BitSet prodFP);

        public abstract int removeMatchedAtomsAndUpdateAAM(IReaction reaction);
    }

    private static final int LARGE_JOB_THRESHOLD = 1000;
    private final static ILoggingTool LOGGER
            = createLoggingTool(GraphMatcher.class);
    private static final ReactionMappingEngine MAPPING_ENGINE
            = SmsdReactionMappingEngine.getInstance();

    static MatcherSettings matcherSettingsFor(IMappingAlgorithm theory,
            int numberOfCyclesEduct, int numberOfCyclesProduct,
            boolean hasPerfectRings) {
        boolean ringFlag = numberOfCyclesEduct > 0 && numberOfCyclesProduct > 0;
        switch (theory) {
            case RINGS:
                return new MatcherSettings(false, false, ringFlag, hasPerfectRings);
            case MAX:
                return new MatcherSettings(false, true, hasPerfectRings, false);
            case MIN:
            case MIXTURE:
            default:
                return new MatcherSettings(false, false, hasPerfectRings, false);
        }
    }

    static final class MatcherSettings implements Serializable {

        private static final long serialVersionUID = 0x2f0f0bbce57fL;
        private final boolean atomType;
        private final boolean bondMatch;
        private final boolean ringMatch;
        private final boolean ringSizeMatch;

        MatcherSettings(boolean atomType, boolean bondMatch,
                boolean ringMatch, boolean ringSizeMatch) {
            this.atomType = atomType;
            this.bondMatch = bondMatch;
            this.ringMatch = ringMatch;
            this.ringSizeMatch = ringSizeMatch;
        }
    }

    private static final class PairJob {

        private final Combination representative;
        private final List<Combination> occurrences;
        private final MatcherSettings settings;
        private final boolean hasPerfectRings;
        private final int numberOfCyclesEduct;
        private final int numberOfCyclesProduct;
        private final String queryStructureKey;
        private final String targetStructureKey;

        PairJob(Combination representative,
                MatcherSettings settings,
                boolean hasPerfectRings,
                int numberOfCyclesEduct,
                int numberOfCyclesProduct,
                String queryStructureKey,
                String targetStructureKey) {
            this.representative = representative;
            this.occurrences = new ArrayList<>();
            this.occurrences.add(representative);
            this.settings = settings;
            this.hasPerfectRings = hasPerfectRings;
            this.numberOfCyclesEduct = numberOfCyclesEduct;
            this.numberOfCyclesProduct = numberOfCyclesProduct;
            this.queryStructureKey = queryStructureKey;
            this.targetStructureKey = targetStructureKey;
        }

        void addOccurrence(Combination occurrence) {
            occurrences.add(occurrence);
        }
    }

    private static void harmonizeForSmsd(IAtomContainer container) {
        if (container == null) {
            return;
        }
        for (IBond bond : container.bonds()) {
            if (bond == null) {
                continue;
            }
            if (bond.getOrder() == null || bond.getOrder() == IBond.Order.UNSET) {
                bond.setOrder(IBond.Order.SINGLE);
            }
            if (bond.isAromatic()) {
                if (bond.getBegin() != null) {
                    bond.getBegin().setIsAromatic(true);
                }
                if (bond.getEnd() != null) {
                    bond.getEnd().setIsAromatic(true);
                }
            }
        }
    }

    /**
     *
     * @param mh
     * @return
     * @throws InterruptedException
     */
    public static Collection<MCSSolution> matcher(Holder mh) throws Exception {
        ExecutorService executor = null;
        Collection<MCSSolution> mcsSolutions = new ArrayList<>();
        long matcherStart = currentTimeMillis();
        String reactionId = mh.getReactionID();
        String algorithmName = mh.getTheory() == null ? "UNKNOWN" : mh.getTheory().name();

        LOGGER.debug("Matcher Class for " + mh.getTheory());
        List<Combination> jobReplicatorList = new ArrayList<>();
        int taskCounter = 0;

        try {
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Integer eductCount = reactionStructureInformation.getEductCount();
            Integer productCount = reactionStructureInformation.getProductCount();
            for (int substrateIndex = 0; substrateIndex < eductCount; substrateIndex++) {
                for (int productIndex = 0; productIndex < productCount; productIndex++) {
                    IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                    IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                    boolean hasAtoms = educt != null && product != null
                            && educt.getAtomCount() > 0
                            && product.getAtomCount() > 0;
                    boolean forceInitial = mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1;
                    boolean needsRefresh = forceInitial
                            || reactionStructureInformation.isEductModified(substrateIndex)
                            || reactionStructureInformation.isProductModified(productIndex);
                    LOGGER.debug("educt atoms " + (educt == null ? 0 : educt.getAtomCount())
                            + ", product atoms " + (product == null ? 0 : product.getAtomCount())
                            + ", needsRefresh " + needsRefresh);
                    if ((hasAtoms || forceInitial) && needsRefresh) {
                        Combination c = new Combination(substrateIndex, productIndex);
                        jobReplicatorList.add(c);
                    }
                }
            }

            LOGGER.debug("jobReplicatorList " + jobReplicatorList.size());

            if (jobReplicatorList.isEmpty()) {
                return unmodifiableCollection(mcsSolutions);
            }

            /*
             * Pre-compute aromaticity and cycle counts ONCE per molecule.
             * Previously this ran for every educt×product pair — O(E*P) redundancy.
             */
            CycleFinder allCycles = Cycles.or(Cycles.all(),
                    Cycles.or(Cycles.relevant(), Cycles.essential()));
            Aromaticity aromaticity = new Aromaticity(daylight(), allCycles);
            CycleFinder shortCycles = Cycles.vertexShort();

            int[] eductCycleCache = new int[eductCount];
            String[] eductStructureKeys = new String[eductCount];
            for (int i = 0; i < eductCount; i++) {
                IAtomContainer educt = reactionStructureInformation.getEduct(i);
                if (educt != null && educt.getAtomCount() > 0) {
                    harmonizeForSmsd(educt);
                    try {
                        aromaticity.apply(educt);
                        eductCycleCache[i] = shortCycles.find(educt).numberOfCycles();
                    } catch (CDKException | RuntimeException ex) {
                        eductCycleCache[i] = 0;
                    }
                    eductStructureKeys[i] = MappingKeyUtil.computeStructureKey(educt);
                }
            }

            int[] productCycleCache = new int[productCount];
            String[] productStructureKeys = new String[productCount];
            for (int j = 0; j < productCount; j++) {
                IAtomContainer product = reactionStructureInformation.getProduct(j);
                if (product != null && product.getAtomCount() > 0) {
                    harmonizeForSmsd(product);
                    try {
                        aromaticity.apply(product);
                        productCycleCache[j] = shortCycles.find(product).numberOfCycles();
                    } catch (CDKException | RuntimeException ex) {
                        productCycleCache[j] = 0;
                    }
                    productStructureKeys[j] = MappingKeyUtil.computeStructureKey(product);
                }
            }

            Map<String, PairJob> pairJobs = new LinkedHashMap<>();
            for (Combination c : jobReplicatorList) {
                int substrateIndex = c.getRowIndex();
                int productIndex = c.getColIndex();
                int numberOfCyclesEduct = eductCycleCache[substrateIndex];
                int numberOfCyclesProduct = productCycleCache[productIndex];
                boolean ringSizeEqual = (numberOfCyclesEduct == numberOfCyclesProduct);
                MatcherSettings settings = matcherSettingsFor(
                        mh.getTheory(),
                        numberOfCyclesEduct,
                        numberOfCyclesProduct,
                        ringSizeEqual);
                String queryStructureKey = eductStructureKeys[substrateIndex] == null ? "" : eductStructureKeys[substrateIndex];
                String targetStructureKey = productStructureKeys[productIndex] == null ? "" : productStructureKeys[productIndex];
                String pairKey = MappingKeyUtil.buildPairKey(
                        queryStructureKey,
                        targetStructureKey,
                        mh.getTheory().name(),
                        settings.atomType,
                        settings.bondMatch,
                        settings.ringMatch,
                        settings.ringSizeMatch);
                PairJob pairJob = pairJobs.get(pairKey);
                if (pairJob == null) {
                    pairJobs.put(pairKey, new PairJob(
                            c,
                            settings,
                            ringSizeEqual,
                            numberOfCyclesEduct,
                            numberOfCyclesProduct,
                            queryStructureKey,
                            targetStructureKey));
                } else {
                    pairJob.addOccurrence(c);
                }
            }

            /*
             * Assign the threads
             *
             * Use Single Thread to computed MCS as muntiple threads lock the calculations!
             */
            int threadsAvailable = Math.max(1, getRuntime().availableProcessors() - 1);
            threadsAvailable = Math.max(1, Math.min(threadsAvailable, pairJobs.size()));

            LOGGER.debug("Candidate pairs " + jobReplicatorList.size()
                    + ", unique structural pairs " + pairJobs.size());
            LOGGER.debug(threadsAvailable + " threads requested for MCS in " + mh.getTheory());

            executor = Executors.newFixedThreadPool(threadsAvailable);
            CompletionService<MCSSolution> callablesQueue = new ExecutorCompletionService<>(executor);

            List<PairJob> jobsToRun = new ArrayList<>();
            List<MCSThread> listOfJobs = new ArrayList<>();
            Map<Combination, PairJob> pairJobsByRepresentative = new HashMap<>();

            int skippedIdentity = 0, skippedRatio = 0, skippedTanimoto = 0;
            List<MCSSolution> directMCSSolutions = new ArrayList<>();

            for (PairJob pairJob : pairJobs.values()) {
                Combination representative = pairJob.representative;
                int substrateIndex = representative.getRowIndex();
                int productIndex = representative.getColIndex();
                IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                pairJobsByRepresentative.put(representative, pairJob);

                /*
                 * PRE-FILTER 1: Identity — if structural keys match, build the
                 * atom mapping directly (atom i → atom i). Do NOT run SMSD: identical
                 * molecules can have multiple valid MCS solutions due to symmetry,
                 * and SMSD may return a non-identity mapping that causes spurious
                 * bond changes in the calculator.
                 */
                if (!pairJob.queryStructureKey.isEmpty()
                        && pairJob.queryStructureKey.equals(pairJob.targetStructureKey)
                        && educt.getAtomCount() == product.getAtomCount()) {
                    try {
                        IAtomContainer eductClone = cloneWithIDs(educt);
                        IAtomContainer productClone = cloneWithIDs(product);
                        AtomAtomMapping identityAAM = new AtomAtomMapping(eductClone, productClone);
                        for (int ai = 0; ai < eductClone.getAtomCount(); ai++) {
                            identityAAM.put(eductClone.getAtom(ai), productClone.getAtom(ai));
                        }
                        MCSSolution identityMCS = new MCSSolution(substrateIndex, productIndex,
                                eductClone, productClone, identityAAM);
                        directMCSSolutions.add(identityMCS);
                        skippedIdentity++;
                        continue;
                    } catch (Exception ex) {
                        LOGGER.debug("Identity shortcut failed, falling back to MCS: " + ex.getMessage());
                    }
                }

                /*
                 * PRE-FILTER 2: Atom count ratio — skip pairs where the smaller
                 * molecule is < 30% of the larger. Such pairs rarely contribute
                 * meaningful mappings and waste MCS computation.
                 */
                int eAtoms = educt.getAtomCount();
                int pAtoms = product.getAtomCount();
                if (eAtoms > 0 && pAtoms > 0) {
                    double ratio = (double) Math.min(eAtoms, pAtoms) / Math.max(eAtoms, pAtoms);
                    if (ratio < 0.3 && Math.min(eAtoms, pAtoms) > 3) {
                        skippedRatio++;
                        continue;
                    }
                }

                /*
                 * PRE-FILTER 3: Tanimoto similarity — skip pairs with very low
                 * fingerprint similarity. These molecules share almost no structure.
                 */
                double tanimoto = mh.getFPSimilarityMatrix().getValue(substrateIndex, productIndex);
                if (tanimoto >= 0 && tanimoto < 0.05 && eAtoms > 5 && pAtoms > 5) {
                    skippedTanimoto++;
                    continue;
                }

                jobsToRun.add(pairJob);
            }

            if (skippedIdentity + skippedRatio + skippedTanimoto > 0) {
                LOGGER.debug("Pre-filter: skipped " + skippedIdentity + " identity, "
                        + skippedRatio + " ratio, " + skippedTanimoto + " tanimoto pairs");
            }

            int invocationIndex = MappingDiagnostics.recordMatcherInvocation(
                    reactionId,
                    algorithmName,
                    jobReplicatorList.size(),
                    pairJobs.size(),
                    skippedIdentity,
                    skippedRatio,
                    skippedTanimoto,
                    jobsToRun.size());

            for (PairJob pairJob : jobsToRun) {
                Combination representative = pairJob.representative;
                int substrateIndex = representative.getRowIndex();
                int productIndex = representative.getColIndex();
                IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                MCSThread mcsThread = new MCSThread(mh.getTheory(),
                        substrateIndex, productIndex, educt, product,
                        reactionId, algorithmName, invocationIndex);
                mcsThread.setHasPerfectRings(pairJob.hasPerfectRings);
                mcsThread.setEductRingCount(pairJob.numberOfCyclesEduct);
                mcsThread.setProductRingCount(pairJob.numberOfCyclesProduct);
                listOfJobs.add(mcsThread);
            }

            if (listOfJobs.size() > LARGE_JOB_THRESHOLD) {
                LOGGER.warn("Large job: " + listOfJobs.size() + " MCS pairs to compute");
            }
            if (!listOfJobs.isEmpty()) {
                for (MCSThread mcsThreadJob : listOfJobs) {
                    callablesQueue.submit(mcsThreadJob);
                    taskCounter++;
                }
            }

            LOGGER.debug("submited " + taskCounter + " jobs");
            Collection<MCSSolution> threadedUniqueMCSSolutions = new ArrayList<>();
            int collected = 0;
            long matcherDeadline = currentTimeMillis() + MATCHER_BUDGET_MS;
            for (int count = 0; count < taskCounter; count++) {
                try {
                    long remaining = matcherDeadline - currentTimeMillis();
                    if (remaining <= 0) {
                        LOGGER.warn("Matcher budget (" + MATCHER_BUDGET_MS
                                + "ms) exhausted — " + (taskCounter - collected)
                                + " remaining pairs skipped");
                        break;
                    }
                    long pollMs = Math.min(remaining, MCS_POLL_TIMEOUT_MS);
                    java.util.concurrent.Future<MCSSolution> future =
                            callablesQueue.poll(pollMs, TimeUnit.MILLISECONDS);
                    if (future == null) {
                        LOGGER.warn("MCS poll timed out after " + pollMs
                                + "ms — " + (taskCounter - collected)
                                + " remaining pairs will be skipped");
                        break;
                    }
                    MCSSolution isomorphism = future.get(); // already complete
                    if (isomorphism != null) {
                        threadedUniqueMCSSolutions.add(isomorphism);
                    }
                    collected++;
                } catch (ExecutionException ex) {
                    collected++;
                    Throwable cause = ex.getCause() != null ? ex.getCause() : ex;
                    LOGGER.error(SEVERE, "MCS worker failed", cause);
                }
            }
            // Add directly-constructed identity mappings (bypassed MCSThread)
            threadedUniqueMCSSolutions.addAll(directMCSSolutions);
            // Shut down the local executor; interrupt any stuck threads
            executor.shutdown();
            if (!executor.awaitTermination(MATCHER_SHUTDOWN_TIMEOUT_MS, TimeUnit.MILLISECONDS)) {
                LOGGER.warn("MCS executor did not shut down cleanly — forcing shutdown");
                executor.shutdownNow();
            }

            LOGGER.debug("==Gathering MCS solution from the Thread==");
            long replayedMappings = 0;
            threadedUniqueMCSSolutions.stream().filter((mcs) -> !(mcs == null)).forEach((MCSSolution mcs) -> {
                Combination representative = new Combination(
                        mcs.getQueryPosition(),
                        mcs.getTargetPosition());
                PairJob pairJob = pairJobsByRepresentative.get(representative);
                if (pairJob == null) {
                    return;
                }
                for (Combination occurrence : pairJob.occurrences) {
                    MCSSolution replicatedMCS = replicateMappingOnContainers(mh, occurrence, mcs);
                    if (replicatedMCS == null) {
                        continue;
                    }
                    LOGGER.debug("======MCSSolution======");
                    LOGGER.debug("MCS " + " I " + occurrence.getRowIndex()
                            + " J " + occurrence.getColIndex()
                            + " Number of Atom Mapped " + mcs.getAtomAtomMapping().getCount()
                            + " Number of Atom Mapped replicatedMCS " + replicatedMCS.getAtomAtomMapping().getCount());
                    mcsSolutions.add(replicatedMCS);
                }
            });
            replayedMappings = mcsSolutions.size();
            MappingDiagnostics.recordMatcherCompletion(
                    reactionId,
                    algorithmName,
                    invocationIndex,
                    replayedMappings,
                    currentTimeMillis() - matcherStart);
            jobReplicatorList.clear();

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        } finally {
            if (executor != null) {
                executor.shutdownNow();
            }
        }
        return unmodifiableCollection(mcsSolutions);
    }

    /**
     *
     * @param mh
     * @param solution
     * @param mcs
     * @return
     */
    static MCSSolution replicateMappingOnContainers(Holder mh, Combination solution, MCSSolution mcs) {
        try {
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            IAtomContainer q = reactionStructureInformation.getEduct(solution.getRowIndex());
            IAtomContainer t = reactionStructureInformation.getProduct(solution.getColIndex());

            int diff1 = q.getAtomCount() - mcs.getQueryContainer().getAtomCount();
            int diff2 = t.getAtomCount() - mcs.getTargetContainer().getAtomCount();

            if (diff1 != 0 && diff2 != 0) {
                LOGGER.debug(NEW_LINE + NEW_LINE + " " + solution.getRowIndex() + ", Diff in ac1 " + diff1);
                LOGGER.debug(solution.getColIndex() + ", Diff in ac2 " + diff2);
                LOGGER.debug(NEW_LINE + "ac1 " + q.getAtomCount());
                LOGGER.debug(NEW_LINE + "ac2 " + t.getAtomCount());
                LOGGER.debug(NEW_LINE + "mac1 " + mcs.getQueryContainer().getAtomCount());
                LOGGER.debug(NEW_LINE + "mac2 " + mcs.getTargetContainer().getAtomCount());
            }

            AtomAtomMapping atomAtomMapping = mcs.getAtomAtomMapping();
            AtomAtomMapping atomAtomMappingNew = new AtomAtomMapping(q, t);

            atomAtomMapping.getMappingsByIndex().forEach((queryIndex, targetIndex) -> {
                if (queryIndex >= 0 && queryIndex < q.getAtomCount()
                        && targetIndex >= 0 && targetIndex < t.getAtomCount()) {
                    atomAtomMappingNew.put(q.getAtom(queryIndex), t.getAtom(targetIndex));
                } else {
                    LOGGER.error(WARNING, "Unexpected atom index while replaying cached mapping");
                }
            });

            LOGGER.debug("------Mapped PAIRS------");
            LOGGER.debug("Query " + q.getAtomCount());
            LOGGER.debug("Target " + t.getAtomCount());
            LOGGER.debug("Mapping Size " + atomAtomMappingNew.getCount());
            return new MCSSolution(solution.getRowIndex(), solution.getColIndex(), q, t, atomAtomMappingNew);
        } catch (IOException | CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return null;
    }



    /**
     *
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class Combination implements Serializable, Comparable<Combination>, Comparator<Combination> {

        private static final long serialVersionUID = 786786786L;
        private final int row;
        private final int col;

        /**
         *
         * @param row
         * @param col
         */
        public Combination(int row, int col) {
            this.row = row;
            this.col = col;
        }

        /**
         *
         * @return
         */
        public int getRowIndex() {
            return row;
        }

        /**
         *
         * @return
         */
        public int getColIndex() {
            return col;
        }

        @Override
        public String toString() {
            return "Combination{" + "row=" + row + ", col=" + col + '}';
        }

        @Override
        public int hashCode() {
            int hash = 7;
            hash = 13 * hash + this.row;
            hash = 13 * hash + this.col;
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
            final Combination other = (Combination) obj;
            if (this.row != other.row) {
                return false;
            }
            return this.col == other.col;
        }

        @Override
        public int compareTo(Combination o) {
            int rowCompare = Integer.compare(this.row, o.row);
            return rowCompare != 0 ? rowCompare : Integer.compare(this.col, o.col);
        }

        @Override
        public int compare(Combination o1, Combination o2) {
            return o1.compareTo(o2);
        }
    }



    /**
     *
     * @author Syed Asad Rahman, BioInception
     * @contact asad.rahman@bioinceptionlabs.com
     */
    public static class GraphMatching extends AbstractGraphMatching implements Serializable {

        private final static ILoggingTool LOGGER = createLoggingTool(GraphMatching.class);
        private static final long serialVersionUID = 0xf06b2d5f9L;
        private final IAtomContainer educt;
        private final IAtomContainer product;
        private IAtomContainer matchedPart = null;
        private Map<IAtom, IAtom> bestAtomMappingList;
        private int fragmentCount = 0;

        /**
         * Creates a new instance of GraphMatching
         *
         * @param reaction_ID
         * @param eductOrg
         * @param productOrg
         * @param suffix
         * @param removeHydrogen
         * @throws Exception
         */
        public GraphMatching(String reaction_ID, IAtomContainer eductOrg, IAtomContainer productOrg, String suffix, boolean removeHydrogen) throws Exception {

            try {

                educt = eductOrg;
                product = productOrg;
                educt.setID(eductOrg.getID());
                product.setID(productOrg.getID());

                if (educt.getAtomCount() > 0 && product.getAtomCount() > 0) {
                    setMatchedPart(cloneWithIDs(educt));
                }
            } catch (CloneNotSupportedException e) {
                throw new CDKException("Error: In GraphMatching Class" + e);
            }

        }

        /**
         *
         * @param holder
         * @param removeHydrogen
         * @param substrateIndex
         * @param productIndex
         * @param eductFP
         * @param prodFP
         * @return
         */
        @Override
        public boolean mcsMatch(Holder holder,
                boolean removeHydrogen,
                Integer substrateIndex,
                Integer productIndex,
                BitSet eductFP,
                BitSet prodFP) {

            if (educt.getAtomCount() <= 0 && product.getAtomCount() <= 0) {
                return false;
            }

            try {
                try {
                    setMCSUpdationFlags(holder, substrateIndex, productIndex);
                } catch (Exception ex) {
                    LOGGER.error(Level.SEVERE, null, ex);
                }
                BestMatch initMCSAtom = holder.getBestMatchContainer();
                if (initMCSAtom.containsKey(substrateIndex, productIndex)) {
                    this.bestAtomMappingList = initMCSAtom.getAtomMatch(substrateIndex, productIndex).getMappingsByAtoms();
                    this.fragmentCount = initMCSAtom.getTotalFragmentCount(substrateIndex, productIndex);
                    if (this.bestAtomMappingList != null && !this.bestAtomMappingList.isEmpty()) {
                        return true;
                    }
                }
            } catch (IOException ex) {
                LOGGER.debug("Files: " + educt.getID() + ", " + product.getID());
                LOGGER.debug(SEVERE, null, ex);
            }
            return false;
        }

        /**
         *
         * @param reaction
         * @return
         */
        @Override
        public int removeMatchedAtomsAndUpdateAAM(IReaction reaction) {
            int delta = 0;

            LOGGER.debug("Before removing Mol Size E: " + educt.getAtomCount()
                    + " , Before removing Mol Size P: " + product.getAtomCount());
            int beforeESize = educt.getAtomCount();
            Map<String, IAtom> eductAtomsById = indexAtomsById(educt);
            Map<String, IAtom> productAtomsById = indexAtomsById(product);
            Map<String, IAtom> matchedAtomsById = indexAtomsById(matchedPart);

            if (bestAtomMappingList != null) {
                for (Map.Entry<IAtom, IAtom> map : bestAtomMappingList.entrySet()) {
                    String eID = map.getKey().getID();
                    IAtom eAtom = eductAtomsById.get(eID);
                    String pID = map.getValue().getID();
                    LOGGER.debug("eID " + eID + ",pID " + pID);
                    IAtom pAtom = productAtomsById.get(pID);

                    if (eAtom != null && pAtom != null) {
                        IMapping im = SilentChemObjectBuilder.getInstance().newInstance(IMapping.class, eAtom, pAtom);
                        reaction.addMapping(im);
                    }
                    if (eAtom != null) {
                        educt.removeAtom(eAtom);
                        eductAtomsById.remove(eID);
                    }
                    if (pAtom != null) {
                        product.removeAtom(pAtom);
                        productAtomsById.remove(pID);
                    }
                    delta = fragmentCount;
                }
            }

            for (IAtom atom : educt.atoms()) {
                IAtom matchedAtom = matchedAtomsById.get(atom.getID());
                if (matchedAtom != null) {
                    matchedPart.removeAtom(matchedAtom);
                    matchedAtomsById.remove(atom.getID());
                }
            }

            LOGGER.debug("After removing Mol Size E: " + educt.getAtomCount()
                    + " , After removing Mol Size P: " + product.getAtomCount());

            if (beforeESize == educt.getAtomCount()) {
                try {
                    throw new CDKException("Failed to remove matched parts between " + educt.getID() + ": "
                            + educt.getAtomCount() + " , " + product.getID() + " : " + product.getAtomCount()
                            + ", Mapping count: " + bestAtomMappingList.size() + "...atom ids did not matched!");
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, "Failed to remove matched parts between " + educt.getID() + ": "
                            + educt.getAtomCount() + " , " + product.getID() + " : " + product.getAtomCount()
                            + ", Mapping count: " + bestAtomMappingList.size() + "...atom ids did not matched!", ex);

                    throw new RuntimeException("Failed to remove matched parts between " + educt.getID() + ": "
                            + educt.getAtomCount() + " , " + product.getID() + " : " + product.getAtomCount()
                            + ", Mapping count: " + bestAtomMappingList.size() + "...atom ids did not matched!", ex);
                }
            }
            return delta;
        }

        private Map<String, IAtom> indexAtomsById(IAtomContainer container) {
            Map<String, IAtom> atomsById = new HashMap<>();
            if (container == null) {
                return atomsById;
            }
            for (IAtom atom : container.atoms()) {
                if (atom.getID() != null) {
                    atomsById.put(atom.getID(), atom);
                }
            }
            return atomsById;
        }

        /**
         *
         * @return
         */
        @Override
        public IAtomContainer getRemainingEduct() {
            return educt;
        }

        /**
         *
         * @return
         */
        @Override
        public IAtomContainer getRemainingProduct() {
            return product;
        }

        /**
         *
         * @return
         */
        protected Map<IAtom, IAtom> getFirstAtomMapping() {
            return unmodifiableMap(bestAtomMappingList);
        }

        /**
         * @return the matchedPart
         */
        @Override
        public IAtomContainer getMatchedPart() {
            return matchedPart;
        }

        /**
         * @param aMatchedPart the matchedPart to set
         */
        private void setMatchedPart(IAtomContainer aMatchedPart) {
            matchedPart = aMatchedPart;
        }
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class MCSSolution implements Serializable {

        private static final long serialVersionUID = 0xc678991ddf0L;
        private final IAtomContainer queryContainer;
        private final int targetPosition;
        private final IAtomContainer targetContainer;
        private final AtomAtomMapping atomatomMapping;
        private final int queryPosition;
        private Integer stereoScore;
        private Integer fragmentSize;
        private Double energy;

        /**
         *
         * @param queryPosition
         * @param targetPosition
         * @param queryContainer
         * @param targetContainer
         * @param aam
         */
        public MCSSolution(int queryPosition, int targetPosition,
                IAtomContainer queryContainer, IAtomContainer targetContainer, AtomAtomMapping aam) {
            this.queryContainer = queryContainer;
            this.targetPosition = targetPosition;
            this.targetContainer = targetContainer;
            this.atomatomMapping = aam;
            this.queryPosition = queryPosition;
            this.energy = null;
            this.fragmentSize = null;
            this.stereoScore = null;
        }

        /**
         * @return the stereoScore, or 0 if null
         */
        public Integer getStereoScore() {
            return stereoScore != null ? stereoScore : 0;
        }

        /**
         * @param stereoScore the stereoScore to set
         */
        public void setStereoScore(Integer stereoScore) {
            this.stereoScore = stereoScore;
        }

        /**
         * @return the fragmentSize, or 0 if null
         */
        public Integer getFragmentSize() {
            return fragmentSize != null ? fragmentSize : 0;
        }

        /**
         * @param fragmentSize the fragmentSize to set
         */
        public void setFragmentSize(Integer fragmentSize) {
            this.fragmentSize = fragmentSize;
        }

        /**
         * @return the energy, or 0.0 if null
         */
        public Double getEnergy() {
            return energy != null ? energy : 0.0;
        }

        /**
         * @param energy the energy to set
         */
        public void setEnergy(Double energy) {
            this.energy = energy;
        }

        /**
         * @return the queryContainer
         */
        public IAtomContainer getQueryContainer() {
            return queryContainer;
        }

        /**
         * @return the targetContainer
         */
        public IAtomContainer getTargetContainer() {
            return targetContainer;
        }

        /**
         * @return the atomatomMapping
         */
        public AtomAtomMapping getAtomAtomMapping() {
            return atomatomMapping;
        }

        /**
         * @return the targetPosition
         */
        public int getTargetPosition() {
            return targetPosition;
        }

        /**
         * @return the queryPosition
         */
        public int getQueryPosition() {
            return queryPosition;
        }
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class MCSThread implements Callable<MCSSolution> {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(MCSThread.class);

        static final String NEW_LINE = getProperty("line.separator");

        /**
         *
         */
        protected final IAtomContainer compound1;

        /**
         *
         */
        protected final IAtomContainer compound2;

        /**
         *
         */
        protected final int queryPosition;

        /**
         *
         */
        protected final int targetPosition;

        /**
         *
         */
        protected final IMappingAlgorithm theory;
        private final String reactionId;
        private final String algorithmName;
        private final int invocationIndex;
        private final Map<String, Integer> compound1SymbolCounts;
        private final Map<String, Integer> compound2SymbolCounts;
        private final boolean moleculesConnected;

        /**
         *
         */
        long startTime;
        private boolean hasRings;
        private int numberOfCyclesEduct;
        private int numberOfCyclesProduct;

        /**
         *
         * @param theory
         * @param queryPosition
         * @param targetPosition
         * @param educt
         * @param product
         * @param bondMatcher
         * @param ringMatcher
         * @param atomMatcher
         * @throws CloneNotSupportedException
         * @throws org.openscience.cdk.exception.CDKException
         */
        MCSThread(IMappingAlgorithm theory, int queryPosition, int targetPosition,
                IAtomContainer educt, IAtomContainer product,
                String reactionId, String algorithmName, int invocationIndex)
                throws CloneNotSupportedException, CDKException {
            this.compound1 = getNewContainerWithIDs(educt);
            this.compound2 = getNewContainerWithIDs(product);
            this.queryPosition = queryPosition;
            this.targetPosition = targetPosition;
            this.theory = theory;
            this.reactionId = reactionId;
            this.algorithmName = algorithmName;
            this.invocationIndex = invocationIndex;
            this.numberOfCyclesEduct = 0;
            this.numberOfCyclesProduct = 0;
            this.compound1SymbolCounts = countAtomsBySymbol(this.compound1);
            this.compound2SymbolCounts = countAtomsBySymbol(this.compound2);
            this.moleculesConnected = isConnected(this.compound1) && isConnected(this.compound2);
        }

        void printMatch(BaseMapping isomorphism) {
            int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                    : isomorphism.getFirstAtomMapping().getCount();

            LOGGER.debug("Q: " + isomorphism.getQuery().getID()
                    + " T: " + isomorphism.getTarget().getID()
                    + " atoms: " + isomorphism.getQuery().getAtomCount()
                    + " atoms: " + isomorphism.getTarget().getAtomCount()
                    + " overlaps: " + overlap);
        }

        @Override
        public MCSSolution call() throws Exception {
            boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;
            int commonAtomUpperBound = commonAtomUpperBound(compound1SymbolCounts, compound2SymbolCounts);

            AtomMatcher am;
            BondMatcher bm;

            LOGGER.debug("in mcsthread call ");

            try {
                if (commonAtomUpperBound == 0) {
                    return emptySolution();
                }

                MCSSolution singleAtomSolution = singleAtomSolution();
                if (singleAtomSolution != null) {
                    return singleAtomSolution;
                }

                /*
                 * IMP: Do not perform substructure matching for disconnected molecules
                 */
                boolean moleculeConnected = moleculesConnected;
                /*
                     Check if MCS matching required or not very IMP step
                 */
                boolean possibleVFmatch12 = isPossibleSubgraphMatch(compound1SymbolCounts, compound2SymbolCounts);
                LOGGER.debug("VF Matcher 1->2 " + possibleVFmatch12);

                boolean possibleVFmatch21 = isPossibleSubgraphMatch(compound2SymbolCounts, compound1SymbolCounts);
                LOGGER.debug("VF Matcher 2->1 " + possibleVFmatch21);

                if (moleculeConnected && possibleVFmatch12) {
                    LOGGER.debug("Substructure 1");
                    this.startTime = currentTimeMillis();

                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());

                    LOGGER.debug("---1.1---");
                    BaseMapping substructure;
                    am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                    substructure = MAPPING_ENGINE.findSubstructure(ac1, ac2, am, bm, true,
                            SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);

                    if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                        am = AtomBondMatcher.atomMatcher(false, ringFlag);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---1.3---");
                        MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                        substructure = MAPPING_ENGINE.findSubstructure(ac1, ac2,
                                am, bm, true, SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);
                    } else if (moleculeConnected && !substructure.isSubgraph()) {
                        am = AtomBondMatcher.atomMatcher(false, false);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---1.2---");
                        MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                        substructure = MAPPING_ENGINE.findSubstructure(ac1, ac2, am, bm, true,
                                SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);
                    }
                    MAPPING_ENGINE.applyDefaultFilters(substructure);
                    if (substructure.isSubgraph()
                            && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
                        LOGGER.debug("Found Substructure 1");
                        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                                substructure.getQuery(), substructure.getTarget(), substructure.getFirstAtomMapping());
                        mcs.setEnergy(substructure.getEnergyScore(0));
                        mcs.setFragmentSize(substructure.getFragmentSize(0));
                        mcs.setStereoScore(substructure.getStereoScore(0));
                        long stopTime = currentTimeMillis();
                        long time = stopTime - startTime;
                        printMatch(substructure);
                        LOGGER.debug("\" Time:\" " + time);
                        return mcs;
                    } else {
                        LOGGER.debug("not a Substructure 1");
                    }
                }

                if (moleculeConnected && !possibleVFmatch12 && possibleVFmatch21) {

                    LOGGER.debug("Substructure 2");
                    this.startTime = currentTimeMillis();

                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());
                    BaseMapping substructure;

                    LOGGER.debug("---2.1---");
                    am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                    substructure = MAPPING_ENGINE.findSubstructure(ac2, ac1, am, bm, true,
                            SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);

                    if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                        am = AtomBondMatcher.atomMatcher(false, ringFlag);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---2.3---");
                        MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                        substructure = MAPPING_ENGINE.findSubstructure(ac2, ac1, am, bm, true,
                                SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);
                    } else if (moleculeConnected && !substructure.isSubgraph()) {
                        am = AtomBondMatcher.atomMatcher(false, false);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---2.2---");
                        MappingDiagnostics.recordSubstructureSearch(reactionId, algorithmName, invocationIndex);
                        substructure = MAPPING_ENGINE.findSubstructure(ac2, ac1, am, bm, true,
                                SINGLE_SUBGRAPH_MATCH, SUBGRAPH_TIMEOUT_MS);
                    }
                    MAPPING_ENGINE.applyDefaultFilters(substructure);

                    if (substructure.isSubgraph()
                            && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {

                        LOGGER.debug("Found Substructure 2");
                        AtomAtomMapping aam = new AtomAtomMapping(substructure.getTarget(), substructure.getQuery());
                        Map<IAtom, IAtom> mappings = substructure.getFirstAtomMapping().getMappingsByAtoms();
                        mappings.keySet().stream().forEach((atom1) -> {
                            IAtom atom2 = mappings.get(atom1);
                            aam.put(atom2, atom1);
                        });
                        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                                substructure.getTarget(), substructure.getQuery(), aam);
                        mcs.setEnergy(substructure.getEnergyScore(0));
                        mcs.setFragmentSize(substructure.getFragmentSize(0));
                        mcs.setStereoScore(substructure.getStereoScore(0));

                        long stopTime = currentTimeMillis();
                        long time = stopTime - startTime;
                        printMatch(substructure);
                        LOGGER.debug("\" Time:\" " + time);
                        return mcs;
                    } else {
                        LOGGER.debug("not a Substructure 2");
                    }
                }

                /*
                 * If substructure matches have failed then call MCS
                 */
                LOGGER.debug("==============================================");
                LOGGER.debug("No Substructure found - switching to MCS");
                LOGGER.debug("Q: " + getCompound1().getID()
                        + NEW_LINE
                        + " T: " + getCompound2().getID()
                        + NEW_LINE
                        + " atomsE: " + compound1.getAtomCount()
                        + " atomsP: " + compound2.getAtomCount());
                LOGGER.debug("==============================================");
                this.startTime = currentTimeMillis();
                MCSSolution mcs = mcs();
                long stopTime = currentTimeMillis();
                long time = stopTime - startTime;
                LOGGER.debug("\"MCS Time:\" " + time);
                return mcs;

            } catch (CDKException | CloneNotSupportedException | RuntimeException ex) {
                LOGGER.error(SEVERE, "Error in generating MCS Solution", ex);
            }
            return null;
        }

        /**
         * Clone molecule preserving IDs.
         * Aromaticity and atom-type perception already done in GraphMatcher
         * before MCSThread is created — do NOT repeat here (was ~25% of total time).
         */
        private IAtomContainer getNewContainerWithIDs(IAtomContainer mol)
                throws CDKException, CloneNotSupportedException {
            if (mol != null && mol.getAtomCount() > 0) {
                IAtomContainer ac = ExtAtomContainerManipulator.cloneWithIDs(mol);

                for (int i = 0; i < ac.getAtomCount(); i++) {
                    String atomID = mol.getAtom(i).getID() == null
                            ? valueOf(i) : mol.getAtom(i).getID();
                    ac.getAtom(i).setID(atomID);
                }
                String containerID = mol.getID() == null ? valueOf(nanoTime()) : mol.getID();
                ac.setID(containerID);
                harmonizeForSmsd(ac);

                return ac;
            }
            return mol;
        }

        private boolean isPossibleSubgraphMatch(Map<String, Integer> queryAtomCounts,
                Map<String, Integer> targetAtomCounts) {
            if (queryAtomCounts.size() > targetAtomCounts.size()) {
                return false;
            }

            // Check all atom types in query exist in target with sufficient count
            for (Map.Entry<String, Integer> entry : queryAtomCounts.entrySet()) {
                Integer targetCount = targetAtomCounts.get(entry.getKey());
                if (targetCount == null || entry.getValue() > targetCount) {
                    return false;
                }
            }

            return true;
        }

        private int commonAtomUpperBound(Map<String, Integer> leftAtomCounts,
                Map<String, Integer> rightAtomCounts) {
            int common = 0;
            for (Map.Entry<String, Integer> entry : leftAtomCounts.entrySet()) {
                Integer rightCount = rightAtomCounts.get(entry.getKey());
                if (rightCount != null) {
                    common += Math.min(entry.getValue(), rightCount);
                }
            }
            return common;
        }

        MCSSolution mcs() throws CDKException, CloneNotSupportedException {

            LOGGER.debug("=============MCS============");
            /*
             * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
             */
            IAtomContainer ac1 = duplicate(getCompound1());
            IAtomContainer ac2 = duplicate(getCompound2());

            // Guard: cannot compute MCS on empty molecules
            if (ac1 == null || ac2 == null || ac1.getAtomCount() == 0 || ac2.getAtomCount() == 0) {
                return null;
            }
                BaseMapping isomorphism;
            MatcherSettings settings = matcherSettingsFor(
                    theory,
                    numberOfCyclesEduct,
                    numberOfCyclesProduct,
                    isHasPerfectRings());

            String key;
            MCSSolution mcs;
            AtomMatcher am;
            BondMatcher bm;
            am = AtomBondMatcher.atomMatcher(settings.atomType, settings.ringSizeMatch);
            bm = AtomBondMatcher.bondMatcher(settings.bondMatch, settings.ringMatch);

            key = generateUniqueKey(settings);
            ThreadSafeCache<String, MCSSolution> mappingCache = ThreadSafeCache.getInstance();
            MCSSolution cachedSolution = mappingCache.get(key);
            if (cachedSolution != null) {
                LOGGER.debug("===={Aladdin} Mapping {Gini}====");
                MappingDiagnostics.recordMcsCacheHit(reactionId, algorithmName, invocationIndex);
                mcs = copyOldSolutionToNew(
                        getQueryPosition(), getTargetPosition(),
                        getCompound1(), getCompound2(),
                        cachedSolution);

            } else {
                SearchEngine.McsOptions mcsOptions = new SearchEngine.McsOptions();
                mcsOptions.timeoutMs = MCS_TIMEOUT_MS;
                mcsOptions.connectedOnly = moleculesConnected;
                mcsOptions.disconnectedMCS = !mcsOptions.connectedOnly;
                mcsOptions.maximizeBonds = settings.bondMatch;
                MappingDiagnostics.recordActualMcsSearch(reactionId, algorithmName, invocationIndex);
                isomorphism = MAPPING_ENGINE.findMcs(ac1, ac2, Algorithm.VFLibMCS, am, bm, mcsOptions);
                mcs = addMCSSolution(key, mappingCache, isomorphism);
            }

            return mcs;

        }

        private Map<String, Integer> countAtomsBySymbol(IAtomContainer container) {
            Map<String, Integer> counts = new HashMap<>();
            if (container == null) {
                return counts;
            }
            for (IAtom atom : container.atoms()) {
                counts.merge(atom.getSymbol(), 1, Integer::sum);
            }
            return counts;
        }

        private IAtomContainer duplicate(IAtomContainer ac) throws CloneNotSupportedException {
            IAtomContainer a = ac.clone();
            a.setID(ac.getID());
            a.setProperties(ac.getProperties());

            for (int i = 0; i < a.getAtomCount(); i++) {
                a.getAtom(i).setID(ac.getAtom(i).getID());
            }

            harmonizeForSmsd(a);

            return a;
        }

        private MCSSolution emptySolution() {
            return new MCSSolution(
                    getQueryPosition(),
                    getTargetPosition(),
                    getCompound1(),
                    getCompound2(),
                    new AtomAtomMapping(getCompound1(), getCompound2()));
        }

        private MCSSolution singleAtomSolution() {
            IAtomContainer query = getCompound1();
            IAtomContainer target = getCompound2();
            if (query == null || target == null) {
                return null;
            }
            if (Math.min(query.getAtomCount(), target.getAtomCount()) != 1) {
                return null;
            }
            IAtom queryAtom = query.getAtomCount() == 1 ? query.getAtom(0) : null;
            IAtom targetAtom = target.getAtomCount() == 1 ? target.getAtom(0) : null;

            if (queryAtom != null) {
                for (IAtom candidate : target.atoms()) {
                    if (queryAtom.getSymbol().equalsIgnoreCase(candidate.getSymbol())) {
                        return singleAtomMapping(query, target, queryAtom, candidate);
                    }
                }
                return emptySolution();
            }

            if (targetAtom != null) {
                for (IAtom candidate : query.atoms()) {
                    if (candidate.getSymbol().equalsIgnoreCase(targetAtom.getSymbol())) {
                        return singleAtomMapping(query, target, candidate, targetAtom);
                    }
                }
                return emptySolution();
            }
            return null;
        }

        private MCSSolution singleAtomMapping(IAtomContainer query, IAtomContainer target,
                IAtom queryAtom, IAtom targetAtom) {
            AtomAtomMapping mapping = new AtomAtomMapping(query, target);
            mapping.put(queryAtom, targetAtom);
            MCSSolution solution = new MCSSolution(
                    getQueryPosition(), getTargetPosition(), query, target, mapping);
            solution.setFragmentSize(1);
            solution.setStereoScore(0);
            solution.setEnergy(0.0);
            return solution;
        }

        /**
         * @return the compound1
         */
        IAtomContainer getCompound1() {
            return compound1;
        }

        /**
         * @return the compound2
         */
        IAtomContainer getCompound2() {
            return compound2;
        }

        /**
         * @return the queryPosition
         */
        int getQueryPosition() {
            return queryPosition;
        }

        /**
         * @return the targetPosition
         */
        int getTargetPosition() {
            return targetPosition;
        }

        void setHasPerfectRings(boolean ring) {
            this.hasRings = ring;
        }

        /**
         * @return the hasRings
         */
        boolean isHasPerfectRings() {
            return hasRings;
        }

        /*
         * Check if fragmented container has single atom
         */
        private boolean isConnected(IAtomContainer compound) {
            LOGGER.debug("isConnected");
            IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound);
            return partitionIntoMolecules.getAtomContainerCount() <= 1;
        }

        void setEductRingCount(int numberOfCyclesEduct) {
            this.numberOfCyclesEduct = numberOfCyclesEduct;
        }

        void setProductRingCount(int numberOfCyclesProduct) {
            this.numberOfCyclesProduct = numberOfCyclesProduct;
        }

        String generateUniqueKey(MatcherSettings settings) {
            return MappingKeyUtil.buildPairKey(
                    compound1,
                    compound2,
                    theory.name(),
                    settings.atomType,
                    settings.bondMatch,
                    settings.ringMatch,
                    settings.ringSizeMatch);
        }

        /*
         * copy old mapping from the cache to new
         */
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

        MCSSolution addMCSSolution(String key, ThreadSafeCache<String, MCSSolution> mappingcache, BaseMapping isomorphism) {

            MAPPING_ENGINE.applyDefaultFilters(isomorphism);
            LOGGER.debug("MCS " + isomorphism.getFirstAtomMapping().getCount());
            /*
             * In case of Complete subgraph, don't use Energy filter
             *
             */
            MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                    isomorphism.getQuery(), isomorphism.getTarget(), isomorphism.getFirstAtomMapping());
            mcs.setEnergy(isomorphism.getEnergyScore(0));
            mcs.setFragmentSize(isomorphism.getFragmentSize(0));
            mcs.setStereoScore(isomorphism.getStereoScore(0));
            long stopTime = currentTimeMillis();
            long time = stopTime - startTime;
            printMatch(isomorphism);
            LOGGER.debug("\" Time:\" " + time);
            MCSSolution cached = mappingcache.putIfAbsent(key, mcs);
            if (cached == mcs) {
                LOGGER.debug("Key " + key);
                LOGGER.debug("mcs size " + mcs.getAtomAtomMapping().getCount());
                LOGGER.debug("mcs map " + mcs.getAtomAtomMapping().getMappingsByIndex());
                LOGGER.debug("\n\n\n ");
                return mcs;
            }
            return copyOldSolutionToNew(
                    getQueryPosition(), getTargetPosition(),
                    getCompound1(), getCompound2(),
                    cached);
        }
    }


}
