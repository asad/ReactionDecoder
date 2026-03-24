/*
 * GraphMatcher - consolidated graph matching classes.
 * Merged: Combination, GraphMatching, MCSSolution, MCSThread into GraphMatcher
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinceptionlabs.reactionblast.mapping.AbstractGraphMatching;
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
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
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.BaseMapping.Algorithm;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.ExtAtomContainerManipulator;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
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

    private static final int LARGE_JOB_THRESHOLD = 1000;
    private final static ILoggingTool LOGGER
            = createLoggingTool(GraphMatcher.class);

    /**
     *
     * @param mh
     * @return
     * @throws InterruptedException
     */
    public static Collection<MCSSolution> matcher(Holder mh) throws Exception {
        ExecutorService executor = null;
        Collection<MCSSolution> mcsSolutions = new ArrayList<>();

        LOGGER.debug("Matcher Class for " + mh.getTheory());
        Set<Combination> jobReplicatorList = new TreeSet<>();
        int taskCounter = 0;

        try {
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Integer eductCount = reactionStructureInformation.getEductCount();
            Integer productCount = reactionStructureInformation.getProductCount();
            for (int substrateIndex = 0; substrateIndex < eductCount; substrateIndex++) {
                for (int productIndex = 0; productIndex < productCount; productIndex++) {
                    IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                    IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                    LOGGER.debug("reactionStructureInformation.getEduct(substrateIndex).getAtomCount() " + reactionStructureInformation.getEduct(substrateIndex).getAtomCount());
                    LOGGER.debug("reactionStructureInformation.getProduct(productIndex).getAtomCount() " + reactionStructureInformation.getProduct(productIndex).getAtomCount());
                    if ((educt != null && product != null)
                            && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                            && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                            || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
//                        if (reactionStructureInformation.isEductModified(substrateIndex)
//                                || reactionStructureInformation.isProductModified(productIndex)) {

                        Combination c = new Combination(substrateIndex, productIndex);
                        jobReplicatorList.add(c);
//                        }
                    }
                }
            }

            LOGGER.debug("jobReplicatorList " + jobReplicatorList.size());

            if (jobReplicatorList.isEmpty()) {
                return unmodifiableCollection(mcsSolutions);
            }

            Map<Combination, Set<Combination>> jobMap = new TreeMap<>();

            for (Combination c : jobReplicatorList) {
                int substrateIndex = c.getRowIndex();
                int productIndex = c.getColIndex();
                IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                IAtomContainer product = reactionStructureInformation.getProduct(productIndex);

                boolean flag = false;
                for (Combination k : jobMap.keySet()) {
                    IAtomContainer eductJob = reactionStructureInformation.getEduct(k.getRowIndex());
                    IAtomContainer productJob = reactionStructureInformation.getProduct(k.getColIndex());

                    if (eductJob == educt
                            && productJob == product) {
                        if (eductJob.getAtomCount() == educt.getAtomCount()
                                && productJob.getAtomCount() == (product.getAtomCount())) {
                            jobMap.get(k).add(c);
                            flag = true;
                            break;
                        }
                    }
                }

                if (!flag) {
                    Set<Combination> set = new TreeSet<>();
                    jobMap.put(c, set);
                }
            }

            /*
             * Assign the threads
             *
             * Use Single Thread to computed MCS as muntiple threads lock the calculations!
             */
//            executor = newSingleThreadExecutor();
            int threadsAvailable = Math.max(1, getRuntime().availableProcessors() - 1);
            if (threadsAvailable > jobMap.size()) {
                threadsAvailable = jobMap.size();
            }

            LOGGER.debug(threadsAvailable + " threads requested for MCS in " + mh.getTheory());

//            executor = Executors.newSingleThreadExecutor();
//            executor = Executors.newCachedThreadPool();
            executor = Executors.newFixedThreadPool(threadsAvailable);
            CompletionService<MCSSolution> callablesQueue = new ExecutorCompletionService<>(executor);

            List<MCSThread> listOfJobs = new ArrayList<>();

            /*
             * Pre-compute aromaticity and cycle counts ONCE per molecule.
             * Previously this ran for every educt×product pair — O(E*P) redundancy.
             */
            CycleFinder allCycles = Cycles.or(Cycles.all(),
                    Cycles.or(Cycles.relevant(), Cycles.essential()));
            Aromaticity aromaticity = new Aromaticity(daylight(), allCycles);
            CycleFinder shortCycles = Cycles.vertexShort();

            Map<Integer, Integer> eductCycleCache = new TreeMap<>();
            for (int i = 0; i < eductCount; i++) {
                IAtomContainer educt = reactionStructureInformation.getEduct(i);
                if (educt != null && educt.getAtomCount() > 0) {
                    aromaticity.apply(educt);
                    eductCycleCache.put(i, shortCycles.find(educt).numberOfCycles());
                } else {
                    eductCycleCache.put(i, 0);
                }
            }

            Map<Integer, Integer> productCycleCache = new TreeMap<>();
            for (int j = 0; j < productCount; j++) {
                IAtomContainer product = reactionStructureInformation.getProduct(j);
                if (product != null && product.getAtomCount() > 0) {
                    aromaticity.apply(product);
                    productCycleCache.put(j, shortCycles.find(product).numberOfCycles());
                } else {
                    productCycleCache.put(j, 0);
                }
            }

            /*
             * Pre-compute canonical SMILES for identity detection.
             */
            SmilesGenerator canonSmigen = new SmilesGenerator(
                    SmiFlavor.Canonical | SmiFlavor.Stereo);
            Map<Integer, String> eductSmiles = new TreeMap<>();
            for (int i = 0; i < eductCount; i++) {
                IAtomContainer e = reactionStructureInformation.getEduct(i);
                if (e != null && e.getAtomCount() > 0) {
                    try { eductSmiles.put(i, canonSmigen.create(e)); }
                    catch (CDKException ex) { eductSmiles.put(i, ""); }
                }
            }
            Map<Integer, String> productSmiles = new TreeMap<>();
            for (int j = 0; j < productCount; j++) {
                IAtomContainer p = reactionStructureInformation.getProduct(j);
                if (p != null && p.getAtomCount() > 0) {
                    try { productSmiles.put(j, canonSmigen.create(p)); }
                    catch (CDKException ex) { productSmiles.put(j, ""); }
                }
            }

            int skippedIdentity = 0, skippedRatio = 0, skippedTanimoto = 0;

            for (Combination c : jobMap.keySet()) {
                int substrateIndex = c.getRowIndex();
                int productIndex = c.getColIndex();
                IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                IAtomContainer product = reactionStructureInformation.getProduct(productIndex);

                /*
                 * PRE-FILTER 1: Identity — if canonical SMILES match, the MCS
                 * will be the full molecule. Still run MCS (needed for correct
                 * atom correspondence), but mark as likely identity for logging.
                 */
                String eSmi = eductSmiles.getOrDefault(substrateIndex, "");
                String pSmi = productSmiles.getOrDefault(productIndex, "");
                if (!eSmi.isEmpty() && eSmi.equals(pSmi)) {
                    skippedIdentity++; // just count, still run MCS for correct mapping
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

                int numberOfCyclesEduct = eductCycleCache.getOrDefault(substrateIndex, 0);
                int numberOfCyclesProduct = productCycleCache.getOrDefault(productIndex, 0);
                boolean ringSizeEqual = (numberOfCyclesEduct == numberOfCyclesProduct);

                MCSThread mcsThread = new MCSThread(mh.getTheory(),
                        substrateIndex, productIndex, educt, product);
                mcsThread.setHasPerfectRings(ringSizeEqual);
                mcsThread.setEductRingCount(numberOfCyclesEduct);
                mcsThread.setProductRingCount(numberOfCyclesProduct);
                listOfJobs.add(mcsThread);
            }

            if (skippedIdentity + skippedRatio + skippedTanimoto > 0) {
                LOGGER.debug("Pre-filter: skipped " + skippedIdentity + " identity, "
                        + skippedRatio + " ratio, " + skippedTanimoto + " tanimoto pairs");
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
            for (int count = 0; count < taskCounter; count++) {
                MCSSolution isomorphism = callablesQueue.take().get();
                threadedUniqueMCSSolutions.add(isomorphism);
            }
            // This will make the executor accept no new threads
            // and finish all existing threads in the queue
            executor.shutdown();
            // Wait until all threads are finished
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

            LOGGER.debug("==Gathering MCS solution from the Thread==");
            threadedUniqueMCSSolutions.stream().filter((mcs) -> !(mcs == null)).map((MCSSolution mcs) -> {
                int queryPosition = mcs.getQueryPosition();
                int targetPosition = mcs.getTargetPosition();
                Combination referenceKey = null;
                for (Combination c : jobMap.keySet()) {
                    if (c.getRowIndex() == queryPosition && c.getColIndex() == targetPosition) {
                        referenceKey = c;
                        MCSSolution replicatedMCS = replicateMappingOnContainers(mh, c, mcs);
                        LOGGER.debug("======MCSSolution======");
                        LOGGER.debug("MCS " + " I " + queryPosition
                                + " J " + targetPosition
                                + " Number of Atom Mapped " + mcs.getAtomAtomMapping().getCount()
                                + " Number of Atom Mapped replicatedMCS " + replicatedMCS.getAtomAtomMapping().getCount());
                        mcsSolutions.add(replicatedMCS);
                    }
                }
                return referenceKey;
            }).filter((removeKey) -> (removeKey != null)).forEach((removeKey) -> {
                jobMap.remove(removeKey);
            });
            jobReplicatorList.clear();

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        } finally {
            if (executor != null) {
                executor.shutdown();
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

            // Build ID→Atom maps for O(1) lookup instead of O(n) linear scan
            Map<String, IAtom> qIdMap = new java.util.HashMap<>();
            for (IAtom a : q.atoms()) {
                if (a.getID() != null) {
                    qIdMap.put(a.getID(), a);
                }
            }
            Map<String, IAtom> tIdMap = new java.util.HashMap<>();
            for (IAtom a : t.atoms()) {
                if (a.getID() != null) {
                    tIdMap.put(a.getID(), a);
                }
            }

            atomAtomMapping.getMappingsByAtoms().forEach((a, b) -> {
                IAtom atomByID1 = a.getID() != null ? qIdMap.get(a.getID()) : null;
                IAtom atomByID2 = b.getID() != null ? tIdMap.get(b.getID()) : null;
                if (atomByID1 != null && atomByID2 != null) {
                    atomAtomMappingNew.put(atomByID1, atomByID2);
                } else {
                    LOGGER.error(WARNING, "UnExpected NULL ATOM FOUND");
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
            String a = this.row + "_" + this.col;
            String b = o.row + "_" + o.col;
            return a.compareTo(b);
        }

        @Override
        public int compare(Combination o1, Combination o2) {
            String a = o1.row + "_" + o1.col;
            String b = o2.row + "_" + o2.col;
            return a.compareTo(b);
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

            if (bestAtomMappingList != null) {
                for (Map.Entry<IAtom, IAtom> map : bestAtomMappingList.entrySet()) {
                    String eID = map.getKey().getID();
                    IAtom eAtom = getAtomByID(educt, eID);
                    String pID = map.getValue().getID();
                    LOGGER.debug("eID " + eID + ",pID " + pID);
                    IAtom pAtom = getAtomByID(product, pID);

                    if (eAtom != null && pAtom != null) {
                        IMapping im = SilentChemObjectBuilder.getInstance().newInstance(IMapping.class, eAtom, pAtom);
                        reaction.addMapping(im);
                    }
                    if (eAtom != null) {
                        educt.removeAtom(eAtom);
                    }
                    if (pAtom != null) {
                        product.removeAtom(pAtom);
                    }
                    delta = fragmentCount;
                }
            }

            for (IAtom atom : educt.atoms()) {
                IAtom matchedAtom = getAtomByID(matchedPart, atom.getID());
                if (matchedAtom != null) {
                    matchedPart.removeAtom(matchedAtom);
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

        private IAtom getAtomByID(IAtomContainer ac, String ID) {
            if (ID == null) {
                return null;
            }
            for (IAtom atom : ac.atoms()) {
                if (ID.equalsIgnoreCase(atom.getID())) {
                    return atom;
                }
            }
            return null;
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

        private final SmilesGenerator smiles;

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
                IAtomContainer educt, IAtomContainer product)
                throws CloneNotSupportedException, CDKException {
            /*
             * create SMILES
             */
            smiles = new SmilesGenerator(
                    //SmiFlavor.Unique|
                    SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.Stereo
                    | SmiFlavor.AtomAtomMap);

            this.compound1 = getNewContainerWithIDs(educt);
            this.compound2 = getNewContainerWithIDs(product);
            this.queryPosition = queryPosition;
            this.targetPosition = targetPosition;
            this.theory = theory;
            this.numberOfCyclesEduct = 0;
            this.numberOfCyclesProduct = 0;
        }

        void printMatch(BaseMapping isomorphism) {
            int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                    : isomorphism.getFirstAtomMapping().getCount();

            try {
                LOGGER.debug("Q: " + isomorphism.getQuery().getID()
                        + " T: " + isomorphism.getTarget().getID()
                        + " atoms: " + isomorphism.getQuery().getAtomCount()
                        + " atoms: " + isomorphism.getTarget().getAtomCount()
                        + " overlaps: " + overlap
                        + " mcs " + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
            } catch (CloneNotSupportedException | CDKException ex) {
                LOGGER.error(Level.SEVERE, "Print MCS ", ex.getMessage());
            }
        }

        @Override
        public MCSSolution call() throws Exception {
            boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

            AtomMatcher am;
            BondMatcher bm;

            LOGGER.debug("in mcsthread call ");

            try {
                /*
                 * IMP: Do not perform substructure matching for disconnected molecules
                 */
                boolean moleculeConnected = isMoleculeConnected(getCompound1(), getCompound2());
                /*
                     Check if MCS matching required or not very IMP step
                 */
                boolean possibleVFmatch12 = isPossibleSubgraphMatch(getCompound1(), getCompound2());
                LOGGER.debug("VF Matcher 1->2 " + possibleVFmatch12);

                boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
                LOGGER.debug("VF Matcher 2->1 " + possibleVFmatch21);

                if (moleculeConnected && possibleVFmatch12) {
                    LOGGER.debug("Substructure 1");
                    this.startTime = currentTimeMillis();

                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());

                    LOGGER.debug("---1.1---");
                    Substructure substructure;
                    am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    substructure = new Substructure(ac1, ac2, am, bm, true);

                    if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                        am = AtomBondMatcher.atomMatcher(false, ringFlag);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---1.3---");
                        substructure = new Substructure(ac1, ac2,
                                am, bm, true);
                    } else if (moleculeConnected && !substructure.isSubgraph()) {
                        am = AtomBondMatcher.atomMatcher(false, false);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---1.2---");
                        substructure = new Substructure(ac1, ac2, am, bm, true);
                    }
                    substructure.setChemFilters(true, true, true);
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
                    Substructure substructure;

                    LOGGER.debug("---2.1---");
                    am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    substructure = new Substructure(ac2, ac1, am, bm, true);

                    if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                        am = AtomBondMatcher.atomMatcher(false, ringFlag);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---2.3---");
                        substructure = new Substructure(ac2, ac1, am, bm, true);
                    } else if (moleculeConnected && !substructure.isSubgraph()) {
                        am = AtomBondMatcher.atomMatcher(false, false);
                        bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                        LOGGER.debug("---2.2---");
                        substructure = new Substructure(ac2, ac1, am, bm, true);
                    }
                    substructure.setChemFilters(true, true, true);

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

            } catch (CDKException | CloneNotSupportedException ex) {
                LOGGER.error(SEVERE, "Error in generating MCS Solution: ", ex.getMessage());
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

                return ac;
            }
            return mol;
        }

        private boolean isPossibleSubgraphMatch(IAtomContainer q, IAtomContainer t) {
            LOGGER.debug("check isPossibleSubgraphMatch " + q.getID() + "," + t.getID());
            Map<String, Integer> atomCount1 = new HashMap<>();
            Map<String, Integer> atomCount2 = new HashMap<>();

            for (IAtom a : q.atoms()) {
                atomCount1.merge(a.getSymbol(), 1, Integer::sum);
            }

            for (IAtom b : t.atoms()) {
                atomCount2.merge(b.getSymbol(), 1, Integer::sum);
            }

            if (atomCount1.size() > atomCount2.size()) {
                return false;
            }

            // Check all atom types in query exist in target with sufficient count
            for (Map.Entry<String, Integer> entry : atomCount1.entrySet()) {
                Integer targetCount = atomCount2.get(entry.getKey());
                if (targetCount == null || entry.getValue() > targetCount) {
                    return false;
                }
            }

            return true;
        }

        private int expectedMaxGraphmatch(IAtomContainer q, IAtomContainer t) {
            /*
             a={c,c,c,o,n}
             b={c,c,c,p}
             expectedMaxGraphmatch=3;
             */
            Map<String, Integer> countQ = new HashMap<>();
            Map<String, Integer> countT = new HashMap<>();

            for (IAtom a : q.atoms()) {
                String hyb = a.getHybridization() == UNSET
                        ? a.getSymbol() : a.getAtomTypeName();
                countQ.merge(hyb, 1, Integer::sum);
            }

            for (IAtom b : t.atoms()) {
                String hyb = b.getHybridization() == UNSET
                        ? b.getSymbol() : b.getAtomTypeName();
                countT.merge(hyb, 1, Integer::sum);
            }

            if (countQ.isEmpty()) {
                return 0;
            }

            // Multiset intersection: min of counts for each common type
            int common = 0;
            for (Map.Entry<String, Integer> entry : countQ.entrySet()) {
                Integer tCount = countT.get(entry.getKey());
                if (tCount != null) {
                    common += Math.min(entry.getValue(), tCount);
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
            Isomorphism isomorphism;
            int expectedMaxGraphmatch = expectedMaxGraphmatch(ac1, ac2);
            boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

            LOGGER.debug("Expected matches " + expectedMaxGraphmatch);

            String key;
            MCSSolution mcs;
            AtomMatcher am;
            BondMatcher bm;
            boolean atomType;
            boolean bondMatch;
            boolean ringMatch;
            boolean ringSizeMatch;

            switch (theory) {
                case RINGS:

                    atomType = false;
                    bondMatch = false;
                    ringMatch = ringFlag;
                    ringSizeMatch = isHasPerfectRings();

                    break;

                case MIN:

                    atomType = false;
                    bondMatch = false;
                    ringMatch = isHasPerfectRings();
                    ringSizeMatch = false;

                    break;

                case MAX:

                    atomType = false;
                    bondMatch = true;
                    ringMatch = isHasPerfectRings();
                    ringSizeMatch = false;

                    break;
                default:

                    atomType = false;
                    bondMatch = false;
                    ringMatch = isHasPerfectRings();
                    ringSizeMatch = false;

                    break;
            }

            am = AtomBondMatcher.atomMatcher(atomType, ringSizeMatch);
            bm = AtomBondMatcher.bondMatcher(bondMatch, ringMatch);

            key = generateUniqueKey(getCompound1().getID(), getCompound2().getID(),
                    compound1.getAtomCount(), compound2.getAtomCount(),
                    compound1.getBondCount(), compound2.getBondCount(),
                    atomType,
                    bondMatch,
                    ringMatch,
                    ringSizeMatch,
                    numberOfCyclesEduct,
                    numberOfCyclesProduct
            );
            if (ThreadSafeCache.getInstance().containsKey(key)) {
                LOGGER.debug("===={Aladdin} Mapping {Gini}====");
                MCSSolution solution = (MCSSolution) ThreadSafeCache.getInstance().get(key);
                mcs = copyOldSolutionToNew(
                        getQueryPosition(), getTargetPosition(),
                        getCompound1(), getCompound2(),
                        solution);

            } else {
                isomorphism = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, am, bm);
                mcs = addMCSSolution(key, ThreadSafeCache.getInstance(), isomorphism);
            }

            return mcs;

        }

        private IAtomContainer duplicate(IAtomContainer ac) throws CloneNotSupportedException {
            IAtomContainer a = ac.clone();
            a.setID(ac.getID());
            a.setProperties(ac.getProperties());

            for (int i = 0; i < a.getAtomCount(); i++) {
                a.getAtom(i).setID(ac.getAtom(i).getID());
            }

            // Fix aromatic bond consistency: if a bond is aromatic but its atoms
            // are not flagged aromatic, downgrade the bond to SINGLE to prevent
            // "Aromatic bond connects non-aromatic atoms" errors in SMSD
            for (IBond bond : a.bonds()) {
                if (bond.isAromatic()) {
                    IAtom begin = bond.getBegin();
                    IAtom end = bond.getEnd();
                    if ((begin != null && !begin.isAromatic())
                            || (end != null && !end.isAromatic())) {
                        bond.setIsAromatic(false);
                        if (bond.getOrder() == null || bond.getOrder() == IBond.Order.UNSET) {
                            bond.setOrder(IBond.Order.SINGLE);
                        }
                    }
                }
            }

            return a;
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
        private boolean isMoleculeConnected(IAtomContainer compound1, IAtomContainer compound2) {
            LOGGER.debug("isMoleculeConnected");
            boolean connected1 = true;

            IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound1);
            if (partitionIntoMolecules.getAtomContainerCount() > 1) {
                connected1 = false;
            }

            boolean connected2 = true;

            partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound2);
            if (partitionIntoMolecules.getAtomContainerCount() > 1) {
                connected2 = false;
            }

            return connected1 && connected2;
        }

        void setEductRingCount(int numberOfCyclesEduct) {
            this.numberOfCyclesEduct = numberOfCyclesEduct;
        }

        void setProductRingCount(int numberOfCyclesProduct) {
            this.numberOfCyclesProduct = numberOfCyclesProduct;
        }

        private static final String CACHED_SMILES = "CACHED_CANONICAL_SMILES";
        private static final SmilesGenerator CANONICAL_SMIGEN = new SmilesGenerator(
                SmiFlavor.Canonical | SmiFlavor.Stereo);

        /**
         * Generate a unique cache key based on canonical SMILES (structure-based,
         * not ID-based). This enables cross-reaction cache hits: if the same
         * molecule pair appears in different reactions, the MCS result is reused.
         */
        String generateUniqueKey(
                String id1, String id2,
                int atomCount1, int atomCount2,
                int bondCount1, int bondCount2,
                boolean atomtypeMatcher,
                boolean bondMatcher,
                boolean ringMatcher,
                boolean hasPerfectRings,
                int numberOfCyclesEduct, int numberOfCyclesProduct) {

            StringBuilder key = new StringBuilder();

            // Use canonical SMILES as the molecular identity (not mol IDs)
            String smi1 = getCanonicalSmiles(compound1);
            String smi2 = getCanonicalSmiles(compound2);
            key.append(smi1).append(">>").append(smi2);

            // Append matcher flags that affect the MCS result
            key.append('|')
                    .append(atomtypeMatcher ? '1' : '0')
                    .append(bondMatcher ? '1' : '0')
                    .append(ringMatcher ? '1' : '0')
                    .append(hasPerfectRings ? '1' : '0');

            return key.toString();
        }

        /**
         * Get or compute canonical SMILES for a molecule. Cached on the molecule
         * to avoid recomputation across calls within the same reaction.
         */
        private String getCanonicalSmiles(IAtomContainer mol) {
            String cached = mol.getProperty(CACHED_SMILES);
            if (cached != null) {
                return cached;
            }
            try {
                cached = CANONICAL_SMIGEN.create(mol);
            } catch (CDKException e) {
                // Fallback: use atom/bond counts + fingerprint hash
                cached = mol.getAtomCount() + ":" + mol.getBondCount() + ":" + mol.hashCode();
            }
            mol.setProperty(CACHED_SMILES, cached);
            return cached;
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

        MCSSolution addMCSSolution(String key, ThreadSafeCache<String, MCSSolution> mappingcache, Isomorphism isomorphism) {

            isomorphism.setChemFilters(true, true, true);
            try {
                LOGGER.debug("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", "
                        + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
            } catch (CloneNotSupportedException | CDKException e) {
                LOGGER.error(SEVERE, "Error in computing MCS ", e.getMessage());
            }
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
            if (!mappingcache.containsKey(key)) {
                LOGGER.debug("Key " + key);
                try {
                    LOGGER.debug("mcs size " + mcs.getAtomAtomMapping().getCount());
                    LOGGER.debug("mcs map " + mcs.getAtomAtomMapping().getMappingsByIndex());
                    LOGGER.debug("mcs " + mcs.getAtomAtomMapping().getCommonFragmentAsSMILES());
                    LOGGER.debug("\n\n\n ");
                } catch (CloneNotSupportedException | CDKException ex) {
                    LOGGER.error(SEVERE, "Unable to create SMILES ", ex.getMessage());
                }
                mappingcache.put(key, mcs);
            }
            return mcs;
        }
    }


}
