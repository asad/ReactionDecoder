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
package uk.ac.ebi.reactionblast.mapping.graph;

import java.io.IOException;
import static java.lang.System.gc;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.unmodifiableCollection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.helper.Debugger;
import static java.util.Collections.synchronizedCollection;
import static java.util.concurrent.Executors.newSingleThreadExecutor;

import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.smiles.SmiFlavor;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class GraphMatcher extends Debugger {

    static final String NEW_LINE = getProperty("line.separator");
    private final static boolean DEBUG = false;
    private final static ILoggingTool LOGGER
            = createLoggingTool(GraphMatcher.class);

    /**
     *
     * @param mh
     * @return
     * @throws InterruptedException
     */
    public synchronized static Collection<MCSSolution> matcher(Holder mh) throws InterruptedException {
        ExecutorService executor = null;
        Collection<MCSSolution> mcsSolutions = synchronizedCollection(new ArrayList<>());

        if (DEBUG) {
            System.out.println("Matcher Class for " + mh.getTheory());
        }
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
                    if ((educt != null && product != null)
                            && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                            && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                            || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
                        if (reactionStructureInformation.isEductModified(substrateIndex)
                                || reactionStructureInformation.isProductModified(productIndex)) {

                            Combination c = new Combination(substrateIndex, productIndex);
                            jobReplicatorList.add(c);
                        }
                    }
                }
            }

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
            executor = newSingleThreadExecutor();

//            int threadsAvailable = getRuntime().availableProcessors() - 1;
//            if (threadsAvailable == 0) {
//                threadsAvailable = 1;
//            }
//
//            if (threadsAvailable > jobMap.size()) {
//                threadsAvailable = jobMap.size();
//            }
//            if (DEBUG) {
//                out.println(threadsAvailable + " threads requested for MCS in " + mh.getTheory());
//            }
//
//            if (DEBUG) {
//                executor = newSingleThreadExecutor();
//            } else {
//                executor = Executors.newFixedThreadPool(threadsAvailable);
//            }
            CompletionService<MCSSolution> callablesQueue = new ExecutorCompletionService<>(executor);

            for (Combination c : jobMap.keySet()) {
                int substrateIndex = c.getRowIndex();
                int productIndex = c.getColIndex();
                IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                IAtomContainer product = reactionStructureInformation.getProduct(productIndex);

                /*
                 Ring matcher is set true if both sides have rings else it set to false (IMP for MCS)
                 */
                boolean ring = false;
                boolean ringSizeEqual = false;

                /*
                 * Report All Cycles
                 * or 
                 * CycleFinder cycles = or(Cycles.all(), Cycles.all());
                 * CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.relevant());
                 */
                Aromaticity aromaticity = new Aromaticity(daylight(),
                        Cycles.or(Cycles.all(),
                                Cycles.or(Cycles.relevant(),
                                        Cycles.essential())));
                /*
                 * Aromatise molecule for escaping CDKtoBeam Aromatic bond error
                 */
                aromaticity.apply(educt);
                aromaticity.apply(product);

                /*
                 * Report short cycyles
                 */
                CycleFinder cycles = Cycles.vertexShort();
                Cycles rings = cycles.find(educt);
                int numberOfCyclesEduct = rings.numberOfCycles();
                rings = cycles.find(product);
                int numberOfCyclesProduct = rings.numberOfCycles();
                if (numberOfCyclesEduct > 0 && numberOfCyclesProduct > 0) {
                    ring = true;
                }

                if (numberOfCyclesEduct == numberOfCyclesProduct) {
                    ringSizeEqual = true;
                }
                if (DEBUG) {
                    try {
                        SmilesGenerator smilesGenerator;
                        System.out.println("SMILES");
                        smilesGenerator = new SmilesGenerator(SmiFlavor.Stereo
                                | SmiFlavor.AtomAtomMap);
                        out.println(educt.getID() + " ED: " + smilesGenerator.create(educt));
                        out.println(product.getID() + " PD: " + smilesGenerator.create(product));
                        out.println("numberOfCyclesEduct " + numberOfCyclesEduct);
                        out.println("numberOfCyclesProduct " + numberOfCyclesProduct);
                        out.println("ringSizeEqual " + ringSizeEqual);
                        out.println("Ring " + ring);
                        out.println("----------------------------------");
                    } catch (CDKException e) {
                        LOGGER.error(SEVERE, null, e);
                    }
                }

                MCSThread mcsThread;
                switch (mh.getTheory()) {

                    case MIN:
                        mcsThread = new MCSThread(mh.getTheory(), substrateIndex, productIndex, educt, product, false, ring, true);
                        mcsThread.setHasPerfectRings(ringSizeEqual);
                        mcsThread.setEductRingCount(numberOfCyclesEduct);
                        mcsThread.setProductRingCount(numberOfCyclesProduct);
                        
                        break;

                    case MAX:
                        mcsThread = new MCSThread(mh.getTheory(), substrateIndex, productIndex, educt, product, false, ring, true);
                        mcsThread.setHasPerfectRings(ringSizeEqual);
                        mcsThread.setEductRingCount(numberOfCyclesEduct);
                        mcsThread.setProductRingCount(numberOfCyclesProduct);
                        break;

                    case MIXTURE:
                        mcsThread = new MCSThread(mh.getTheory(), substrateIndex, productIndex, educt, product, false, ring, false);
                        mcsThread.setHasPerfectRings(ringSizeEqual);
                        mcsThread.setEductRingCount(numberOfCyclesEduct);
                        mcsThread.setProductRingCount(numberOfCyclesProduct);
                        break;

                    case RINGS:
                        /*
                         * don't use ring matcher if there are no rings in the molecule
                         * else mappings with be skewed
                         * bond=false;
                         * ring =true;
                         * atom type=true;
                         * Ex: R05219
                         */
                        mcsThread = new MCSThread(mh.getTheory(), substrateIndex, productIndex, educt, product, false, ring, true);
                        mcsThread.setHasPerfectRings(ringSizeEqual);
                        mcsThread.setEductRingCount(numberOfCyclesEduct);
                        mcsThread.setProductRingCount(numberOfCyclesProduct);
                        break;

                    default:
                        mcsThread = null;
                        break;
                }
                if (mcsThread != null) {
                    callablesQueue.submit(mcsThread);
                    taskCounter++;
                }
            }

            Collection<MCSSolution> threadedUniqueMCSSolutions = synchronizedCollection(new ArrayList<>());
            for (int count = 0; count < taskCounter; count++) {
                MCSSolution isomorphism = callablesQueue.take().get();
                threadedUniqueMCSSolutions.add(isomorphism);
            }

//                List<Future<MCSSolution>> invokeAll = executor.invokeAll(callablesQueue);
//                for (Iterator<Future<MCSSolution>> it = invokeAll.iterator(); it.hasNext();) {
//                    Future<MCSSolution> callable = it.next();
//                    MCSSolution isomorphism = callable.get();
//                    if (callable.isDone()) {
//                        mcsSolutions.add(isomorphism);
//                    }
//                }
            // This will make the executor accept no new threads
            // and finish all existing threads in the queue
            executor.shutdown();
            // Wait until all threads are finish
            while (!executor.isTerminated()) {
            }

            if (DEBUG) {
                out.println("==Gathering MCS solution from the Thread==");
            }
            threadedUniqueMCSSolutions.stream().filter((mcs) -> !(mcs == null)).map((MCSSolution mcs) -> {
                int queryPosition = mcs.getQueryPosition();
                int targetPosition = mcs.getTargetPosition();
//                if (DEBUG) {
//                    System.out.println("");
//                    out.println("MCS " + " I " + queryPosition
//                            + " J " + targetPosition
//                            + " Number of Atom Mapped " + mcs.getAtomAtomMapping().getCount());
//                }
                Combination referenceKey = null;
                for (Combination c : jobMap.keySet()) {
                    if (c.getRowIndex() == queryPosition && c.getColIndex() == targetPosition) {
                        referenceKey = c;
                        MCSSolution replicatedMCS = replicateMappingOnContainers(mh, c, mcs);
                        if (DEBUG) {
                            System.out.println("======MCSSolution======");
                            out.println("MCS " + " I " + queryPosition
                                    + " J " + targetPosition
                                    + " Number of Atom Mapped " + mcs.getAtomAtomMapping().getCount()
                                    + " Number of Atom Mapped replicatedMCS " + replicatedMCS.getAtomAtomMapping().getCount());
                        }
                        mcsSolutions.add(replicatedMCS);
                    }
                }
                return referenceKey;
            }).filter((removeKey) -> (removeKey != null)).forEach((removeKey) -> {
                jobMap.remove(removeKey);
            });
            jobReplicatorList.clear();
            gc();

        } catch (IOException | CDKException | ExecutionException | InterruptedException | CloneNotSupportedException ex) {
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

            if (DEBUG) {
                if (diff1 != 0 && diff2 != 0) {

                    out.println(NEW_LINE + NEW_LINE + " " + solution.getRowIndex() + ", Diff in ac1 " + diff1);
                    out.println(solution.getColIndex() + ", Diff in ac2 " + diff2);
                    out.println(NEW_LINE + "ac1 " + q.getAtomCount());
                    out.println(NEW_LINE + "ac2 " + t.getAtomCount());

                    out.println(NEW_LINE + "mac1 " + mcs.getQueryContainer().getAtomCount());
                    out.println(NEW_LINE + "mac2 " + mcs.getTargetContainer().getAtomCount());
                }
            }

            AtomAtomMapping atomAtomMapping = mcs.getAtomAtomMapping();
            AtomAtomMapping atomAtomMappingNew = new AtomAtomMapping(q, t);
            atomAtomMapping.getMappingsByAtoms().keySet().stream().forEach((a) -> {
                IAtom atomByID1 = getAtomByID(q, a);
                IAtom b = atomAtomMapping.getMappingsByAtoms().get(a);
                IAtom atomByID2 = getAtomByID(t, b);
//                if (DEBUG) {
//                    out.println("atomByID1 " + atomByID1.getID() + " atomByID2 " + atomByID2.getID());
//                }
                if (atomByID1 != null && atomByID2 != null) {
                    atomAtomMappingNew.put(atomByID1, atomByID2);
                } else {
                    LOGGER.error(WARNING, "UnExpected NULL ATOM FOUND");
                    System.err.println("WARNING: " + "UnExpected NULL ATOM FOUND");
                }
            });

            if (DEBUG) {
                System.out.println("------Mapped PAIRS------");
                System.out.println("Query " + q.getAtomCount());
                System.out.println("Target " + t.getAtomCount());
                System.out.println("Mapping Size " + atomAtomMappingNew.getCount());
            }
            return new MCSSolution(solution.getRowIndex(), solution.getColIndex(), q, t, atomAtomMappingNew);
        } catch (IOException | CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return null;
    }

    private static IAtom getAtomByID(IAtomContainer ac, IAtom atom) {
        if (atom.getID() == null) {
            return null;
        }
        for (IAtom a : ac.atoms()) {
            if (a.getID().equals(atom.getID())) {
                return a;
            }
        }
        return null;
    }

    /**
     *
     * @param educt
     * @param product
     * @param matchesBond bond type matching
     * @param ringMatcher ring matcher
     * @return
     * @throws CDKException
     */
    static boolean isSubgraph(IAtomContainer educt, IAtomContainer product, boolean matchesBond, boolean ringMatcher) throws CDKException {
        if (educt.getAtomCount() <= product.getAtomCount()) {
            Substructure smsd = new Substructure(educt, product, matchesBond, ringMatcher, false, false);
            return smsd.isSubgraph();
        } else {
            Substructure smsd = new Substructure(product, educt, matchesBond, ringMatcher, false, false);
            return smsd.isSubgraph();
        }
    }
}
