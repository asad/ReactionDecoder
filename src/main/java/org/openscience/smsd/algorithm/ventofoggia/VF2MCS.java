/* Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.ventofoggia;

import java.io.IOException;
import static java.lang.Runtime.getRuntime;
import java.util.*;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.graph.algorithm.VentoFoggia;
import org.openscience.smsd.helper.Mappings;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib
 * {@link org.openscience.smsd.algorithm.ventofoggia12.VF2MCS} and reports MCS
 * between run source and target graphs. Then these solutions are extended using
 * McGregor {@link org.openscience.smsd.algorithm.mcgregor.McGregor} algorithm
 * where ever required.
 *
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public final class VF2MCS extends BaseMCS implements IResults {

    private final List<AtomAtomMapping> allAtomMCS;
    private final static ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(VF2MCS.class);
    private final boolean DEBUG = false;

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     * @param am
     * @param bm
     * @throws org.openscience.cdk.exception.CDKException
     */
    public VF2MCS(IAtomContainer source,
            IAtomContainer target,
            AtomMatcher am, BondMatcher bm) throws CDKException {
        super(source, target, am, bm);
        boolean timeoutVF = searchVFCDKMappings();

        if (DEBUG) {
            System.out.println("time for VF search " + timeoutVF);
        }

        /*
         * An extension is triggered if its mcs solution is smaller than reactant and product. An enrichment is
         * triggered if its mcs solution is equal to reactant or product size.
         *
         *
         */
        int maxVFMappingSize = allLocalMCS.iterator().hasNext()
                ? allLocalMCS.iterator().next().size() : 0;
        if (DEBUG) {
            System.out.println("maxVFMappingSize " + maxVFMappingSize);
        }

        /*
         * Atleast two atoms are unmapped else you will get bug due to unmapped single atoms
         */
        if (timeoutVF || (maxVFMappingSize != (source.getAtomCount())
                && maxVFMappingSize != (target.getAtomCount()))) {

            List<Map<Integer, Integer>> mcsVFSeeds = new ArrayList<>();

            /*
             * Copy VF based MCS solution in the seed
             */
            int counter = 0;
            for (Map<Integer, Integer> vfMapping : allLocalMCS) {
                mcsVFSeeds.add(counter, vfMapping);
                counter++;
            }

            /*
             * Clean VF mapping data
             */
            allLocalMCS.clear();
            allLocalAtomAtomMapping.clear();

            long startTimeSeeds = System.nanoTime();
            /*
             *   Assign the threads
             */
            int threadsAvailable = getRuntime().availableProcessors() - 1;
            if (threadsAvailable == 0) {
                threadsAvailable = 1;
            } else if (threadsAvailable > 2) {
                threadsAvailable = 2;
            }

//            ExecutorService executor = Executors.newCachedThreadPool();
            ExecutorService executor = Executors.newFixedThreadPool(threadsAvailable);
//            ExecutorService executor = Executors.newSingleThreadExecutor();
            CompletionService<List<AtomAtomMapping>> cs = new ExecutorCompletionService<>(executor);

            /*
             * Reduce the target size by removing bonds which do not share 
             * similar Hybridization 
             */
            IAtomContainer targetClone = null;
            try {
                targetClone = target.clone();
                Set<IBond> bondRemovedT = new HashSet<>();
                for (IBond b1 : targetClone.bonds()) {
                    boolean flag = false;
                    for (IBond b2 : source.bonds()) {
                        if (AtomBondMatcher.matchAtomAndBond(b1, b2, atomMatcher, bondMatcher, true)) {
                            flag = true;
                            break;
                        }
                    }
                    if (!flag) {
                        bondRemovedT.add(b1);
                    }
                }

                if (DEBUG) {
                    System.out.println("Bond to be removed " + bondRemovedT.size());
                }
                for (IBond b : bondRemovedT) {
                    targetClone.removeBond(b);
                }

            } catch (CloneNotSupportedException ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }

            /*
             * CDK MCS faulter on disconnected molecules
             */
            //boolean moleculeConnected = isMoleculeConnected(source, targetClone);
            int jobCounter = 0;

            if (targetClone != null) {
                if (source.getBondCount() > 0
                        && targetClone.getBondCount() > 0) {
                    if (DEBUG) {
                        System.out.println(" CALLING UIT ");
                    }
                    MCSSeedGenerator mcsSeedGeneratorUIT
                            = new MCSSeedGenerator(source, targetClone,
                                    Algorithm.CDKMCS, atomMatcher, bondMatcher);
                    cs.submit(mcsSeedGeneratorUIT);
                    jobCounter++;
                }
            }

            if (DEBUG) {
                System.out.println(" CALLING MCSPLUS ");
            }
            MCSSeedGenerator mcsSeedGeneratorKoch
                    = new MCSSeedGenerator(source, targetClone,
                            Algorithm.MCSPlus, atomMatcher, bondMatcher);
            cs.submit(mcsSeedGeneratorKoch);
            jobCounter++;

            /*
             * Generate the UIT based MCS seeds
             */
            Set<Map<Integer, Integer>> mcsSeeds = new HashSet<>();
            /*
             * Collect the results
             */
            for (int i = 0; i < jobCounter; i++) {
                List<AtomAtomMapping> chosen;
                try {
                    chosen = cs.take().get();
                    chosen.stream().map((mapping) -> {
                        Map<Integer, Integer> map = new TreeMap<>();
                        map.putAll(mapping.getMappingsByIndex());
                        return map;
                    }).forEach((map) -> {
                        mcsSeeds.add(map);
                    });
                } catch (Exception ex) {
                    if (DEBUG) {
                        ex.printStackTrace();
                    }
                    LOGGER.error(Level.SEVERE, null, ex);
                }
            }
            executor.shutdown();
            /*
             Wait until all threads are finish
             */

            while (!executor.isTerminated()) {
            }
            System.gc();

            long stopTimeSeeds = System.nanoTime();
            if (DEBUG) {
                System.out.println("time taken for seeds: "
                        + TimeUnit.MILLISECONDS.convert((stopTimeSeeds - startTimeSeeds),
                                TimeUnit.NANOSECONDS) + " ms.");
            }
            /*
             * Store largest MCS seeds generated from MCSPlus and UIT
             */
            int solutionSize = 0;
            counter = 0;
            List<Map<Integer, Integer>> cleanedMCSSeeds = new ArrayList<>();

            if (DEBUG) {
                System.out.println("merging  UIT & KochCliques");
            }

            if (!mcsSeeds.isEmpty()) {
                for (Map<Integer, Integer> map : mcsSeeds) {
                    if (DEBUG) {
                        System.out.println("potential seed MCSPlus, UIT " + map.size());
                    }
                    if (map.size() > solutionSize) {
                        solutionSize = map.size();
                        cleanedMCSSeeds.clear();
                        counter = 0;
                    }
                    if (!map.isEmpty()
                            && map.size() == solutionSize
                            && !super.isCliquePresent(map, cleanedMCSSeeds)) {
                        if (DEBUG) {
                            System.out.println("seed MCS, UIT " + cleanedMCSSeeds.size());
                        }
                        cleanedMCSSeeds.add(counter, map);
                        counter++;
                    }
                }
            }

            /*
             * Add seeds from VF MCS
             */
            mcsVFSeeds.stream().filter((map) -> (!map.isEmpty()
                    && !super.isCliquePresent(map, cleanedMCSSeeds))).forEach((_item) -> {
                cleanedMCSSeeds.addAll(mcsVFSeeds);
            });
            /*
             * Sort biggest clique to smallest
             */
            Collections.sort(cleanedMCSSeeds, new Map1ValueComparator(SortOrder.DESCENDING));

            /*
             * Extend the seeds using McGregor
             */
            try {
                super.extendCliquesWithMcGregor(cleanedMCSSeeds);
            } catch (CDKException | IOException ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }

            /*
             * Clear previous seeds
             */
            mcsSeeds.clear();
            cleanedMCSSeeds.clear();

            /*
             * Integerate the solutions
             */
            solutionSize = 0;
            counter = 0;
            this.allAtomMCS = new ArrayList<>();

            /*
             * Store solutions from VF MCS only
             */
            if (!allLocalAtomAtomMapping.isEmpty()) {
                for (AtomAtomMapping atomMCSMap : allLocalAtomAtomMapping) {
                    if (atomMCSMap.getCount() > solutionSize) {
                        solutionSize = atomMCSMap.getCount();
                        allAtomMCS.clear();
                        counter = 0;
                    }
                    if (!atomMCSMap.isEmpty()
                            && atomMCSMap.getCount() == solutionSize) {
                        allAtomMCS.add(counter, atomMCSMap);
                        counter++;
                    }
                }
            }

            /*
             * Clear the local solution after storing it into mcs solutions
             */
            allLocalMCS.clear();
            allLocalAtomAtomMapping.clear();

        } else {
            if (DEBUG) {
                System.out.println("IS A Subgraph ");
            }
            /*
             * Store solutions from VF MCS only
             */
            int solutionSize = 0;
            int counter = 0;
            this.allAtomMCS = new ArrayList<>();
            /*
             * Store solutions from VF MCS only
             */
            if (!allLocalAtomAtomMapping.isEmpty()) {
                for (AtomAtomMapping atomMCSMap : allLocalAtomAtomMapping) {
                    if (atomMCSMap.getCount() > solutionSize) {
                        solutionSize = atomMCSMap.getCount();
                        allAtomMCS.clear();
                        counter = 0;
                    }
                    if (!atomMCSMap.isEmpty()
                            && atomMCSMap.getCount() == solutionSize) {
                        allAtomMCS.add(counter, atomMCSMap);
                        counter++;
                    }
                }
            }

            /*
             * Clear the local solution after storing it into mcs solutions
             */
            allLocalMCS.clear();
            allLocalAtomAtomMapping.clear();
        }
    }

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     * @param am
     * @param bm
     */
    public VF2MCS(IQueryAtomContainer source, IAtomContainer target, AtomMatcher am, BondMatcher bm) {
        super((IQueryAtomContainer) source, target, am, bm);
        boolean timeoutVF = searchVFCDKMappings();

//        System.out.println("time for VF search " + timeoutVF);

        /*
         * An extension is triggered if its mcs solution is smaller than reactant and product. An enrichment is
         * triggered if its mcs solution is equal to reactant or product size.
         *
         *
         */
        if (!timeoutVF) {

            List<Map<Integer, Integer>> mcsVFSeeds = new ArrayList<>();

            /*
             * Copy VF based MCS solution in the seed
             */
            int counter = 0;
            for (Map<Integer, Integer> vfMapping : allLocalMCS) {
                mcsVFSeeds.add(counter, vfMapping);
                counter++;
            }

            /*
             * Clean VF mapping data
             */
            allLocalMCS.clear();
            allLocalAtomAtomMapping.clear();

            long startTimeSeeds = System.nanoTime();

            ExecutorService executor = Executors.newCachedThreadPool();
            CompletionService<List<AtomAtomMapping>> cs = new ExecutorCompletionService<>(executor);

            /*
             * Reduce the target size by removing bonds which do not share 
             * similar Hybridization 
             */
            IAtomContainer targetClone = null;
            try {
                targetClone = target.clone();
                Set<IBond> bondRemovedT = new HashSet<>();
                for (IBond b1 : source.bonds()) {
                    IQueryBond bond = (IQueryBond) b1;
                    IQueryAtom a1 = (IQueryAtom) b1.getAtom(0);
                    IQueryAtom a2 = (IQueryAtom) b1.getAtom(1);
                    for (IBond b2 : targetClone.bonds()) {
                        boolean matches = bond.matches(b2);
                        if (a1.matches(b2.getAtom(0)) && a2.matches(b2.getAtom(1)) && !matches) {
                            bondRemovedT.add(b2);
                        } else if (a2.matches(b2.getAtom(0)) && a1.matches(b2.getAtom(1)) && !matches) {
                            bondRemovedT.add(b2);
                        }
                    }
                }

//                System.out.println("Bond to be removed " + bondRemovedT.size());
                for (IBond b : bondRemovedT) {
                    targetClone.removeBond(b);
                }

            } catch (CloneNotSupportedException ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }

            MCSSeedGenerator mcsSeedGeneratorUIT
                    = new MCSSeedGenerator((IQueryAtomContainer) source, targetClone, Algorithm.CDKMCS);
            MCSSeedGenerator mcsSeedGeneratorKoch
                    = new MCSSeedGenerator((IQueryAtomContainer) source, targetClone, Algorithm.MCSPlus);

            int jobCounter = 0;
            cs.submit(mcsSeedGeneratorUIT);
            jobCounter++;
            cs.submit(mcsSeedGeneratorKoch);
            jobCounter++;

            /*
             * Generate the UIT based MCS seeds
             */
            Set<Map<Integer, Integer>> mcsSeeds = new HashSet<>();
            /*
             * Collect the results
             */
            for (int i = 0; i < jobCounter; i++) {
                List<AtomAtomMapping> chosen;
                try {
                    chosen = cs.take().get();
                    chosen.stream().map((mapping) -> {
                        Map<Integer, Integer> map = new TreeMap<>();
                        map.putAll(mapping.getMappingsByIndex());
                        return map;
                    }).forEach((map) -> {
                        mcsSeeds.add(map);
                    });
                } catch (InterruptedException | ExecutionException ex) {
                    LOGGER.error(Level.SEVERE, null, ex);
                }
            }
            executor.shutdown();
            // Wait until all threads are finish
            while (!executor.isTerminated()) {
            }
            System.gc();

//            long stopTimeSeeds = System.nanoTime();
//            System.out.println("done seeds " + (stopTimeSeeds - startTimeSeeds));
            /*
             * Store largest MCS seeds generated from MCSPlus and UIT
             */
            int solutionSize = 0;
            counter = 0;
            List<Map<Integer, Integer>> cleanedMCSSeeds = new ArrayList<>();
//            System.out.println("mergin  UIT & KochCliques");
            if (!mcsSeeds.isEmpty()) {
                for (Map<Integer, Integer> map : mcsSeeds) {
                    if (map.size() > solutionSize) {
                        solutionSize = map.size();
                        cleanedMCSSeeds.clear();
                        counter = 0;
                    }
                    if (!map.isEmpty()
                            && map.size() == solutionSize
                            && !super.hasClique(map, cleanedMCSSeeds)) {
                        cleanedMCSSeeds.add(counter, map);
                        counter++;
                    }
                }
            }
            for (Map<Integer, Integer> map : mcsVFSeeds) {
                if (!map.isEmpty()
                        && map.size() >= solutionSize
                        && !super.hasClique(map, cleanedMCSSeeds)) {
                    cleanedMCSSeeds.add(counter, map);
                    counter++;
                }
            }
            /*
             * Sort biggest clique to smallest
             */
            Collections.sort(cleanedMCSSeeds, new Map1ValueComparator(SortOrder.DESCENDING));

            /*
             * Extend the seeds using McGregor
             */
            try {
                super.extendCliquesWithMcGregor(cleanedMCSSeeds);
            } catch (CDKException | IOException ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }

            /*
             * Clear previous seeds
             */
            mcsSeeds.clear();
            cleanedMCSSeeds.clear();

            /*
             * Integerate the solutions
             */
            solutionSize = 0;
            counter = 0;
            this.allAtomMCS = new ArrayList<>();

            /*
             * Store solutions from VF MCS only
             */
            if (!allLocalAtomAtomMapping.isEmpty()) {
                for (AtomAtomMapping atomMCSMap : allLocalAtomAtomMapping) {
                    if (atomMCSMap.getCount() > solutionSize) {
                        solutionSize = atomMCSMap.getCount();
                        allAtomMCS.clear();
                        counter = 0;
                    }
                    if (!atomMCSMap.isEmpty()
                            && atomMCSMap.getCount() == solutionSize) {
                        allAtomMCS.add(counter, atomMCSMap);
                        counter++;
                    }
                }
            }

            /*
             * Clear the local solution after storing it into mcs solutions
             */
            allLocalMCS.clear();
            allLocalAtomAtomMapping.clear();

        } else {

            /*
             * Store solutions from VF MCS only
             */
            int solSize = 0;
            int counter = 0;
            this.allAtomMCS = new ArrayList<>();
            if (!allLocalAtomAtomMapping.isEmpty()) {
                for (AtomAtomMapping atomMCSMap : allLocalAtomAtomMapping) {
                    if (atomMCSMap.getCount() > solSize) {
                        solSize = atomMCSMap.getCount();
                        allAtomMCS.clear();
                        counter = 0;
                    }
                    if (!atomMCSMap.isEmpty()
                            && atomMCSMap.getCount() == solSize) {
                        allAtomMCS.add(counter, atomMCSMap);
                        counter++;
                    }
                }
            }
        }
    }

//    /*
//     * Note: VF MCS will search for cliques which will match the types. Mcgregor will extend the cliques depending of
//     * the bond type (sensitive and insensitive).
//     */
//    protected synchronized boolean searchVFMappings() {
////        System.out.println("searchVFMappings ");
//        IQuery queryCompiler;
//        IMapper mapper;
//
//        if (!(source instanceof IQueryAtomContainer)
//                && !(target instanceof IQueryAtomContainer)) {
//            countR = getReactantMol().getAtomCount();
//            countP = getProductMol().getAtomCount();
//        }
//
//        if (source instanceof IQueryAtomContainer) {
//            queryCompiler = new QueryCompiler((IQueryAtomContainer) source).compile();
//            mapper = new VFMCSMapper(queryCompiler);
//            List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
//            if (maps != null) {
//                vfLibSolutions.addAll(maps);
//            }
//            setVFMappings(true, queryCompiler);
//
//        } else if (countR <= countP) {//isBondMatchFlag()
//            queryCompiler = new QueryCompiler(this.source, true, isMatchRings(), isMatchAtomType()).compile();
//            mapper = new VFMCSMapper(queryCompiler);
//            List<Map<INode, IAtom>> map = mapper.getMaps(this.target);
//            if (map != null) {
//                vfLibSolutions.addAll(map);
//            }
//            setVFMappings(true, queryCompiler);
//        } else {
//            queryCompiler = new QueryCompiler(this.target, true, isMatchRings(), isMatchAtomType()).compile();
//            mapper = new VFMCSMapper(queryCompiler);
//            List<Map<INode, IAtom>> map = mapper.getMaps(this.source);
//            if (map != null) {
//                vfLibSolutions.addAll(map);
//            }
//            setVFMappings(false, queryCompiler);
//        }
//        return mapper.isTimeout();
//    }
    /*
     * Note: CDK VF will search for core hits. Mcgregor will extend the cliques depending of the bond type (sensitive and
     * insensitive).
     */
    private synchronized boolean searchVFCDKMappings() {
        if (DEBUG) {
            System.out.println("searchVFCDKMappings ");
        }

        if (!(source instanceof IQueryAtomContainer)
                && !(target instanceof IQueryAtomContainer)) {

            countR = getReactantMol().getAtomCount();
            countP = getProductMol().getAtomCount();
        }
        if (source instanceof IQueryAtomContainer) {
            VentoFoggia findSubstructure = VentoFoggia.findSubstructure(source, atomMatcher, bondMatcher); // create pattern
            Mappings matchAll = findSubstructure.matchAll((IQueryAtomContainer) target);
            Iterable<Map<IAtom, IAtom>> toAtomMap = matchAll.limit(10).toAtomMap();
            for (Map<IAtom, IAtom> map : toAtomMap) {
                vfLibSolutions.add(map);
            }
            setVFMappings(true);
        } else if (countR <= countP) {

            VentoFoggia findSubstructure = VentoFoggia.findSubstructure(source, atomMatcher, bondMatcher); // create pattern
            Mappings matchAll = findSubstructure.matchAll(target);
            Iterable<Map<IAtom, IAtom>> toAtomMap = matchAll.limit(10).toAtomMap();
            for (Map<IAtom, IAtom> map : toAtomMap) {
                vfLibSolutions.add(map);
            }
            setVFMappings(true);
        } else if (countR > countP) {

            VentoFoggia findSubstructure = VentoFoggia.findSubstructure(target, atomMatcher, bondMatcher); // create pattern
            Mappings matchAll = findSubstructure.matchAll(source);
            Iterable<Map<IAtom, IAtom>> toAtomMap = matchAll.limit(10).toAtomMap();
            for (Map<IAtom, IAtom> map : toAtomMap) {
                vfLibSolutions.add(map);
            }
            setVFMappings(false);
        }

        if (DEBUG) {
            System.out.println("Sol count " + vfLibSolutions.size());
            System.out.println("Sol size " + (vfLibSolutions.iterator().hasNext() ? vfLibSolutions.iterator().next().size() : 0));
        }
        return !vfLibSolutions.isEmpty();
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(getReactantMol(), getProductMol());
    }

    /*
     * Check if fragmented container has single atom
     */
    synchronized boolean isMoleculeConnected(IAtomContainer compound1, IAtomContainer compound2) {

        boolean connected1 = true;

        IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound1);
        for (IAtomContainer a : partitionIntoMolecules.atomContainers()) {

            if (a.getAtomCount() == 1) {
                connected1 = false;
            }
        }

        boolean connected2 = true;

        partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound2);
        for (IAtomContainer a : partitionIntoMolecules.atomContainers()) {

            if (a.getAtomCount() == 1) {
                connected2 = false;
            }
        }
        return connected1 & connected2;
    }
}
