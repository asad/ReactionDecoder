/*
 * Copyright (C) 2014 Syed Asad Rahman <asad at ebi.ac.uk>.
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
package org.openscience.smsd.mcss;

import static java.lang.Math.ceil;
import static java.lang.Runtime.getRuntime;
import static java.lang.System.gc;
import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.sort;
import static java.util.Collections.unmodifiableCollection;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import static java.util.concurrent.Executors.newFixedThreadPool;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.tools.AtomContainerComparator;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.removeHydrogens;

/**
 *
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public class MCSS {

    static final String NEW_LINE = getProperty("line.separator");
    private final static ILoggingTool LOGGER = createLoggingTool(MCSS.class);
    private final Collection<IAtomContainer> calculateMCSS;
    private final boolean matchBonds;
    private final boolean matchRings;
    private final boolean matchAtomType;

    /**
     *
     * @param jobList
     * @param jobType
     * @param numberOfThreads
     */
    public MCSS(List<IAtomContainer> jobList, JobType jobType, int numberOfThreads) {
        this(jobList, jobType, numberOfThreads, true, true, true);
    }

    /**
     *
     * @param jobList
     * @param jobType
     * @param numberOfThreads
     * @param matchBonds
     * @param matchRings
     * @param matchAtomType
     */
    public MCSS(
            List<IAtomContainer> jobList,
            JobType jobType,
            int numberOfThreads,
            boolean matchBonds,
            boolean matchRings,
            boolean matchAtomType) {
        int threadsAvailable = getRuntime().availableProcessors() - 1;

        LOGGER.debug("Demand threads: " + numberOfThreads);
        LOGGER.debug(", Available threads: " + threadsAvailable);
        if (numberOfThreads > 0 && threadsAvailable >= numberOfThreads) {
            threadsAvailable = numberOfThreads;
        } else if (threadsAvailable <= 0) {
            threadsAvailable = 1;
        }
        LOGGER.debug(", Assigned threads: " + threadsAvailable + NEW_LINE);
        /*
         * Remove hydrogen from the molecules
         **/
        List<IAtomContainer> selectedJobs = new ArrayList<>(jobList.size());
        jobList.stream().forEach((ac) -> {
            selectedJobs.add(removeHydrogens(ac));
        });
        /*
         * Sort the molecules in the ascending order of atom size and atom type
         */
        Comparator<IAtomContainer> comparator = new AtomContainerComparator();
        sort(selectedJobs, comparator);
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        this.matchAtomType = matchAtomType;
        /*
         * Call the MCS
         */
        calculateMCSS = calculateMCSS(selectedJobs, jobType, threadsAvailable);
        selectedJobs.clear();
    }

    private synchronized Collection<IAtomContainer> calculateMCSS(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        List<IAtomContainer> newMCSSList;
        if (nThreads == 1) {
            newMCSSList = new LinkedList<>(submitSingleThreadedJob(mcssList, jobType, nThreads));
        } else {
            /*
             * Calling recursive MCS
             */
            newMCSSList = new LinkedList<>(submitMultiThreadedJob(mcssList, jobType, nThreads));
            while (newMCSSList.size() > 1) {
                if (newMCSSList.size() > 2) {
                    newMCSSList = new LinkedList<>(submitMultiThreadedJob(newMCSSList, jobType, nThreads));
                } else {
                    newMCSSList = new LinkedList<>(submitMultiThreadedJob(newMCSSList, jobType, 1));
                }
            }
        }
        if (!mcssList.isEmpty() && !newMCSSList.isEmpty()) {
            IAtomContainer inTheList = mcssList.get(mcssList.size() - 1);
            if (inTheList == newMCSSList.iterator().next()) {
                return new LinkedBlockingQueue<>();
            }
        }
        return newMCSSList;
    }

    /**
     * @return the calculateMCSS
     */
    public synchronized Collection<IAtomContainer> getCalculateMCSS() {
        return unmodifiableCollection(calculateMCSS);
    }

    private synchronized LinkedBlockingQueue<IAtomContainer> submitSingleThreadedJob(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        LinkedBlockingQueue<IAtomContainer> solutions = new LinkedBlockingQueue<>();
        MCSSThread task = new MCSSThread(mcssList, jobType, 1);
        LinkedBlockingQueue<IAtomContainer> results = task.call();
        if (results != null) {
            solutions.addAll(results);
        }
        return solutions;
    }

    private synchronized LinkedBlockingQueue<IAtomContainer> submitMultiThreadedJob(List<IAtomContainer> mcssList, JobType jobType, int nThreads) {
        int taskNumber = 1;
        LinkedBlockingQueue<IAtomContainer> solutions = new LinkedBlockingQueue<>();
        LinkedBlockingQueue<Callable<LinkedBlockingQueue<IAtomContainer>>> callablesQueue = new LinkedBlockingQueue<>();
        ExecutorService threadPool = newFixedThreadPool(nThreads);
        int step = (int) ceil(mcssList.size() / nThreads);
        if (step < 2) {
            step = 2; // Can't have a step size of less than 2
        }
        for (int i = 0; i < mcssList.size(); i += step) {
            int endPoint = i + step;
            if (endPoint > mcssList.size()) {
                endPoint = mcssList.size();
            }
            List<IAtomContainer> subList = new ArrayList<>(mcssList.subList(i, endPoint));
            if (subList.size() > 1) {
                MCSSThread mcssJobThread = new MCSSThread(subList, jobType, taskNumber, matchBonds, matchRings,
                        matchAtomType);
                callablesQueue.add(mcssJobThread);
                taskNumber++;
            } else {
                solutions.add(subList.get(0));
            }
        }
        try {
            /*
             * Wait all the threads to finish
             */
            List<Future<LinkedBlockingQueue<IAtomContainer>>> futureList = threadPool.invokeAll(callablesQueue);
            /*
             * Collect the results
             */
            for (Future<LinkedBlockingQueue<IAtomContainer>> callable : futureList) {
                LinkedBlockingQueue<IAtomContainer> mapping = callable.get();
                if (callable.isDone() && mapping != null) {
                    solutions.addAll(mapping);
                } else {
                    LOGGER.warn("WARNING: InComplete job in AtomMappingTool: ");
                }
            }
            threadPool.shutdown();
            // Wait until all threads are finish
            while (!threadPool.isTerminated()) {
            }
            gc();
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            LOGGER.error(e);
        } finally {
            threadPool.shutdown();
        }

        return solutions;
    }

    /**
     *
     * @return
     */
    public synchronized String getTitle() {
        return "Calculating Maximum Commmon Substrutures (MCSS) using SMSD";
    }
}
