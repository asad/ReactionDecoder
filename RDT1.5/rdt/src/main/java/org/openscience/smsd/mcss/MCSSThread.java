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

import static java.lang.Integer.MAX_VALUE;
import static java.util.Calendar.getInstance;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import static org.openscience.smsd.interfaces.Algorithm.DEFAULT;
import static org.openscience.smsd.mcss.JobType.MULTIPLE;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.removeHydrogens;

/**
 * 
 * 
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public class MCSSThread implements Callable<LinkedBlockingQueue<IAtomContainer>> {

    private final static ILoggingTool logger
            = createLoggingTool(MCSSThread.class);
    private static final Logger LOG = getLogger(MCSSThread.class.getName());
    private final List<IAtomContainer> mcssList;
    private final JobType jobType;
    private final int taskNumber;
    private final boolean matchBonds;
    private final boolean matchRings;
    private final boolean matchAtomType;

    /**
     *
     * @param mcssList
     * @param jobType MULTIPLE/SINGLE
     * @param taskNumber
     */
    public MCSSThread(List<IAtomContainer> mcssList, JobType jobType, int taskNumber) {
        this(mcssList, jobType, taskNumber, true, true, true);
    }

    /**
     *
     * @param mcssList
     * @param jobType
     * @param taskNumber
     * @param matchBonds
     * @param matchRings
     * @param matchAtomType
     */
    public MCSSThread(List<IAtomContainer> mcssList, JobType jobType, int taskNumber, boolean matchBonds, boolean matchRings,
            boolean matchAtomType) {
        this.mcssList = mcssList;
        this.jobType = jobType;
        this.taskNumber = taskNumber;
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        this.matchAtomType = matchAtomType;

    }

    @Override
    public synchronized LinkedBlockingQueue<IAtomContainer> call() {
        if (this.jobType.equals(MULTIPLE)) {
            return multiSolution();
        } else {
            return singleSolution();
        }
    }
    /*
     * MULTIPLE Fragments of MCS are returned if present
     */

    private synchronized LinkedBlockingQueue<IAtomContainer> multiSolution() {
        /*
         * Store final solution here
         */
        LinkedBlockingQueue<IAtomContainer> mcss = new LinkedBlockingQueue<>();

        logger.debug("Calling MCSSTask " + taskNumber + " with " + mcssList.size() + " items");
        long startTime = getInstance().getTimeInMillis();
        IAtomContainer querySeed = mcssList.get(0);
        long calcTime = startTime;

        ConcurrentLinkedQueue<IAtomContainer> seeds = new ConcurrentLinkedQueue<>();
        try {
            /*
             * Local Seeds
             */
            Set<Fragment> localSeeds = new TreeSet<>();
            int minSeedSize = querySeed.getAtomCount();

            for (int index = 1; index < mcssList.size(); index++) {
                IAtomContainer target = mcssList.get(index);
                Collection<Fragment> fragmentsFromMCS;
                BaseMapping comparison;

                comparison = new Isomorphism(querySeed, target, DEFAULT, matchBonds, matchRings, matchAtomType);
                comparison.setChemFilters(true, true, true);
                fragmentsFromMCS = getMCSS(comparison);

                logger.debug("comparison for task " + taskNumber + " has " + fragmentsFromMCS.size()
                        + " unique matches of size " + comparison.getFirstAtomMapping().getCount());
                logger.debug("MCSS for task " + taskNumber + " has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
                logger.debug("Target for task " + taskNumber + " has " + target.getAtomCount() + " atoms, and " + target.getBondCount() + " bonds");
                long endCalcTime = getInstance().getTimeInMillis();
                logger.debug("Task " + taskNumber + " index " + index + " took " + (endCalcTime - calcTime) + "ms");
                calcTime = endCalcTime;

                if (fragmentsFromMCS.isEmpty()) {
                    localSeeds.clear();
                    break;
                }
                Iterator<Fragment> iterator = fragmentsFromMCS.iterator();
                /*
                 * Store rest of the unique hits
                 */
                while (iterator.hasNext()) {
                    Fragment fragment = iterator.next();
                    if (minSeedSize > fragment.getContainer().getAtomCount()) {
                        localSeeds.clear();
                        minSeedSize = fragment.getContainer().getAtomCount();
                    }
                    if (minSeedSize == fragment.getContainer().getAtomCount()) {
                        localSeeds.add(fragment);
                    }
                }
            }
            /*
             * Add all the Maximum Unique Substructures
             */
            if (!localSeeds.isEmpty()) {
                for (Fragment f : localSeeds) {
                    seeds.add(f.getContainer());
                }
                localSeeds.clear();
            }

            logger.debug("No of Potential MULTIPLE " + seeds.size());

            /*
             * Choose only cleaned MULTIPLE Substructures
             */
            minSeedSize = MAX_VALUE;

            while (!seeds.isEmpty()) {
                IAtomContainer fragmentMCS = seeds.poll();
                localSeeds = new TreeSet<>();
                logger.debug("Potential MULTIPLE " + getMCSSSmiles(fragmentMCS));
                Collection<Fragment> fragmentsFromMCS;
                for (IAtomContainer target : mcssList) {
                    Isomorphism comparison = new Isomorphism(fragmentMCS, target, DEFAULT, matchBonds, matchRings, matchAtomType);
                    comparison.setChemFilters(true, true, true);
                    fragmentsFromMCS = getMCSS(comparison);

                    /*
                     * Only true MCSS is added
                     */
                    if (fragmentsFromMCS == null || fragmentsFromMCS.isEmpty()) {
                        localSeeds.clear();
                        break;
                    }
                    Iterator<Fragment> iterator = fragmentsFromMCS.iterator();
                    /*
                     * Store rest of the unique hits
                     */
                    while (iterator.hasNext()) {
                        Fragment fragment = iterator.next();
                        if (minSeedSize > fragment.getContainer().getAtomCount()) {
                            localSeeds.clear();
                            minSeedSize = fragment.getContainer().getAtomCount();
                        }
                        if (minSeedSize == fragment.getContainer().getAtomCount()) {
                            localSeeds.add(fragment);
                        }
                    }
                    /*
                     * Top solution
                     */
                    fragmentMCS = localSeeds.iterator().next().getContainer();
                }

                /*
                 * Add all the Maximum Unique Substructures
                 */
                if (!localSeeds.isEmpty()) {
                    for (Fragment f : localSeeds) {
                        mcss.add(f.getContainer());
                    }
                    localSeeds.clear();
                }

            }
        } catch (CDKException e) {
            logger.error("ERROR IN MCS Thread: ", e);
        }
        long endTime = getInstance().getTimeInMillis();
        logger.debug("Done: task " + taskNumber + " took " + (endTime - startTime) + "ms");
        logger.debug(" and mcss has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
        return mcss;
    }
    /*
     * SINGLE Fragment of MCS is returned if present.
     */

    private synchronized LinkedBlockingQueue<IAtomContainer> singleSolution() {

        logger.debug("Calling MCSSTask " + taskNumber + " with " + mcssList.size() + " items");
        LinkedBlockingQueue<IAtomContainer> mcss = new LinkedBlockingQueue<>();
        long startTime = getInstance().getTimeInMillis();
        IAtomContainer querySeed = mcssList.get(0);
        long calcTime = startTime;

        try {
            for (int index = 1; index < mcssList.size(); index++) {
                IAtomContainer target = removeHydrogens(mcssList.get(index));
                Collection<Fragment> fragmentsFomMCS;
                BaseMapping comparison;

                comparison = new Isomorphism(querySeed, target, DEFAULT, matchBonds, matchRings, matchAtomType);
                comparison.setChemFilters(true, true, true);
                fragmentsFomMCS = getMCSS(comparison);

                logger.debug("comparison for task " + taskNumber + " has " + fragmentsFomMCS.size()
                        + " unique matches of size " + comparison.getFirstAtomMapping().getCount());
                logger.debug("MCSS for task " + taskNumber + " has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
                logger.debug("Target for task " + taskNumber + " has " + target.getAtomCount() + " atoms, and " + target.getBondCount() + " bonds");
                long endCalcTime = getInstance().getTimeInMillis();
                logger.debug("Task " + taskNumber + " index " + index + " took " + (endCalcTime - calcTime) + "ms");
                calcTime = endCalcTime;

                if (fragmentsFomMCS.isEmpty()) {
                    break;
                }
                querySeed = fragmentsFomMCS.iterator().next().getContainer();
            }

            if (querySeed != null) {
                mcss.add(querySeed);
                long endTime = getInstance().getTimeInMillis();
                logger.debug("Done: task " + taskNumber + " took " + (endTime - startTime) + "ms");
                logger.debug(" and mcss has " + querySeed.getAtomCount() + " atoms, and " + querySeed.getBondCount() + " bonds");
            }
        } catch (Exception e) {
            logger.error("ERROR IN MCS Thread: ", e);
        }
        return mcss;
    }

    private synchronized Collection<Fragment> getMCSS(BaseMapping comparison) {
        Set<Fragment> matchList = new HashSet<>();
        for (AtomAtomMapping mapping : comparison.getAllAtomMapping()) {
            IAtomContainer match;
            try {
                match = mapping.getCommonFragment();
                try {
                    matchList.add(new Fragment(match));
                } catch (CDKException ex) {
                    logger.error("ERROR IN MCS Thread: ", ex);
                }
            } catch (CloneNotSupportedException ex) {
                logger.error("ERROR IN MCS Thread: ", ex);
            }
        }
        return matchList;
    }

    /**
     * Return SMILES
     *
     * @param ac
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    public synchronized String getMCSSSmiles(IAtomContainer ac) throws CDKException {
        SmilesGenerator g = new SmilesGenerator().aromatic();
        return g.create(ac);
    }

    /**
     * @return the taskNumber
     */
    public synchronized int getTaskNumber() {
        return taskNumber;
    }
}
