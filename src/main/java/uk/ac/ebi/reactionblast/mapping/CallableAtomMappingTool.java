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
package uk.ac.ebi.reactionblast.mapping;

import java.io.File;
import static java.io.File.separator;
import java.io.FileWriter;
import java.io.Serializable;
import static java.lang.Runtime.getRuntime;
import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.gc;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import static java.util.Collections.synchronizedMap;
import static java.util.Collections.unmodifiableMap;
import java.util.EnumMap;
import java.util.Map;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MAX;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIN;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIXTURE;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000RXNWriter;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * @Copyright Syed Asad Rahman (C) 2004-2018
 */
public class CallableAtomMappingTool implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static boolean DEBUG = false;
    private final static ILoggingTool LOGGER
            = createLoggingTool(CallableAtomMappingTool.class);
    private static final long serialVersionUID = 0x29e2adb1716b13eL;

    /**
     * Creates mapping PDFs for all the processed reaction mappings
     *
     * @param reactor
     * @param outputDirectoryName
     * @param outFileName
     * @throws Exception
     */
    public static synchronized void writeMappingRXN(Reactor reactor, String outputDirectoryName, String outFileName) throws Exception {
        String reactionID = reactor.getReactionWithAtomAtomMapping().getID();
        IReaction mappedReaction = reactor.getReactionWithAtomAtomMapping();
        if (reactionID == null) {
            reactionID = valueOf(currentTimeMillis());
            reactor.getReactionWithAtomAtomMapping().setID(reactionID);
        }

        String outputFile = outputDirectoryName;
        outputFile += separator + outFileName;
        try (MDLV2000RXNWriter rxnW = new MDLV2000RXNWriter(new FileWriter(new File(outputFile)))) {
            rxnW.write(mappedReaction);
        }
    }
    private Map<IMappingAlgorithm, Reactor> solution = null;

    /**
     * Takes a standardizer to standardize the reaction before mapping.
     *
     * @param reaction
     * @param standardizer
     * @param removeHydrogen
     * @throws Exception
     */
    public CallableAtomMappingTool(
            IReaction reaction,
            IStandardizer standardizer,
            boolean removeHydrogen) throws Exception {
        solution = synchronizedMap(new EnumMap<>(IMappingAlgorithm.class));
        generateAtomAtomMapping(reaction, standardizer, removeHydrogen);
    }

    private synchronized void generateAtomAtomMapping(
            IReaction reaction,
            IStandardizer standardizer,
            boolean removeHydrogen) {
        ExecutorService executor = null;
        if (DEBUG) {
            executor = Executors.newSingleThreadExecutor();
        } else {
            int threadsAvailable = getRuntime().availableProcessors() - 1;
            if (threadsAvailable == 0) {
                threadsAvailable = 1;
            }

            if (threadsAvailable > 4) {
                threadsAvailable = 4;
            }
            executor = Executors.newFixedThreadPool(threadsAvailable);
        }
        int jobCounter = 0;
        try {
            CompletionService<Reactor> cs = new ExecutorCompletionService<>(executor);

            /*
             * MAX Algorithm
             */
            LOGGER.info(NEW_LINE + "|++++++++++++++++++++++++++++|");
            LOGGER.info("a) Global Model: ");
            if (DEBUG) {
                out.println(NEW_LINE + "-----------------------------------" + NEW_LINE);
                out.println(NEW_LINE + "STEP 1: Global Model Standardize Reactions" + NEW_LINE);
            }
            IReaction cleanedReaction1 = null;
            try {
                cleanedReaction1 = standardizer.standardize(reaction);
            } catch (Exception e) {
                LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
                LOGGER.error(e);
            }
            if (DEBUG) {
                out.println(NEW_LINE + "STEP a: Calling Mapping Models" + NEW_LINE);
            }
            MappingThread maxThread = new MappingThread("IMappingAlgorithm.MAX", cleanedReaction1, MAX, removeHydrogen);
            cs.submit(maxThread);
            jobCounter++;
            /*
             * MIN Algorithm
             */
            LOGGER.info(NEW_LINE + "|++++++++++++++++++++++++++++|");
            LOGGER.info("b) Local Model: ");
            if (DEBUG) {
                out.println(NEW_LINE + "-----------------------------------" + NEW_LINE);
                out.println(NEW_LINE + "STEP b: Local Model Standardize Reactions" + NEW_LINE);
            }
            IReaction cleanedReaction2 = null;
            try {
                cleanedReaction2 = standardizer.standardize(reaction);
            } catch (Exception e) {
                LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
                LOGGER.error(e);
            }
            MappingThread minThread = new MappingThread("IMappingAlgorithm.MIN", cleanedReaction2, MIN, removeHydrogen);
            cs.submit(minThread);
            jobCounter++;
            /*
             * MIXTURE Algorithm
             */
            LOGGER.info(NEW_LINE + "|++++++++++++++++++++++++++++|");
            LOGGER.info("c) Mixture Model: ");
            if (DEBUG) {
                out.println(NEW_LINE + "-----------------------------------" + NEW_LINE);
                out.println(NEW_LINE + "STEP c: Mixture Model Standardize Reactions" + NEW_LINE);
            }
            IReaction cleanedReaction3 = null;
            try {
                cleanedReaction3 = standardizer.standardize(reaction);
            } catch (Exception e) {
                LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
                LOGGER.error(e);
            }
            MappingThread maxMixtureThread = new MappingThread("IMappingAlgorithm.MIXTURE", cleanedReaction3, MIXTURE, removeHydrogen);
            cs.submit(maxMixtureThread);
            jobCounter++;
//
//            /*
//             * RINGS Minimization
//             */
//            LOGGER.info(NEW_LINE + "|++++++++++++++++++++++++++++|");
//            LOGGER.info("d) Rings Model: ");
//            if (DEBUG) {
//                out.println(NEW_LINE + "-----------------------------------" + NEW_LINE);
//                out.println(NEW_LINE + "STEP d: Rings Model Standardize Reactions" + NEW_LINE);
//            }
//            IReaction cleanedReaction4 = null;
//            try {
//                cleanedReaction4 = standardizer.standardize(reaction);
//            } catch (Exception e) {
//                LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
//                LOGGER.error(e);
//            }
//            MappingThread ringThread = new MappingThread("IMappingAlgorithm.RINGS", cleanedReaction4, RINGS, removeHydrogen);
//            cs.submit(ringThread);
//            jobCounter++;

            /*
             * Collect the results
             */
            for (int i = 0; i < jobCounter; i++) {
                Reactor chosen = cs.take().get();
                putSolution(chosen.getAlgorithm(), chosen);
            }
            executor.shutdown();
            /*
             Wait until all threads are finish
             *
             */
            while (!executor.isTerminated()) {
            }
            gc();
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            LOGGER.error(e);
        } finally {
            executor.shutdown();
        }

        LOGGER.info("!!!!Atom-Atom Mapping Done!!!!");
    }

    /**
     * @return the solution
     */
    public Map<IMappingAlgorithm, Reactor> getSolutions() {
        return unmodifiableMap(solution);
    }

    /**
     * @param solution the solution to set
     */
    private void putSolution(IMappingAlgorithm choice, Reactor reactor) {
        this.solution.put(choice, reactor);
    }

}
