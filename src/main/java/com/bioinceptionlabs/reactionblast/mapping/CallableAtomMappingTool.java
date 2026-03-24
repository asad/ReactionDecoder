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
package com.bioinceptionlabs.reactionblast.mapping;

import java.io.File;
import static java.io.File.separator;
import java.io.FileWriter;
import java.io.Serializable;
import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.getProperty;
import static java.util.Collections.unmodifiableMap;
import java.util.EnumMap;
import java.util.Map;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import com.bioinceptionlabs.reactionblast.interfaces.IStandardizer;
import com.bioinceptionlabs.reactionblast.mapping.cache.ThreadSafeCache;
import com.bioinceptionlabs.reactionblast.mapping.graph.MCSSolution;
import com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm.MAX;
import static com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm.MIN;
import static com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm.MIXTURE;
import static com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;
import com.bioinceptionlabs.reactionblast.tools.rxnfile.MDLV2000RXNWriter;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 * @Copyright Syed Asad Rahman (C) 2004-2020
 */
public class CallableAtomMappingTool implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
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
    public static void writeMappingRXN(Reactor reactor, String outputDirectoryName, String outFileName) throws Exception {
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
     * @param checkComplex will check complex mapping like ring systems
     * @throws Exception
     */
    public CallableAtomMappingTool(
            IReaction reaction,
            IStandardizer standardizer,
            boolean removeHydrogen,
            boolean checkComplex) throws Exception {
        solution = new EnumMap<>(IMappingAlgorithm.class);
        generateAtomAtomMapping(reaction, standardizer, removeHydrogen, checkComplex);
    }

    private void generateAtomAtomMapping(
            IReaction reaction,
            IStandardizer standardizer,
            boolean removeHydrogen,
            boolean checkComplex) {
        /*
         * Mapping cache initialized
         */
        ThreadSafeCache<String, MCSSolution> mappingcache = ThreadSafeCache.getInstance();

        /*
         * Standardize the reaction ONCE and clone for each algorithm.
         * Previously this was done 4 times independently — major overhead.
         */
        IReaction standardizedReaction = null;
        try {
            standardizedReaction = standardizer.standardize(reaction);
        } catch (Exception e) {
            LOGGER.debug("ERROR: in AtomMappingTool standardization: " + e.getMessage());
            LOGGER.error(e);
        }
        if (standardizedReaction == null) {
            LOGGER.error("Failed to standardize reaction — cannot proceed with mapping");
            return;
        }

        int numAlgorithms = checkComplex ? 4 : 3;
        ExecutorService executor = Executors.newFixedThreadPool(numAlgorithms);
        int jobCounter = 0;
        try {
            CompletionService<Reactor> cs = new ExecutorCompletionService<>(executor);

            /*
             * MIN Algorithm
             */
            LOGGER.debug("Submitting MIN algorithm");
            IReaction cleanedReaction1 = cloneReaction(standardizedReaction);
            cs.submit(new MappingThread("IMappingAlgorithm.MIN", cleanedReaction1, MIN, removeHydrogen));
            jobCounter++;

            /*
             * MAX Algorithm
             */
            LOGGER.debug("Submitting MAX algorithm");
            IReaction cleanedReaction2 = cloneReaction(standardizedReaction);
            cs.submit(new MappingThread("IMappingAlgorithm.MAX", cleanedReaction2, MAX, removeHydrogen));
            jobCounter++;

            /*
             * MIXTURE Algorithm
             */
            LOGGER.debug("Submitting MIXTURE algorithm");
            IReaction cleanedReaction3 = cloneReaction(standardizedReaction);
            cs.submit(new MappingThread("IMappingAlgorithm.MIXTURE", cleanedReaction3, MIXTURE, removeHydrogen));
            jobCounter++;

            if (checkComplex) {
                /*
                 * RINGS Algorithm
                 */
                LOGGER.debug("Submitting RINGS algorithm");
                IReaction cleanedReaction4 = cloneReaction(standardizedReaction);
                cs.submit(new MappingThread("IMappingAlgorithm.RINGS", cleanedReaction4, RINGS, removeHydrogen));
                jobCounter++;
            }

            /*
             * Collect the results — all algorithms run in parallel
             */
            for (int i = 0; i < jobCounter; i++) {
                Reactor chosen = cs.take().get();
                putSolution(chosen.getAlgorithm(), chosen);
            }
            executor.shutdown();
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            LOGGER.debug("======DONE CallableAtomMappingTool=======");
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            LOGGER.error(e);
        } finally {
            executor.shutdown();
        }
        LOGGER.debug("!!!!Atom-Atom Mapping Done!!!!");
        /*
         * Mapping cache cleared
         */
        mappingcache.cleanup();
    }

    /**
     * Deep-clone a reaction so each algorithm gets an independent copy.
     */
    private IReaction cloneReaction(IReaction reaction) {
        try {
            return (IReaction) reaction.clone();
        } catch (CloneNotSupportedException e) {
            LOGGER.error("Failed to clone reaction: " + e.getMessage());
            return reaction;
        }
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
