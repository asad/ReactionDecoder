/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.lang.System.currentTimeMillis;
import static java.lang.System.gc;
import static java.lang.System.out;
import static java.util.Collections.unmodifiableMap;
import java.util.EnumMap;
import java.util.Map;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MAX;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIN;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000RXNWriter;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIX;
import java.util.logging.Level;
import uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool;
import static java.util.Collections.synchronizedMap;
import uk.ac.ebi.reactionblast.mapping.container.CDKReactionBuilder;
import static java.lang.String.valueOf;
import static java.util.concurrent.Executors.newCachedThreadPool;
import static java.util.concurrent.Executors.newSingleThreadExecutor;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * @Copyright Syed Asad Rahman (C) 2004-2011
 */
public class CallableAtomMappingTool implements Serializable {

    private final static boolean DEBUG = false;
    private final static ILoggingTool logger
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
        solution = synchronizedMap(new EnumMap<IMappingAlgorithm, Reactor>(IMappingAlgorithm.class));

        logger.info("\n|++++++++++++++++++++++++++++|");
        logger.info(" Standardize Reaction");

        IReaction cleanedReaction;
        try {
            cleanedReaction = standardizer.standardize(reaction);
            cleanedReaction = CDKReactionBuilder.preprocessStandardizedReaction(cleanedReaction);
            generateAtomAtomMapping(cleanedReaction, removeHydrogen);
        } catch (Exception e) {
            logger.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            logger.error(e);
        }
    }

    private synchronized void generateAtomAtomMapping(
            IReaction cleanedReaction, boolean removeHydrogen) {
        ExecutorService executor = null;
        if (DEBUG) {
            executor = newSingleThreadExecutor();
        } else {
            executor = newCachedThreadPool();
        }
        int jobCounter = 0;
        try {
            CompletionService<Reactor> cs = new ExecutorCompletionService<>(executor);

            if (DEBUG) {
                out.println("\nSTEP 2: Calling Mapping Models\n");
            }

            /*
             * MAX Algorithm
             */
            logger.info("\n|++++++++++++++++++++++++++++|");
            logger.info("a) Global Model: ");

            MappingThread maxThread = new MappingThread("IMappingAlgorithm.MAX", ExtReactionManipulatorTool.deepClone(cleanedReaction), MAX, removeHydrogen);
            cs.submit(maxThread);
            jobCounter++;

            /*
             * MIN Algorithm
             */
            logger.info("\n|++++++++++++++++++++++++++++|");
            logger.info("c) Local Model: ");

            MappingThread minThread = new MappingThread("IMappingAlgorithm.MIN", ExtReactionManipulatorTool.deepClone(cleanedReaction), MIN, removeHydrogen);
            cs.submit(minThread);
            jobCounter++;

            /*
             * MIX Algorithm
             */
            logger.info("\n|++++++++++++++++++++++++++++|");
            logger.info("b) Mixture Model: ");

            MappingThread maxMixtureThread = new MappingThread("IMappingAlgorithm.MIX", ExtReactionManipulatorTool.deepClone(cleanedReaction), MIX, removeHydrogen);
            cs.submit(maxMixtureThread);
            jobCounter++;

            /*
             * RINGS Minimization
             */
            logger.info("\n|++++++++++++++++++++++++++++|");
            logger.info("d) Rings Model: ");

            MappingThread ringThread = new MappingThread("IMappingAlgorithm.RINGS", ExtReactionManipulatorTool.deepClone(cleanedReaction), RINGS, removeHydrogen);
            cs.submit(ringThread);
            jobCounter++;

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
            logger.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            logger.error(e);
        } catch (CloneNotSupportedException ex) {
            Logger.getLogger(CallableAtomMappingTool.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            executor.shutdown();
        }

        logger.info("!!!!Atom-Atom Mapping Done!!!!");
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
