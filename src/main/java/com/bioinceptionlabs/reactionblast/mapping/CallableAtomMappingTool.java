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

import java.io.Serializable;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.getProperty;
import static java.util.Collections.unmodifiableMap;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import com.bioinceptionlabs.reactionblast.mapping.ThreadSafeCache;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MAX;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MIN;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MIXTURE;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.RINGS;

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
    private static final int MAPPING_PARALLELISM
            = Math.max(2, Math.min(3, Runtime.getRuntime().availableProcessors()));
    private static final ExecutorService MAPPING_EXECUTOR
            = Executors.newFixedThreadPool(MAPPING_PARALLELISM, new MappingThreadFactory());

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
            StandardizeReaction standardizer,
            boolean removeHydrogen,
            boolean checkComplex) throws Exception {
        solution = new EnumMap<>(IMappingAlgorithm.class);
        generateAtomAtomMapping(reaction, standardizer, removeHydrogen, checkComplex);
    }

    /**
     * Funnel architecture: pick the cheapest high-value first algorithm,
     * check quality, only run remaining algorithms if that first pass is
     * insufficient. RINGS leads only for ring-containing reactions; otherwise
     * MIN leads on acyclic reactions.
     *
     * Quality gate: if the first pass produces a near-complete mapping, accept
     * it immediately and skip the rest of the algorithm family.
     */
    private void generateAtomAtomMapping(
            IReaction reaction,
            StandardizeReaction standardizer,
            boolean removeHydrogen,
            boolean checkComplex) {
        long mappingStart = currentTimeMillis();
        /*
         * Standardize the reaction ONCE.
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

        MappingDiagnostics.resetReaction(standardizedReaction.getID());

        if (isIdentityReaction(standardizedReaction)) {
            try {
                Reactor minResult = new MappingThread(
                        "IMappingAlgorithm.MIN", standardizedReaction, MIN, removeHydrogen).call();
                putSolution(MIN, minResult);
            } catch (InterruptedException | ExecutionException e) {
                LOGGER.debug("MIN identity phase failed: " + e.getMessage());
                LOGGER.error(e);
            } catch (Exception e) {
                LOGGER.debug("MIN identity phase failed: " + e.getMessage());
                LOGGER.error(e);
            } finally {
                ThreadSafeCache.getInstance().cleanup();
            }
            return;
        }

        /*
         * Phase 1: Run RINGS first if checkComplex is true (most common case).
         * RINGS handles ring-containing molecules best and covers ~75% of
         * drug-like / organic reactions.
         *
         * Skip funnel for large multi-substrate reactions (>5 molecules)
         * where RINGS alone is unlikely to succeed.
         */
        int totalMolecules = standardizedReaction.getReactantCount()
                + standardizedReaction.getProductCount();
        boolean hasRings = hasRingSystems(standardizedReaction);
        IMappingAlgorithm firstPass = (checkComplex && hasRings) ? RINGS : MIN;
        if (totalMolecules <= 5) {
            try {
                Reactor firstPassResult = new MappingThread(
                        "IMappingAlgorithm." + firstPass.name(),
                        standardizedReaction, firstPass, removeHydrogen).call();
                putSolution(firstPass, firstPassResult);

                if (isMappingAcceptable(firstPassResult)) {
                    LOGGER.debug(firstPass + " mapping accepted — skipping remaining algorithms");
                    ThreadSafeCache.getInstance().cleanup();
                    return;
                }
                LOGGER.debug(firstPass + " mapping insufficient — running remaining algorithms");
            } catch (InterruptedException | ExecutionException e) {
                LOGGER.debug(firstPass + " phase failed: " + e.getMessage());
                LOGGER.error(e);
            } catch (Exception e) {
                LOGGER.debug(firstPass + " phase failed: " + e.getMessage());
                LOGGER.error(e);
            }
        }

        /*
         * Phase 2: Run remaining algorithms in parallel (only if RINGS wasn't enough).
         * If funnel was skipped (large reaction), run all 4 algorithms.
         */
        boolean minAlreadyRun = solution.containsKey(MIN);
        boolean ringsAlreadyRun = solution.containsKey(RINGS);
        IMappingAlgorithm[] remaining;
        if (minAlreadyRun && ringsAlreadyRun) {
            remaining = new IMappingAlgorithm[]{MAX};
        } else if (minAlreadyRun) {
            remaining = new IMappingAlgorithm[]{MAX, RINGS};
        } else if (ringsAlreadyRun) {
            remaining = new IMappingAlgorithm[]{MIN, MAX};
        } else {
            remaining = new IMappingAlgorithm[]{MIN, MAX, RINGS};
        }

        try {
            CompletionService<Reactor> cs = new ExecutorCompletionService<>(MAPPING_EXECUTOR);
            int jobCounter = 0;
            for (IMappingAlgorithm algo : remaining) {
                LOGGER.debug("Submitting " + algo.description());
                cs.submit(new MappingThread("IMappingAlgorithm." + algo.name(),
                        standardizedReaction, algo, removeHydrogen));
                jobCounter++;
            }
            for (int i = 0; i < jobCounter; i++) {
                Reactor chosen = cs.take().get();
                putSolution(chosen.getAlgorithm(), chosen);
            }
            LOGGER.debug("======DONE CallableAtomMappingTool=======");
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
            LOGGER.error(e);
        } finally {
            if (standardizedReaction != null && standardizedReaction.getID() != null) {
                MappingDiagnostics.recordMappingPhase(
                        standardizedReaction.getID(),
                        currentTimeMillis() - mappingStart);
            }
            LOGGER.debug("!!!!Atom-Atom Mapping Done!!!!");
            ThreadSafeCache.getInstance().cleanup();
        }
    }

    /**
     * Quality gate for funnel architecture.
     * Checks if a mapping result is "good enough" to skip remaining algorithms.
     *
     * Criteria:
     * 1. Reactor must not be null and must have a valid mapped reaction
     * 2. All non-hydrogen reactant atoms must be mapped
     * 3. Reaction must NOT be an identity/transporter (reactants ≡ products)
     *    — identity reactions need MIN algorithm for correct zero-change detection
     * 4. The reaction must have distinct reactants and products (not a no-op)
     *
     * This is a conservative gate — it accepts the RINGS result only when
     * the mapping is complete and the reaction involves actual bond changes.
     */
    private boolean isMappingAcceptable(Reactor reactor) {
        if (reactor == null) {
            return false;
        }
        try {
            IReaction mapped = reactor.getReactionWithAtomAtomMapping();
            if (mapped == null) {
                return false;
            }

            if (!isAtomBalanced(mapped)) {
                LOGGER.debug("Unbalanced reaction detected — need full pipeline");
                return false;
            }

            // Check if this is an identity/transporter reaction (reactants ≡ products).
            // These need the full pipeline because MIN correctly detects zero change.
            if (isIdentityReaction(mapped)) {
                LOGGER.debug("Identity/transporter reaction detected — need full pipeline");
                return false;
            }

            // Check that all non-H atoms in reactants have been mapped
            int totalReactantAtoms = 0;
            int mappedReactantAtoms = 0;
            for (IAtomContainer ac : mapped.getReactants().atomContainers()) {
                for (IAtom atom : ac.atoms()) {
                    if (!"H".equals(atom.getSymbol())) {
                        totalReactantAtoms++;
                        if (atom.getProperty(org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING) != null) {
                            Object mapNum = atom.getProperty(org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING);
                            if (mapNum instanceof Integer && (Integer) mapNum > 0) {
                                mappedReactantAtoms++;
                            }
                        }
                    }
                }
            }

            if (totalReactantAtoms == 0) {
                return false;
            }

            double mappingCoverage = (double) mappedReactantAtoms / totalReactantAtoms;
            LOGGER.debug("RINGS mapping coverage: " + mappedReactantAtoms + "/" + totalReactantAtoms
                    + " (" + String.format("%.1f%%", mappingCoverage * 100) + ")");

            // Accept if ≥95% of atoms are mapped (allowing small gaps for reagent atoms)
            return mappingCoverage >= 0.95;

        } catch (Exception e) {
            LOGGER.debug("Error checking mapping quality: " + e.getMessage());
            return false;
        }
    }

    private boolean isAtomBalanced(IReaction reaction) {
        return countAtoms(reaction.getReactants()).equals(countAtoms(reaction.getProducts()));
    }

    private Map<String, Integer> countAtoms(org.openscience.cdk.interfaces.IAtomContainerSet molSet) {
        Map<String, Integer> counts = new HashMap<>();
        for (IAtomContainer container : molSet.atomContainers()) {
            for (IAtom atom : container.atoms()) {
                if ("H".equals(atom.getSymbol())) {
                    continue;
                }
                counts.merge(atom.getSymbol(), 1, Integer::sum);
            }
        }
        return counts;
    }

    /**
     * Check if a reaction is an identity/transporter (reactants ≡ products).
     * Uses canonical SMILES comparison of each reactant-product pair.
     */
    private boolean isIdentityReaction(IReaction reaction) {
        if (reaction.getReactantCount() != reaction.getProductCount()) {
            return false;
        }
        try {
            org.openscience.cdk.smiles.SmilesGenerator sg = new org.openscience.cdk.smiles.SmilesGenerator(
                    org.openscience.cdk.smiles.SmiFlavor.Canonical);
            java.util.Set<String> reactantSmiles = new java.util.TreeSet<>();
            java.util.Set<String> productSmiles = new java.util.TreeSet<>();
            for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                reactantSmiles.add(sg.create(ac));
            }
            for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                productSmiles.add(sg.create(ac));
            }
            return reactantSmiles.equals(productSmiles);
        } catch (Exception e) {
            return false;
        }
    }

    private boolean hasRingSystems(IReaction reaction) {
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            if (hasRingSystems(ac)) {
                return true;
            }
        }
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            if (hasRingSystems(ac)) {
                return true;
            }
        }
        return false;
    }

    private boolean hasRingSystems(IAtomContainer container) {
        if (container == null || container.getAtomCount() < 3 || container.getBondCount() < 3) {
            return false;
        }
        for (IAtom atom : container.atoms()) {
            if (atom.isInRing() || atom.isAromatic()) {
                return true;
            }
        }
        try {
            return Cycles.sssr(container).numberOfCycles() > 0;
        } catch (Exception e) {
            return false;
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

    /**
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    static class MappingThread implements java.util.concurrent.Callable<Reactor> {

        private static final ILoggingTool MT_LOGGER = createLoggingTool(MappingThread.class);

        private final IReaction cleanedReaction;
        private final IMappingAlgorithm algorithm;
        private final boolean removeHydrogen;

        MappingThread(String message, IReaction cleanedReaction,
                IMappingAlgorithm algorithm, boolean removeHydrogen) {
            this.cleanedReaction = cleanedReaction;
            this.algorithm = algorithm;
            this.removeHydrogen = removeHydrogen;
            MT_LOGGER.info("|++++++++++++++++++++++++++++|");
            MT_LOGGER.info("|Atom Atom Mapping Tool Initialized for " + message);
        }

        @Override
        public Reactor call() throws Exception {
            try {
                Reactor reactor;
                reactor = new Reactor(cleanedReaction, removeHydrogen, algorithm);
                MT_LOGGER.info("|Done " + reactor.getAlgorithm() + " |");
                return reactor;
            } catch (Exception ex) {
                throw ex;
            }
        }
    }

    private static final class MappingThreadFactory implements ThreadFactory {

        @Override
        public Thread newThread(Runnable runnable) {
            Thread thread = new Thread(runnable, "rdt-mapping");
            thread.setDaemon(true);
            return thread;
        }
    }

}
