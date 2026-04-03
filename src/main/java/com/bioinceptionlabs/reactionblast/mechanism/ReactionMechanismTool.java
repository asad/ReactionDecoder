/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.Serializable;
import static java.lang.Integer.MIN_VALUE;

import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.unmodifiableCollection;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAllAtomContainers;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAtomCount;
import org.openscience.smsd.BondEnergies;
import static org.openscience.smsd.BondEnergies.getInstance;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.IFeature;
import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import com.bioinceptionlabs.reactionblast.mapping.CallableAtomMappingTool;
import com.bioinceptionlabs.reactionblast.mapping.Reactor;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.USER_DEFINED;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import static java.lang.Integer.parseInt;
import static java.lang.Math.abs;
import static java.lang.System.getProperty;

import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getAtomArray;
import org.openscience.smsd.ExtAtomContainerManipulator;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ReactionMechanismTool implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private static final String BENCHMARK_ATOM_ID = "benchmarkAtomId";
    private static final String SOURCE_ATOM_ID = "sourceAtomId";
    private static final int DEFAULT_FULL_SCORING_CANDIDATES = 3;
    private static final int MAX_FULL_SCORING_CANDIDATES = 4;
    private final static ILoggingTool LOGGER
            = createLoggingTool(ReactionMechanismTool.class);
    private static final long serialVersionUID = 07342630505L;
    private MappingSolution selectedMapping;
    private Collection<MappingSolution> allSolutions;
    private final boolean accept_no_change;

    // ---- Toolkit-agnostic constructors (ReactionGraph) ----

    /**
     * Toolkit-agnostic constructor. Pass a ReactionGraph from any toolkit.
     *
     * @param reactionGraph toolkit-agnostic reaction graph
     * @param forcedMapping overwrite any existing mapping
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     * @param checkComplex check complex mapping like ring systems
     * @param accept_no_change accept no bond change (transporter reactions)
     * @param standardizer reaction standardizer
     * @throws Exception
     */
    public ReactionMechanismTool(
            com.bioinceptionlabs.reactionblast.model.ReactionGraph reactionGraph,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D,
            boolean checkComplex,
            boolean accept_no_change,
            StandardizeReaction standardizer) throws Exception {
        this(com.bioinceptionlabs.reactionblast.cdk.CDKAdapter.toCDK(reactionGraph),
                forcedMapping, generate2D, generate3D, checkComplex, accept_no_change, standardizer);
    }

    /**
     * Toolkit-agnostic constructor with defaults.
     *
     * @param reactionGraph toolkit-agnostic reaction graph
     * @param forcedMapping overwrite any existing mapping
     * @param checkComplex check complex mapping like ring systems
     * @throws Exception
     */
    public ReactionMechanismTool(
            com.bioinceptionlabs.reactionblast.model.ReactionGraph reactionGraph,
            boolean forcedMapping,
            boolean checkComplex) throws Exception {
        this(com.bioinceptionlabs.reactionblast.cdk.CDKAdapter.toCDK(reactionGraph),
                forcedMapping, true, false, checkComplex, false, new StandardizeReaction());
    }

    // ---- CDK constructors (backward compatible) ----

    /**
     *
     * @param reaction
     * @param forcedMapping force re-mapping of the reactions
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     * @param checkComplex check complex mapping like rings systems
     * rearrangement
     * @throws Exception
     */
    public ReactionMechanismTool(IReaction reaction,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D,
            boolean checkComplex) throws Exception {
        this(reaction,
                forcedMapping,
                generate2D,
                generate3D,
                checkComplex,
                false,
                new StandardizeReaction());
    }

    /**
     *
     * @param reaction
     * @param forcedMapping force re-mapping of the reactions
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     *
     * @param checkComplex check complex mapping like rings systems
     * @param accept_no_change accept no bond change, transporter reactions
     * rearrangement
     * @throws Exception
     */
    public ReactionMechanismTool(IReaction reaction,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D,
            boolean checkComplex,
            boolean accept_no_change
    ) throws Exception {
        this(reaction,
                forcedMapping,
                generate2D,
                generate3D,
                checkComplex,
                accept_no_change,
                new StandardizeReaction());
    }

    /**
     *
     * @param reaction CDK reaction object
     * @param forcedMapping overwrite any existing mapping
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     * @param checkComplex check complex mapping like rings systems
     * @param accept_no_change accept no bond change, transporter reactions
     * @param standardizer standardize reaction
     * @throws CDKException
     * @throws AssertionError
     * @throws Exception
     */
    @SuppressWarnings("deprecation")
    public ReactionMechanismTool(IReaction reaction,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D,
            boolean checkComplex,
            boolean accept_no_change,
            StandardizeReaction standardizer) throws CDKException, AssertionError, Exception {
        if (reaction == null) {
            throw new IllegalArgumentException("Reaction cannot be null");
        }
        this.allSolutions = new ArrayList<>();
        this.selectedMapping = null;
        this.accept_no_change = accept_no_change;//transporter reactions

        if (reaction.getReactantCount() == 0 || reaction.getProductCount() == 0) {
            LOGGER.warn("Reaction has no reactants or no products: {0}", reaction.getID());
        }

        /*
         * IMP: Set all null hydrogen counts to 0, else CDKToBeam cries out loudly
         */
        for (IAtomContainer a : reaction.getReactants().atomContainers()) {
            ExtAtomContainerManipulator.setNullHCountToZero(a);
        }
        /*
         * IMP: Set all null hydrogen counts to 0, else CDKToBeam cries out loudly
         */
        for (IAtomContainer a : reaction.getProducts().atomContainers()) {
            ExtAtomContainerManipulator.setNullHCountToZero(a);
        }

        if (!isBalanced(reaction)) {
            LOGGER.info("Atoms not balanced in the input reaction: {0}; "
                    + "unbalanced reaction may result in erroneous bond change assumptions!", reaction.getID());
            if (!forcedMapping) {
                return;
            }
        }
        if (!forcedMapping && reaction.getFlag(MAPPED)
                && getAtomCount(reaction.getReactants())
                == reaction.getMappingCount()) {
            try {
                LOGGER.info("Using user defined mappings!");
                /*
                 Set Atom IDs
                 */
                for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                    for (IAtom a : ac.atoms()) {
                        a.setID(a.getProperty(ATOM_ATOM_MAPPING) + "");
                    }
                }
                for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                    for (IAtom a : ac.atoms()) {
                        a.setID(a.getProperty(ATOM_ATOM_MAPPING) + "");
                    }
                }
                /*
                 Set mapped flags
                 */
                for (IMapping map : reaction.mappings()) {
                    if (map.getChemObject(0) != null && map.getChemObject(1) != null) {
                        map.getChemObject(0).setFlag(MAPPED, true);
                        map.getChemObject(1).setFlag(MAPPED, true);
                    }
                }
                boolean selected = isMappingSolutionAcceptable(null, USER_DEFINED,
                        reaction, generate2D, generate3D);
                LOGGER.info("is solution: " + USER_DEFINED + " selected: " + selected);
            } catch (Exception e) {
                LOGGER.error(SEVERE, null, e);
                throw new CDKException(NEW_LINE + "ERROR: Unable to calculate bond changes: " + e.getMessage());
            }
        } else {
            try {
                boolean onlyCoreMappingByMCS = true;
                CallableAtomMappingTool amt = new CallableAtomMappingTool(reaction, standardizer,
                        onlyCoreMappingByMCS, checkComplex);
                Map<IMappingAlgorithm, Reactor> solutions = amt.getSolutions();
                List<EvaluationCandidate> orderedSolutions = orderSolutionsForEvaluation(solutions);
                List<EvaluationCandidate> candidates = collectCandidatesForEvaluation(orderedSolutions);

                LOGGER.debug("!!!!Calculating Best Mapping Model!!!!");
                for (MappingSolution mappingSolution : computeMappingSolutions(candidates,
                        generate2D, generate3D)) {
                    LOGGER.debug("===considerMappingSolution===");
                    boolean selected = considerMappingSolution(mappingSolution);
                    LOGGER.debug("is solution: " + mappingSolution.getAlgorithmID()
                            + " selected: " + selected);
                }
            } catch (Exception e) {
                LOGGER.error(SEVERE, "Bond change calculation error", e);
                throw new Exception(NEW_LINE + "ERROR: Unable to calculate bond changes: " + e.getMessage(), e);
            }
            if (this.selectedMapping != null) {
                LOGGER.info("Selected algorithm: " + this.selectedMapping.getAlgorithmID().description()
                        + " (bonds=" + this.selectedMapping.getTotalBondChanges()
                        + ", energy=" + this.selectedMapping.getBondEnergySum()
                        + ", fragments=" + this.selectedMapping.getTotalFragmentChanges() + ")");
            }
            LOGGER.debug("=====DONE REACTION MECH TOOL=====");
        }
    }

    /**
     * Check if a reaction is balanced (same heavy atom counts on both sides).
     * Hydrogens are excluded since they are often implicit.
     */
    private boolean isBalanced(IReaction r) {
        Map<String, Integer> reactantAtoms = countHeavyAtoms(r.getReactants());
        Map<String, Integer> productAtoms = countHeavyAtoms(r.getProducts());

        if (!reactantAtoms.equals(productAtoms)) {
            LOGGER.debug("Number of atom(s) on the Left side "
                    + reactantAtoms.values().stream().mapToInt(Integer::intValue).sum()
                    + " =/= Number of atom(s) on the Right side "
                    + productAtoms.values().stream().mapToInt(Integer::intValue).sum());
            LOGGER.debug(reactantAtoms + " =/= " + productAtoms);
            return false;
        }
        return true;
    }

    private Map<String, Integer> countHeavyAtoms(IAtomContainerSet containers) {
        Map<String, Integer> counts = new TreeMap<>();
        for (IAtomContainer mol : containers.atomContainers()) {
            for (IAtom a : mol.atoms()) {
                if (!"H".equals(a.getSymbol())) {
                    counts.merge(a.getSymbol(), 1, Integer::sum);
                }
            }
        }
        return counts;
    }

    private boolean isMappingSolutionAcceptable(Reactor reactor,
            IMappingAlgorithm ma,
            IReaction reaction,
            boolean generate2D,
            boolean generate3D
    ) throws Exception {
        if (reactor != null && reactor.getMappingCount() > 500) {
            LOGGER.warn("Large mapping: " + reactor.getMappingCount()
                    + " atoms — bond change computation may be slow");
        }
        boolean chosen = false;
        try {
            BondChangeCalculator bcc;
            int fragmentDeltaChanges;
            if (reactor == null && ma.equals(USER_DEFINED)) {
                bcc = new BondChangeCalculator(reaction);
                bcc.computeBondChanges(generate2D, generate3D);
                fragmentDeltaChanges = bcc.getTotalFragmentCount();
                int bondChange = (int) getTotalBondChange(bcc.getFormedCleavedWFingerprint());
                bondChange += getTotalBondChange(bcc.getOrderChangesWFingerprint());
                int stereoChanges = (int) getTotalBondChange(bcc.getStereoChangesWFingerprint());
                boolean skipHydrogenRealtedBondChanges = true;
                int bondBreakingEnergy = getTotalBondChangeEnergy(bcc.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);
                int totalSmallestFragmentCount = bcc.getTotalSmallestFragmentSize();
                int totalCarbonBondChanges = getTotalCarbonBondChange(bcc.getFormedCleavedWFingerprint());
                int localScore = bondChange + fragmentDeltaChanges;
                MappingSolution mappingSolution = new MappingSolution(
                        bcc,
                        ma,
                        bcc.getReaction(),
                        reactor,
                        bondBreakingEnergy,
                        totalCarbonBondChanges,
                        bondChange,
                        fragmentDeltaChanges,
                        stereoChanges,
                        totalSmallestFragmentCount,
                        localScore,
                        bcc.getEnergyDelta()
                );

                chosen = true;
                mappingSolution.setChosen(chosen);
                this.selectedMapping = mappingSolution;
                this.allSolutions.add(mappingSolution);
            } else {
                if (reactor == null) {
                    throw new CDKException("Reactor is NULL");
                }

                bcc = new BondChangeCalculator(reactor.getReactionWithAtomAtomMapping());
                bcc.computeBondChanges(generate2D, generate3D);
                fragmentDeltaChanges = bcc.getTotalFragmentCount() + reactor.getDelta();

                int bondCleavedFormed = (int) getTotalBondChange(bcc.getFormedCleavedWFingerprint());
                int bondChange = bondCleavedFormed;
                bondChange += getTotalBondChange(bcc.getOrderChangesWFingerprint());
                int stereoChanges = (int) getTotalBondChange(bcc.getStereoChangesWFingerprint());
                boolean skipHydrogenRealtedBondChanges = true;
                int bondBreakingEnergy = getTotalBondChangeEnergy(bcc.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);
                int totalSmallestFragmentCount = bcc.getTotalSmallestFragmentSize();
                int totalCarbonBondChanges = getTotalCarbonBondChange(bcc.getFormedCleavedWFingerprint());
                int localScore = bondChange + fragmentDeltaChanges;
                LOGGER.info(
                        "Score: " + fragmentDeltaChanges + " : " + bondChange);
                LOGGER.info(
                        ", Energy Barrier: " + bondBreakingEnergy);
                LOGGER.info(
                        ", Energy Delta: " + bcc.getEnergyDelta());

                bcc.getReaction().setFlag(MAPPED, true);

                MappingSolution mappingSolution = new MappingSolution(
                        bcc,
                        ma,
                        bcc.getReaction(),
                        reactor,
                        bondBreakingEnergy,
                        totalCarbonBondChanges,
                        bondChange,
                        fragmentDeltaChanges,
                        stereoChanges,
                        totalSmallestFragmentCount,
                        localScore,
                        bcc.getEnergyDelta());

                if (ma == null) {
                    throw new CDKException("Model is pointing to NULL");
                }
                LOGGER.info("MA: " + ma.description());
                boolean changeFeasible = isChangeFeasible(mappingSolution);
                if (changeFeasible) {
                    chosen = changeFeasible;
                    mappingSolution.setChosen(changeFeasible);
                    this.selectedMapping = mappingSolution;
                }
                this.allSolutions.add(mappingSolution);
            }
        } catch (Exception e) {
            LOGGER.error(SEVERE, "isMappingSolutionAcceptable failed", e);
            throw new Exception(NEW_LINE + "ERROR: Unable to calculate bond changes: " + e.getMessage(), e);
        }
        return chosen;
    }

    /**
     * Determines if a new mapping solution should replace the current best.
     *
     * Chemical rationale (Occam's razor for reaction mapping):
     * 1. Prefer fewer total bond changes (simplest mechanism)
     * 2. Among equal bond changes, prefer fewer fragment changes (less molecular rearrangement)
     * 3. Among equal fragments, prefer lower bond energy change (thermodynamic favorability)
     * 4. Use stereo changes and carbon bond changes as tiebreakers
     *
     * Special cases: transporters (no bond changes), identity reactions,
     * and reactions with only stereo changes are handled first.
     */
    private boolean isChangeFeasible(MappingSolution ms) {

        /*
        * This condition is valuble to trace graph isomorphism as only min algorithm checks this change. the idea is to
        * assume a change if rest of the algorithm detects no change.
        *
        * TODO: check what is the impact if this logic if there are only stereo changes in a reaction.
         */
        if (this.selectedMapping != null) {
            LOGGER.debug(NEW_LINE + " selectedMapping.getAlgorithmID().description() " + selectedMapping.getAlgorithmID().description());
            LOGGER.debug(" selectedMapping.getTotalStereoChanges() " + selectedMapping.getTotalStereoChanges());
            LOGGER.debug(" selectedMapping.getTotalBondChanges() " + selectedMapping.getTotalBondChanges());
            LOGGER.debug(" selectedMapping.getSmallestFragmentCount() " + selectedMapping.getSmallestFragmentCount());
            LOGGER.debug(" selectedMapping.getBondEnergyChange() " + selectedMapping.getBondEnergySum());
            LOGGER.debug(" selectedMapping.getTotalFragmentChanges() " + selectedMapping.getTotalFragmentChanges());
            LOGGER.debug(" ms.getTotalChanges() " + selectedMapping.getTotalChanges());
            LOGGER.debug(" Total Carbon Bond Changes " + selectedMapping.getTotalCarbonBondChanges());
        }
        LOGGER.debug(NEW_LINE + " ms.getAlgorithmID().description() " + ms.getAlgorithmID().description());
        LOGGER.debug(" ms.getTotalStereoChanges() " + ms.getTotalStereoChanges());
        LOGGER.debug(" ms.getTotalBondChanges() " + ms.getTotalBondChanges());
        LOGGER.debug(" ms.getSmallestFragmentCount() " + ms.getSmallestFragmentCount());
        LOGGER.debug(" ms.getBondEnergyChange() " + ms.getBondEnergySum());
        LOGGER.debug(" ms.getTotalFragmentChanges() " + ms.getTotalFragmentChanges());
        LOGGER.debug(" ms.getTotalChanges() " + ms.getTotalChanges());
        LOGGER.debug(" Total Carbon Bond Changes " + ms.getTotalCarbonBondChanges());

        /*
         * only transporter reactions where we expect no bond change
         */
        if (this.selectedMapping != null
                && this.accept_no_change == true
                && ms.getTotalBondChanges() == 0
                && ms.getTotalStereoChanges() == 0) {
            LOGGER.debug("CASE: Transporter");
            return true;
        }

        /*
         * This is to skip reaction where the no change is detected.
         * Example: R02996
         */
        if (this.selectedMapping != null
                && ms.getTotalBondChanges() == 0
                && ms.getTotalStereoChanges() == 0) {
            return false;
        }

        /*
        * if this is the first solution then accept it
         */
        if (this.selectedMapping == null) {
            LOGGER.info("Condition Default " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition Default");
            return true;
        } else if (ms.getBondEnergySum() == 0.
                && ms.getTotalFragmentChanges() == 0
                && ms.getTotalBondChanges() == 0
                && this.selectedMapping.getTotalStereoChanges() >= ms.getTotalStereoChanges()) {
            LOGGER.info("Condition 1 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 1");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                || this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum())) {
            LOGGER.info("Condition 2 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 2");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > 0
                && ms.getTotalFragmentChanges() > 0) {
            LOGGER.info("Condition 3 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 3");
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() >= ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() >= ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() >= ms.getTotalCarbonBondChanges()) {
            /* Example reaction R05069*/
            LOGGER.info("Condition 4 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 4");
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            LOGGER.info("Condition 5 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 5");
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() == ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() >= ms.getTotalCarbonBondChanges()) {
            LOGGER.info("Condition 6 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 6");
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            LOGGER.info("Condition 7 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 7");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()) {
            /* && this.selectedMapping.getBondEnergyChange() > ms.getBondEnergyChange()) {*/
            LOGGER.info("Condition 8 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 8");
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()) {
            LOGGER.info("Condition 9 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 9");
            return true;
        } else if (this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalStereoChanges() > ms.getTotalStereoChanges()) {
            LOGGER.info("Condition 10 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 10");
            return true;
        } else if (this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()) {
            LOGGER.info("Condition 11 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 11");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() < ms.getTotalBondChanges()
                && this.selectedMapping.getBondEnergySum() < ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: R00652 N Vs O exchange*/
            LOGGER.info("Condition 12 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 12");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: Rhea 22881 N Vs O exchange*/
            LOGGER.info("Condition 13 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 13");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() < ms.getSmallestFragmentCount()) {
            /*
            * Rhea  reaction RHEA:20301 bigger fragment preferred 
             */
            LOGGER.info("Condition 14 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 14");
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            /*This condition is for reactions like:
            R05421 (O-P over O-C)
            **/
            LOGGER.info("Condition 15 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 15");
            return true;
        } else if (hasEquivalentSelectionScore(this.selectedMapping, ms)
                && hasPreferredCanonicalMapping(ms, this.selectedMapping)) {
            LOGGER.info("Condition 16 " + ms.getAlgorithmID().description());
            LOGGER.debug("CASE: Condition 16");
            return true;
        }
        LOGGER.debug("CASE: FAILED");
        return false;
    }

    private double getTotalBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return total;
    }

    private int getTotalCarbonBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().filter((key) -> (key.getPattern().contains("C-C")
                || key.getPattern().contains("C=C")
                || key.getPattern().contains("C#C")
                || key.getPattern().contains("C%C")
                || key.getPattern().contains("C@C"))).map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return (int) total;
    }

    private int getTotalBondChangeEnergy(IPatternFingerprinter fingerprint, boolean skipHydrogen) {
        int total = 0;
        try {
            BondEnergies be = getInstance();
            for (IFeature feature : fingerprint.getFeatures()) {
                double val = feature.getWeight();
                String key = feature.getPattern();
                if (val > 0) {
//                    LOGGER.debug("BOND BROKEN/FORMED: " + key + " : " + val);
                    if (key.contains("-") || key.contains("%") || key.contains("@")) {
                        String[] temp = null;
                        if (key.contains("-")) {
                            temp = key.split("-");
                        } else if (key.contains("%")) {
                            temp = key.split("%");
                        } else if (key.contains("@")) {
                            temp = key.split("@");
                        }
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        if (key.contains("%")) {
                            /*
                             * Aromatic bond (~1.5 order): use average of single + double
                             * bond energies as a chemically sound approximation.
                             */
                            int eSingle = be.getEnergies(temp[0], temp[1], SINGLE);
                            int eDouble = be.getEnergies(temp[0], temp[1], DOUBLE);
                            int energy;
                            if (eSingle > 0 && eDouble > 0) {
                                energy = (eSingle + eDouble) / 2;
                            } else if (eSingle > 0) {
                                energy = eSingle;
                            } else {
                                energy = eDouble > 0 ? eDouble : 0;
                            }
                            if (energy > 0) {
                                total += val * energy;
                            }
                        } else {
                            int energy = be.getEnergies(temp[0], temp[1], SINGLE);
                            if (energy > 0) {
                                total += val * energy;
                            }
                        }
                    } else if (key.contains("=")) {
                        String[] temp = key.split("=");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], DOUBLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("#")) {
                        String[] temp = key.split("#");
                        int energy = be.getEnergies(temp[0], temp[1], TRIPLE);
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("$")) {
                        String[] temp = key.split("\\$");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], QUADRUPLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    }
                }
            }
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return abs(total);
    }

    /**
     *
     * @return
     */
    public String getMappingDescription() {
        return this.selectedMapping.toString();
    }

    /**
     *
     * @return
     */
    public MappingSolution getSelectedSolution() {
        return this.selectedMapping;
    }

    /**
     *
     * @return
     */
    public Collection<MappingSolution> getAllSolutions() {
        return unmodifiableCollection(this.allSolutions);
    }

    /**
     * Get the mapped reaction as a toolkit-agnostic ReactionGraph.
     *
     * @return ReactionGraph with atom-atom mapping, or null if no solution
     */
    public com.bioinceptionlabs.reactionblast.model.ReactionGraph getMappedReactionGraph() {
        if (selectedMapping == null || selectedMapping.getReactor() == null) {
            return null;
        }
        try {
            IReaction mapped = selectedMapping.getReactor().getReactionWithAtomAtomMapping();
            return mapped != null ? com.bioinceptionlabs.reactionblast.cdk.CDKAdapter.fromCDK(mapped) : null;
        } catch (Exception e) {
            LOGGER.error(SEVERE, "Failed to get mapped reaction graph", e);
            return null;
        }
    }

    private List<EvaluationCandidate> orderSolutionsForEvaluation(
            Map<IMappingAlgorithm, Reactor> solutions) {
        List<EvaluationCandidate> ordered = snapshotCandidates(solutions);
        ordered.sort(evaluationCandidateComparator(isIdentityLike(ordered)));
        return ordered;
    }

    private List<EvaluationCandidate> collectCandidatesForEvaluation(
            List<EvaluationCandidate> orderedSolutions) {
        List<EvaluationCandidate> candidates = new ArrayList<>();
        Map<String, EvaluationCandidate> uniqueCandidates = new LinkedHashMap<>();

        for (EvaluationCandidate candidate : orderedSolutions) {
            if (!candidate.coverage.isComplete() || !candidate.coverage.isBalancedMapped()) {
                LOGGER.debug("Unmapped atoms present in this reaction" + "(" + candidate.algorithm + ") algorithm.");
            }
            if (shouldSkipInferiorCoverage(candidate.coverage)) {
                LOGGER.debug("Skipping " + candidate.algorithm + " scoring due to inferior mapping coverage");
                continue;
            }

            String dedupeKey = candidate.coverage.getMappedAtoms()
                    + ":" + candidate.coverage.getUnmappedAtoms()
                    + ":" + candidate.signature;
            if (uniqueCandidates.containsKey(dedupeKey)) {
                LOGGER.debug("Skipping duplicate mapping candidate from " + candidate.algorithm
                        + " equivalent to " + uniqueCandidates.get(dedupeKey).algorithm);
                continue;
            }

            uniqueCandidates.put(dedupeKey, candidate);
            candidates.add(candidate);
        }
        return limitCandidatesForFullScoring(candidates, isIdentityLike(candidates));
    }

    private List<EvaluationCandidate> snapshotCandidates(Map<IMappingAlgorithm, Reactor> solutions) {
        List<EvaluationCandidate> candidates = new ArrayList<>(solutions.size());
        for (Map.Entry<IMappingAlgorithm, Reactor> entry : solutions.entrySet()) {
            IMappingAlgorithm algorithm = entry.getKey();
            Reactor reactor = entry.getValue();
            if (reactor == null) {
                LOGGER.warn("Reactor is NULL");
                continue;
            }
            try {
                IReaction mappedReaction = reactor.getReactionWithAtomAtomMapping();
                MappingCoverage coverage = summarizeCoverage(mappedReaction);
                String signature = canonicalMappingSignature(mappedReaction);
                QuickScore quickScore = estimateQuickScore(mappedReaction, reactor);
                candidates.add(new EvaluationCandidate(
                        algorithm, reactor, mappedReaction, coverage, signature, quickScore));
            } catch (Exception ex) {
                LOGGER.debug("Skipping " + algorithm + " due to snapshot failure: " + ex.getMessage());
            }
        }
        return candidates;
    }

    private List<MappingSolution> computeMappingSolutions(List<EvaluationCandidate> candidates,
            boolean generate2D, boolean generate3D) throws Exception {
        if (candidates.isEmpty()) {
            return new ArrayList<>();
        }
        if (candidates.size() == 1) {
            List<MappingSolution> single = new ArrayList<>(1);
            single.add(computeMappingSolution(candidates.get(0), generate2D, generate3D));
            return single;
        }

        int threadCount = Math.min(candidates.size(),
                Math.max(1, Runtime.getRuntime().availableProcessors() - 1));
        ExecutorService executor = Executors.newFixedThreadPool(threadCount);
        try {
            List<Future<MappingSolution>> futures = new ArrayList<>(candidates.size());
            for (EvaluationCandidate candidate : candidates) {
                futures.add(executor.submit(
                        () -> computeMappingSolution(candidate, generate2D, generate3D)));
            }

            List<MappingSolution> evaluated = new ArrayList<>(candidates.size());
            for (Future<MappingSolution> future : futures) {
                evaluated.add(future.get());
            }
            return evaluated;
        } finally {
            executor.shutdownNow();
        }
    }

    @SuppressWarnings("deprecation")
    private MappingSolution computeMappingSolution(EvaluationCandidate candidate,
            boolean generate2D, boolean generate3D) throws Exception {
        Reactor reactor = candidate.reactor;
        if (reactor == null) {
            throw new CDKException("Reactor is NULL");
        }

        if (reactor.getMappingCount() > 500) {
            LOGGER.warn("Large mapping: " + reactor.getMappingCount()
                    + " atoms — bond change computation may be slow");
        }

        BondChangeCalculator bcc = new BondChangeCalculator(candidate.mappedReaction);
        bcc.computeBondChanges(generate2D, generate3D);
        int fragmentDeltaChanges = bcc.getTotalFragmentCount() + reactor.getDelta();

        int bondCleavedFormed = (int) getTotalBondChange(bcc.getFormedCleavedWFingerprint());
        int bondChange = bondCleavedFormed
                + (int) getTotalBondChange(bcc.getOrderChangesWFingerprint());
        int stereoChanges = (int) getTotalBondChange(bcc.getStereoChangesWFingerprint());
        boolean skipHydrogenRealtedBondChanges = true;
        int bondBreakingEnergy = getTotalBondChangeEnergy(
                bcc.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);
        int totalSmallestFragmentCount = bcc.getTotalSmallestFragmentSize();
        int totalCarbonBondChanges = getTotalCarbonBondChange(
                bcc.getFormedCleavedWFingerprint());
        int localScore = bondChange + fragmentDeltaChanges;

        LOGGER.info("Score: " + fragmentDeltaChanges + " : " + bondChange);
        LOGGER.info(", Energy Barrier: " + bondBreakingEnergy);
        LOGGER.info(", Energy Delta: " + bcc.getEnergyDelta());

        bcc.getReaction().setFlag(MAPPED, true);

        return new MappingSolution(
                bcc,
                candidate.algorithm,
                bcc.getReaction(),
                reactor,
                bondBreakingEnergy,
                totalCarbonBondChanges,
                bondChange,
                fragmentDeltaChanges,
                stereoChanges,
                totalSmallestFragmentCount,
                localScore,
                bcc.getEnergyDelta());
    }

    private Comparator<EvaluationCandidate> evaluationCandidateComparator(boolean identityLike) {
        return Comparator
                .<EvaluationCandidate>comparingInt(candidate -> candidate.coverage.isComplete() ? 0 : 1)
                .thenComparingInt(candidate -> candidate.coverage.isBalancedMapped() ? 0 : 1)
                .thenComparingInt(candidate -> candidate.quickScore.totalScore())
                .thenComparingInt(candidate -> candidate.quickScore.bondChangeEstimate)
                .thenComparingInt(candidate -> candidate.quickScore.orderChangeEstimate)
                .thenComparingInt(candidate -> candidate.quickScore.fragmentPenalty)
                .thenComparingInt(candidate -> candidate.quickScore.unmappedBondPenalty)
                .thenComparingInt(candidate -> candidate.quickScore.carbonBondChangeEstimate)
                .thenComparingInt(candidate -> -candidate.quickScore.mappedBondCount)
                .thenComparingInt(candidate -> -candidate.coverage.getMappedAtoms())
                .thenComparingInt(candidate -> candidate.coverage.getUnmappedAtoms())
                .thenComparingInt(candidate -> algorithmPriority(candidate.algorithm, identityLike))
                .thenComparing(candidate -> candidate.signature);
    }

    private boolean isIdentityLike(List<EvaluationCandidate> candidates) {
        return candidates.stream()
                .map(candidate -> candidate.mappedReaction)
                .filter(mappedReaction -> mappedReaction != null)
                .findFirst()
                .map(this::looksLikeIdentityReaction)
                .orElse(false);
    }

    private List<EvaluationCandidate> limitCandidatesForFullScoring(
            List<EvaluationCandidate> candidates,
            boolean identityLike) {
        if (candidates.size() <= DEFAULT_FULL_SCORING_CANDIDATES) {
            return candidates;
        }

        List<EvaluationCandidate> ranked = new ArrayList<>(candidates);
        ranked.sort(evaluationCandidateComparator(identityLike));

        int limit = Math.min(DEFAULT_FULL_SCORING_CANDIDATES, ranked.size());
        if (hasAmbiguousTopTier(ranked)) {
            limit = Math.min(MAX_FULL_SCORING_CANDIDATES, ranked.size());
        }

        List<EvaluationCandidate> retained = new ArrayList<>(ranked.subList(0, limit));
        if (limit < ranked.size()) {
            QuickScore cutoff = ranked.get(limit - 1).quickScore;
            for (int index = limit; index < ranked.size() && retained.size() < MAX_FULL_SCORING_CANDIDATES; index++) {
                EvaluationCandidate candidate = ranked.get(index);
                if (candidate.quickScore.isNear(cutoff)) {
                    retained.add(candidate);
                }
            }
        }

        if (retained.size() < candidates.size()) {
            LOGGER.debug("Reduced full bond-change scoring from "
                    + candidates.size() + " to " + retained.size() + " candidate(s)");
        }
        return retained;
    }

    private boolean hasAmbiguousTopTier(List<EvaluationCandidate> ranked) {
        if (ranked.size() < 2) {
            return false;
        }

        EvaluationCandidate best = ranked.get(0);
        EvaluationCandidate challenger = ranked.get(1);
        return best.coverage.isComplete() == challenger.coverage.isComplete()
                && best.coverage.isBalancedMapped() == challenger.coverage.isBalancedMapped()
                && best.quickScore.hasEquivalentCoreScore(challenger.quickScore)
                && !best.signature.equals(challenger.signature);
    }

    private QuickScore estimateQuickScore(IReaction reaction, Reactor reactor) {
        if (reaction == null) {
            return new QuickScore(Integer.MAX_VALUE, Integer.MAX_VALUE,
                    Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 0);
        }

        Map<String, BondDescriptor> reactantBonds = collectMappedBondDescriptors(reaction.getReactants());
        Map<String, BondDescriptor> productBonds = collectMappedBondDescriptors(reaction.getProducts());
        Set<String> allBondKeys = new TreeSet<>();
        allBondKeys.addAll(reactantBonds.keySet());
        allBondKeys.addAll(productBonds.keySet());

        int bondChangeEstimate = 0;
        int orderChangeEstimate = 0;
        int carbonBondChangeEstimate = 0;
        for (String bondKey : allBondKeys) {
            BondDescriptor reactantBond = reactantBonds.get(bondKey);
            BondDescriptor productBond = productBonds.get(bondKey);
            if (reactantBond == null || productBond == null) {
                bondChangeEstimate++;
                if ((reactantBond != null && reactantBond.carbonOnly)
                        || (productBond != null && productBond.carbonOnly)) {
                    carbonBondChangeEstimate++;
                }
                continue;
            }
            if (!reactantBond.sameType(productBond)) {
                orderChangeEstimate++;
                if (reactantBond.carbonOnly && productBond.carbonOnly) {
                    carbonBondChangeEstimate++;
                }
            }
        }

        int fragmentPenalty = reactor != null ? Math.max(0, reactor.getDelta()) : 0;
        int unmappedBondPenalty = countUnmappedBondPenalty(reaction.getReactants())
                + countUnmappedBondPenalty(reaction.getProducts());
        int mappedBondCount = reactantBonds.size() + productBonds.size();
        return new QuickScore(
                bondChangeEstimate,
                orderChangeEstimate,
                carbonBondChangeEstimate,
                fragmentPenalty,
                unmappedBondPenalty,
                mappedBondCount);
    }

    private Map<String, BondDescriptor> collectMappedBondDescriptors(IAtomContainerSet containers) {
        Map<String, BondDescriptor> descriptors = new LinkedHashMap<>();
        for (IAtomContainer container : containers.atomContainers()) {
            for (IBond bond : container.bonds()) {
                if (bond == null) {
                    continue;
                }
                IAtom begin = bond.getBegin();
                IAtom end = bond.getEnd();
                if (begin == null || end == null) {
                    continue;
                }
                if ("H".equals(begin.getSymbol()) || "H".equals(end.getSymbol())) {
                    continue;
                }

                int beginMap = getAtomMapNumber(begin);
                int endMap = getAtomMapNumber(end);
                if (beginMap <= 0 || endMap <= 0) {
                    continue;
                }

                String key = beginMap < endMap
                        ? beginMap + ":" + endMap
                        : endMap + ":" + beginMap;
                descriptors.put(key, new BondDescriptor(
                        toBondOrderValue(bond),
                        bond.isAromatic(),
                        "C".equals(begin.getSymbol()) && "C".equals(end.getSymbol())));
            }
        }
        return descriptors;
    }

    private int countUnmappedBondPenalty(IAtomContainerSet containers) {
        int penalty = 0;
        for (IAtomContainer container : containers.atomContainers()) {
            for (IBond bond : container.bonds()) {
                if (bond == null) {
                    continue;
                }
                IAtom begin = bond.getBegin();
                IAtom end = bond.getEnd();
                if (begin == null || end == null) {
                    continue;
                }
                if ("H".equals(begin.getSymbol()) && "H".equals(end.getSymbol())) {
                    continue;
                }
                if (getAtomMapNumber(begin) <= 0 || getAtomMapNumber(end) <= 0) {
                    penalty++;
                }
            }
        }
        return penalty;
    }

    private int toBondOrderValue(IBond bond) {
        if (bond == null || bond.getOrder() == null) {
            return 0;
        }
        switch (bond.getOrder()) {
            case SINGLE:
                return 1;
            case DOUBLE:
                return 2;
            case TRIPLE:
                return 3;
            case QUADRUPLE:
                return 4;
            default:
                return 0;
        }
    }

    private boolean considerMappingSolution(MappingSolution mappingSolution) throws Exception {
        if (mappingSolution == null) {
            return false;
        }
        if (mappingSolution.getAlgorithmID() == null) {
            throw new CDKException("Model is pointing to NULL");
        }

        LOGGER.info("MA: " + mappingSolution.getAlgorithmID().description());
        boolean changeFeasible = isChangeFeasible(mappingSolution);
        if (changeFeasible) {
            if (this.selectedMapping != null) {
                this.selectedMapping.setChosen(false);
            }
            mappingSolution.setChosen(true);
            this.selectedMapping = mappingSolution;
        }
        this.allSolutions.add(mappingSolution);
        return changeFeasible;
    }

    private int algorithmPriority(IMappingAlgorithm algorithm, boolean identityLike) {
        if (algorithm == USER_DEFINED) {
            return -1;
        }
        if (identityLike) {
            switch (algorithm) {
                case MIN:
                    return 0;
                case RINGS:
                    return 1;
                case MAX:
                    return 2;
                case MIXTURE:
                    return 3;
                default:
                    return 4;
            }
        }
        switch (algorithm) {
            case RINGS:
                return 0;
            case MIN:
                return 1;
            case MAX:
                return 2;
            case MIXTURE:
                return 3;
            default:
                return 4;
        }
    }

    private boolean looksLikeIdentityReaction(IReaction reaction) {
        if (reaction == null || reaction.getReactantCount() != reaction.getProductCount()) {
            return false;
        }
        try {
            SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Canonical);
            List<String> reactants = new ArrayList<>();
            List<String> products = new ArrayList<>();
            for (IAtomContainer reactant : reaction.getReactants().atomContainers()) {
                reactants.add(smilesGenerator.create(reactant));
            }
            for (IAtomContainer product : reaction.getProducts().atomContainers()) {
                products.add(smilesGenerator.create(product));
            }
            reactants.sort(String::compareTo);
            products.sort(String::compareTo);
            return reactants.equals(products);
        } catch (CDKException e) {
            return false;
        }
    }

    private MappingCoverage summarizeCoverage(IReaction reaction) {
        if (reaction == null) {
            return new MappingCoverage(0, 0, 0, 0);
        }
        return new MappingCoverage(
                getTotalNonHydrogenAtomCount(reaction.getReactants()),
                getTotalNonHydrogenAtomCount(reaction.getProducts()),
                getMappedNonHydrogenAtomCount(reaction.getReactants()),
                getMappedNonHydrogenAtomCount(reaction.getProducts()));
    }

    private boolean shouldSkipInferiorCoverage(MappingCoverage candidateCoverage) {
        if (selectedMapping == null || selectedMapping.getReaction() == null) {
            return false;
        }
        MappingCoverage selectedCoverage = summarizeCoverage(selectedMapping.getReaction());
        if (selectedCoverage.isComplete() && selectedCoverage.isBalancedMapped()) {
            return !candidateCoverage.isComplete()
                    || !candidateCoverage.isBalancedMapped()
                    || candidateCoverage.getMappedAtoms() < selectedCoverage.getMappedAtoms();
        }
        return false;
    }

    private boolean hasEquivalentSelectionScore(MappingSolution selected, MappingSolution candidate) {
        if (selected == null || candidate == null) {
            return false;
        }

        MappingCoverage selectedCoverage = summarizeCoverage(selected.getReaction());
        MappingCoverage candidateCoverage = summarizeCoverage(candidate.getReaction());

        return selected.getTotalBondChanges() == candidate.getTotalBondChanges()
                && selected.getTotalFragmentChanges() == candidate.getTotalFragmentChanges()
                && selected.getTotalStereoChanges() == candidate.getTotalStereoChanges()
                && selected.getSmallestFragmentCount() == candidate.getSmallestFragmentCount()
                && selected.getTotalCarbonBondChanges() == candidate.getTotalCarbonBondChanges()
                && selected.getTotalChanges() == candidate.getTotalChanges()
                && Double.compare(selected.getBondEnergySum(), candidate.getBondEnergySum()) == 0
                && Double.compare(selected.getEnergyDelta(), candidate.getEnergyDelta()) == 0
                && selectedCoverage.getMappedAtoms() == candidateCoverage.getMappedAtoms()
                && selectedCoverage.isComplete() == candidateCoverage.isComplete()
                && selectedCoverage.isBalancedMapped() == candidateCoverage.isBalancedMapped();
    }

    private boolean hasPreferredCanonicalMapping(MappingSolution candidate, MappingSolution selected) {
        String candidateSignature = canonicalMappingSignature(candidate.getReaction());
        String selectedSignature = canonicalMappingSignature(selected.getReaction());
        return !candidateSignature.isEmpty()
                && (selectedSignature.isEmpty() || candidateSignature.compareTo(selectedSignature) < 0);
    }

    private String canonicalMappingSignature(IReaction reaction) {
        if (reaction == null) {
            return "";
        }

        Map<Integer, String> reactantPositions = new TreeMap<>();
        Map<Integer, String> productPositions = new TreeMap<>();
        collectMappedAtomPositions(reaction.getReactants(), "R", reactantPositions);
        collectMappedAtomPositions(reaction.getProducts(), "P", productPositions);

        StringBuilder signature = new StringBuilder();
        for (Map.Entry<Integer, String> entry : reactantPositions.entrySet()) {
            String productPosition = productPositions.get(entry.getKey());
            if (productPosition == null) {
                continue;
            }
            if (signature.length() > 0) {
                signature.append('|');
            }
            signature.append(entry.getValue()).append('>').append(productPosition);
        }
        return signature.toString();
    }

    private void collectMappedAtomPositions(IAtomContainerSet containers, String side,
            Map<Integer, String> positions) {
        int moleculeIndex = 0;
        for (IAtomContainer molecule : containers.atomContainers()) {
            for (int atomIndex = 0; atomIndex < molecule.getAtomCount(); atomIndex++) {
                IAtom atom = molecule.getAtom(atomIndex);
                int mappingNumber = getAtomMapNumber(atom);
                if (mappingNumber > 0) {
                    positions.put(mappingNumber, getStableAtomPosition(atom, side, moleculeIndex, atomIndex));
                }
            }
            moleculeIndex++;
        }
    }

    private String getStableAtomPosition(IAtom atom, String side, int moleculeIndex, int atomIndex) {
        if (atom == null) {
            return side + ":" + moleculeIndex + ":" + atomIndex;
        }
        Object benchmarkAtomId = atom.getProperty(BENCHMARK_ATOM_ID);
        if (benchmarkAtomId != null) {
            return benchmarkAtomId.toString();
        }
        Object sourceAtomId = atom.getProperty(SOURCE_ATOM_ID);
        if (sourceAtomId != null) {
            return sourceAtomId.toString();
        }
        return side + ":" + moleculeIndex + ":" + atomIndex;
    }

    private int getAtomMapNumber(IAtom atom) {
        if (atom == null) {
            return 0;
        }
        if (atom.getMapIdx() > 0) {
            return atom.getMapIdx();
        }

        Object atomAtomMapping = atom.getProperty(ATOM_ATOM_MAPPING);
        if (atomAtomMapping instanceof Integer && (Integer) atomAtomMapping > 0) {
            return (Integer) atomAtomMapping;
        }
        if (atomAtomMapping != null) {
            try {
                int parsed = parseInt(atomAtomMapping.toString());
                if (parsed > 0) {
                    return parsed;
                }
            } catch (NumberFormatException ignore) {
                // Fall through to the legacy property.
            }
        }

        Object legacyMapNumber = atom.getProperty("molAtomMapNumber");
        if (legacyMapNumber instanceof Integer && (Integer) legacyMapNumber > 0) {
            return (Integer) legacyMapNumber;
        }
        if (legacyMapNumber != null) {
            try {
                int parsed = parseInt(legacyMapNumber.toString());
                return parsed > 0 ? parsed : 0;
            } catch (NumberFormatException ignore) {
                return 0;
            }
        }
        return 0;
    }

    @SuppressWarnings("deprecation")
    private int getMappedNonHydrogenAtomCount(IAtomContainerSet mol) {
        List<IAtomContainer> allAtomContainers = getAllAtomContainers(mol);
        int count = 0;
        for (IAtomContainer ac : allAtomContainers) {
            IAtom[] atomArray = getAtomArray(ac);
            for (IAtom atom : atomArray) {
                if (atom.getSymbol().equalsIgnoreCase("H")) {
                    continue;
                }
                Object atomAtomMapping = atom.getProperty(ATOM_ATOM_MAPPING);
                if (atom.getFlag(MAPPED)) {
                    count++;
                } else if (atomAtomMapping instanceof Integer && (Integer) atomAtomMapping > 0) {
                    count++;
                } else if (atomAtomMapping != null) {
                    try {
                        if (parseInt(atomAtomMapping.toString()) > 0) {
                            count++;
                        }
                    } catch (NumberFormatException ignore) {
                        // Non-numeric mapping markers do not count toward coverage.
                    }
                }
            }
        }
        return count;
    }

    private int getTotalNonHydrogenAtomCount(IAtomContainerSet mol) {
        int count = 0;
        List<IAtomContainer> allAtomContainers = getAllAtomContainers(mol);
        for (IAtomContainer ac : allAtomContainers) {
            IAtom[] atomArray = getAtomArray(ac);
            for (IAtom atom : atomArray) {
                if (!atom.getSymbol().equalsIgnoreCase("H")) {
                    count++;
                }
            }
        }
        return count;
    }

    private static final class MappingCoverage {

        private final int reactantAtoms;
        private final int productAtoms;
        private final int mappedReactantAtoms;
        private final int mappedProductAtoms;

        private MappingCoverage(int reactantAtoms, int productAtoms,
                int mappedReactantAtoms, int mappedProductAtoms) {
            this.reactantAtoms = reactantAtoms;
            this.productAtoms = productAtoms;
            this.mappedReactantAtoms = mappedReactantAtoms;
            this.mappedProductAtoms = mappedProductAtoms;
        }

        private boolean isComplete() {
            return mappedReactantAtoms == reactantAtoms
                    && mappedProductAtoms == productAtoms;
        }

        private boolean isBalancedMapped() {
            return mappedReactantAtoms == mappedProductAtoms;
        }

        private int getMappedAtoms() {
            return mappedReactantAtoms + mappedProductAtoms;
        }

        private int getUnmappedAtoms() {
            return (reactantAtoms - mappedReactantAtoms)
                    + (productAtoms - mappedProductAtoms);
        }
    }

    private static final class EvaluationCandidate {

        private final IMappingAlgorithm algorithm;
        private final Reactor reactor;
        private final IReaction mappedReaction;
        private final MappingCoverage coverage;
        private final String signature;
        private final QuickScore quickScore;

        private EvaluationCandidate(IMappingAlgorithm algorithm, Reactor reactor,
                IReaction mappedReaction, MappingCoverage coverage,
                String signature, QuickScore quickScore) {
            this.algorithm = algorithm;
            this.reactor = reactor;
            this.mappedReaction = mappedReaction;
            this.coverage = coverage;
            this.signature = signature;
            this.quickScore = quickScore;
        }
    }

    private static final class QuickScore {

        private final int bondChangeEstimate;
        private final int orderChangeEstimate;
        private final int carbonBondChangeEstimate;
        private final int fragmentPenalty;
        private final int unmappedBondPenalty;
        private final int mappedBondCount;

        private QuickScore(int bondChangeEstimate, int orderChangeEstimate,
                int carbonBondChangeEstimate, int fragmentPenalty,
                int unmappedBondPenalty,
                int mappedBondCount) {
            this.bondChangeEstimate = bondChangeEstimate;
            this.orderChangeEstimate = orderChangeEstimate;
            this.carbonBondChangeEstimate = carbonBondChangeEstimate;
            this.fragmentPenalty = fragmentPenalty;
            this.unmappedBondPenalty = unmappedBondPenalty;
            this.mappedBondCount = mappedBondCount;
        }

        private int totalScore() {
            if (bondChangeEstimate == Integer.MAX_VALUE) {
                return Integer.MAX_VALUE;
            }
            return bondChangeEstimate
                    + orderChangeEstimate
                    + carbonBondChangeEstimate
                    + fragmentPenalty;
        }

        private boolean isNear(QuickScore other) {
            if (other == null) {
                return false;
            }
            return Math.abs(totalScore() - other.totalScore()) <= 1
                    && Math.abs(bondChangeEstimate - other.bondChangeEstimate) <= 1
                    && Math.abs(orderChangeEstimate - other.orderChangeEstimate) <= 1;
        }

        private boolean hasEquivalentCoreScore(QuickScore other) {
            return other != null
                    && bondChangeEstimate == other.bondChangeEstimate
                    && orderChangeEstimate == other.orderChangeEstimate
                    && carbonBondChangeEstimate == other.carbonBondChangeEstimate
                    && fragmentPenalty == other.fragmentPenalty
                    && unmappedBondPenalty == other.unmappedBondPenalty;
        }
    }

    private static final class BondDescriptor {

        private final int order;
        private final boolean aromatic;
        private final boolean carbonOnly;

        private BondDescriptor(int order, boolean aromatic, boolean carbonOnly) {
            this.order = order;
            this.aromatic = aromatic;
            this.carbonOnly = carbonOnly;
        }

        private boolean sameType(BondDescriptor other) {
            return other != null
                    && order == other.order
                    && aromatic == other.aromatic;
        }
    }
}
