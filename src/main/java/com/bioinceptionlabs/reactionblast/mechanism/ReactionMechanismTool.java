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
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAllAtomContainers;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAtomCount;
import org.openscience.smsd.BondEnergies;
import static org.openscience.smsd.BondEnergies.getInstance;
import com.bioinceptionlabs.reactionblast.fingerprints.IFeature;
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

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ReactionMechanismTool implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static ILoggingTool LOGGER
            = createLoggingTool(ReactionMechanismTool.class);
    private static final long serialVersionUID = 07342630505L;
    private MappingSolution selectedMapping;
    private Collection<MappingSolution> allSolutions;
    private final boolean accept_no_change;

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

                LOGGER.debug("!!!!Calculating Best Mapping Model!!!!");
                boolean selected;
                for (IMappingAlgorithm algorithm : solutions.keySet()) {

                    Reactor reactor = solutions.get(algorithm);

                    if (reactor == null) {
                        LOGGER.warn("Reactor is NULL");
                        return;
                    }

                    int atomCountR = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getReactants());
                    int atomCountP = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getProducts());

                    if (atomCountR != atomCountP) {
                        //LOGGER.warn("ERROR in Mapping - Unmapped atoms present in the reaction: "
                        //        + NEW_LINE + reactor.toString());
                        LOGGER.warn("Unmapped atoms present in this reaction" + "(" + algorithm + ") algorithm.");
//                        throw new AssertionError(newline + "Unmapped atoms present in the reaction mapped by AAM "
//                                + "(" + algorithm + ") algorithm." + newline);
                    }
                    LOGGER.debug("===isMappingSolutionAcceptable===");
                    selected = isMappingSolutionAcceptable(solutions.get(algorithm),
                            algorithm,
                            reactor.getReactionWithAtomAtomMapping(),
                            generate2D,
                            generate3D);
                    LOGGER.debug("is solution: " + algorithm + " selected: " + selected);

                    // Early exit: if this solution has minimal bond changes (≤2)
                    // and zero fragment changes, it's likely optimal — skip remaining algorithms
                    if (selected && this.selectedMapping != null
                            && this.selectedMapping.getTotalBondChanges() <= 2
                            && this.selectedMapping.getTotalFragmentChanges() == 0) {
                        LOGGER.debug("Early exit: optimal mapping found by " + algorithm);
                        break;
                    }
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
            LOGGER.warn("Number of atom(s) on the Left side "
                    + reactantAtoms.values().stream().mapToInt(Integer::intValue).sum()
                    + " =/= Number of atom(s) on the Right side "
                    + productAtoms.values().stream().mapToInt(Integer::intValue).sum());
            LOGGER.warn(reactantAtoms + " =/= " + productAtoms);
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

    private int getNonHydrogenMappingAtomCount(IAtomContainerSet mol) {
        int count = MIN_VALUE;
        List<IAtomContainer> allAtomContainers = getAllAtomContainers(mol);
        for (IAtomContainer ac : allAtomContainers) {
            IAtom[] atomArray = getAtomArray(ac);
            for (IAtom atom : atomArray) {
                if (atom.getSymbol().equalsIgnoreCase("H")) {
                    continue;
                }
                if (atom.getID() != null && parseInt(atom.getID()) > count) {
                    count = parseInt(atom.getID());
                }
            }
        }
        return count;
    }
}
