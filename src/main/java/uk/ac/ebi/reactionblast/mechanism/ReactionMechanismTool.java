/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mechanism;

import java.io.Serializable;
import static java.lang.Integer.MIN_VALUE;
import static java.lang.System.gc;
import static java.lang.System.out;
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
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAllAtomContainers;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAtomCount;
import org.openscience.smsd.tools.BondEnergies;
import static org.openscience.smsd.tools.BondEnergies.getInstance;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mapping.CallableAtomMappingTool;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.USER_DEFINED;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import static java.lang.Integer.parseInt;
import static java.lang.Math.abs;
import static java.lang.System.getProperty;
import static java.util.Collections.synchronizedList;

import org.openscience.cdk.smiles.SmiFlavor;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getAtomArray;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionMechanismTool implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static boolean DEBUG = false;
    private final static ILoggingTool LOGGER
            = createLoggingTool(ReactionMechanismTool.class);
    private static final long serialVersionUID = 07342630505L;
    private MappingSolution selectedMapping;
    private Collection<MappingSolution> allSolutions;

    /**
     *
     * @param reaction
     * @param forcedMapping force re-mapping of the reactions
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     * @throws Exception
     */
    public ReactionMechanismTool(IReaction reaction, boolean forcedMapping, boolean generate2D, boolean generate3D) throws Exception {
        this(reaction, forcedMapping, generate2D, generate3D, new StandardizeReaction());
    }

    /**
     *
     * @param reaction
     * @param forcedMapping
     * @param generate2D deduce stereo on 2D
     * @param generate3D deduce stereo on 3D
     * @param standardizer
     * @throws CDKException
     * @throws AssertionError
     * @throws Exception
     */
    public ReactionMechanismTool(IReaction reaction, boolean forcedMapping,
            boolean generate2D, boolean generate3D, IStandardizer standardizer) throws CDKException, AssertionError, Exception {
        this.allSolutions = synchronizedList(new ArrayList<>());
        this.selectedMapping = null;

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
                boolean selected = isMappingSolutionAcceptable(null, USER_DEFINED, reaction, generate2D, generate3D);
                LOGGER.info("is solution: " + USER_DEFINED + " selected: " + selected);
            } catch (Exception e) {
                String ls = getProperty("line.separator");
                throw new CDKException(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
            }
        } else {
            try {
                if (DEBUG) {
                    SmilesGenerator withAtomClasses = new SmilesGenerator(
                            SmiFlavor.Unique
                            | SmiFlavor.Stereo
                            | SmiFlavor.AtomAtomMap);
                    LOGGER.debug("Input reaction mapped " + withAtomClasses.create(reaction));
                }

                boolean onlyCoreMappingByMCS = true;
                CallableAtomMappingTool amt
                        = new CallableAtomMappingTool(reaction, standardizer, onlyCoreMappingByMCS);
                Map<IMappingAlgorithm, Reactor> solutions = amt.getSolutions();

                LOGGER.info("!!!!Calculating Best Mapping Model!!!!");
                boolean selected;
                for (IMappingAlgorithm algorithm : solutions.keySet()) {

                    Reactor reactor = solutions.get(algorithm);

                    if (reactor == null) {
                        System.out.println("Reactor is NULL");
                        return;
                    }

                    if (DEBUG) {

                        SmilesGenerator withAtomClasses = new SmilesGenerator(
                                SmiFlavor.Unique
                                | SmiFlavor.Stereo
                                | SmiFlavor.AtomAtomMap);
                        out.println("reaction mapped " + withAtomClasses.create(reactor.getReactionWithAtomAtomMapping()));
                    }

                    int atomCountR = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getReactants());
                    int atomCountP = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getProducts());

                    if (atomCountR != atomCountP) {
                        LOGGER.warn("ERROR in Mapping " + reactor.toString());
                        LOGGER.warn("Unmapped atoms present in this reaction" + "(" + algorithm + ") algorithm.");
//                        throw new AssertionError(newline + "Unmapped atoms present in the reaction mapped by AAM "
//                                + "(" + algorithm + ") algorithm." + newline);
                    }
                    selected = isMappingSolutionAcceptable(solutions.get(algorithm), algorithm, reactor.getReactionWithAtomAtomMapping(), generate2D, generate3D);
                    LOGGER.info("is solution: " + algorithm + " selected: " + selected);
                }
                gc();
            } catch (Exception e) {
                String ls = getProperty("line.separator");
                throw new Exception(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
            }
//            System.out.println(this.getMappingDescription());
        }
    }

    private boolean isBalanced(IReaction r) {

        Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
        Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();

        int leftHandAtomCount = 0;
        for (IAtomContainer q : r.getReactants().atomContainers()) {
            for (IAtom a : q.atoms()) {
                if (a.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                    atomUniqueCounter1.put(a.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                    atomUniqueCounter1.put(a.getSymbol(), counter);
                }
                leftHandAtomCount++;
            }
            if (DEBUG) {
                try {
                    out.println("Q=mol " + generic().create(q));
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }

        int rightHandAtomCount = 0;
        for (IAtomContainer t : r.getProducts().atomContainers()) {
            for (IAtom b : t.atoms()) {
                if (b.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                    atomUniqueCounter2.put(b.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                    atomUniqueCounter2.put(b.getSymbol(), counter);
                }
                rightHandAtomCount++;
            }
            if (DEBUG) {
                try {
                    out.println("T=mol " + generic().create(t));
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + leftHandAtomCount);
            out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }

        if (leftHandAtomCount != rightHandAtomCount) {
            LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                    + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
            LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                    + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
            LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            return false;
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            out.println("atomUniqueCounter2 " + atomUniqueCounter2);
        }
        return atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet());
    }

    private synchronized boolean isMappingSolutionAcceptable(Reactor reactor,
            IMappingAlgorithm ma,
            IReaction reaction,
            boolean generate2D,
            boolean generate3D) throws Exception {

        boolean chosen = false;
        try {
            BondChangeCalculator bcc;
            int fragmentDeltaChanges;
            if (reactor == null && ma.equals(USER_DEFINED)) {
                bcc = new BondChangeCalculator(reaction, generate2D, generate3D);
                fragmentDeltaChanges = 0;
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
                        bcc.getEnergyDelta());

                chosen = true;
                mappingSolution.setChosen(chosen);
                this.selectedMapping = mappingSolution;
                this.allSolutions.add(mappingSolution);
            } else {
                if (reactor == null) {
                    throw new CDKException("Reactor is NULL");
                }

                bcc = new BondChangeCalculator(reactor.getReactionWithAtomAtomMapping(), generate2D, generate3D);

                fragmentDeltaChanges = reactor.getDelta();

                int bondCleavedFormed = (int) getTotalBondChange(bcc.getFormedCleavedWFingerprint());
                int bondChange = bondCleavedFormed;
                bondChange += getTotalBondChange(bcc.getOrderChangesWFingerprint());
                int stereoChanges = (int) getTotalBondChange(bcc.getStereoChangesWFingerprint());
                boolean skipHydrogenRealtedBondChanges = true;
                int bondBreakingEnergy = getTotalBondChangeEnergy(bcc.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);
                int totalSmallestFragmentCount = bcc.getTotalSmallestFragmentSize();
                int totalCarbonBondChanges = getTotalCarbonBondChange(bcc.getFormedCleavedWFingerprint());

                LOGGER.info(
                        "Score: " + fragmentDeltaChanges + " : " + bondChange);
                LOGGER.info(
                        ", Energy Barrier: " + bondBreakingEnergy);
                LOGGER.info(
                        ", Energy Delta: " + bcc.getEnergyDelta());

                int localScore = bondChange + fragmentDeltaChanges;
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
                if (isChangeFeasible(mappingSolution)) {
                    chosen = true;
                    mappingSolution.setChosen(chosen);
                    this.selectedMapping = mappingSolution;
                }
                this.allSolutions.add(mappingSolution);
            }
        } catch (Exception e) {
            String ls = getProperty("line.separator");
            throw new Exception(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
        }
        return chosen;
    }

    /*
    * if bond changes are lesser than stored bond changes then update the flag or if stereo changes are lesser than
    * stores stereo changes
     */
    private synchronized boolean isChangeFeasible(MappingSolution ms) {

        /*
        * This condition is valuble to trace graph isomorphism as only min algorithm checks this change. the idea is to
        * assume a change if rest of the algorithm detects no change.
        *
        * TODO: check what is the impact if this logic if there are only stereo changes in a reaction.
         */
        if (DEBUG) {

            if (this.selectedMapping != null) {
                out.println(NEW_LINE + " selectedMapping.getAlgorithmID().description() " + selectedMapping.getAlgorithmID().description());
                out.println(" selectedMapping.getTotalBondChanges() " + selectedMapping.getTotalBondChanges());
                out.println(" selectedMapping.getSmallestFragmentCount() " + selectedMapping.getSmallestFragmentCount());
                out.println(" selectedMapping.getBondEnergyChange() " + selectedMapping.getBondEnergySum());
                out.println(" selectedMapping.getTotalFragmentChanges() " + selectedMapping.getTotalFragmentChanges());
                out.println(" Total Carbon Bond Changes " + selectedMapping.getTotalCarbonBondChanges());
            }
            out.println(NEW_LINE + " ms.getAlgorithmID().description() " + ms.getAlgorithmID().description());
            out.println(" ms.getTotalBondChanges() " + ms.getTotalBondChanges());
            out.println(" ms.getSmallestFragmentCount() " + ms.getSmallestFragmentCount());
            out.println(" ms.getBondEnergyChange() " + ms.getBondEnergySum());
            out.println(" ms.getTotalFragmentChanges() " + ms.getTotalFragmentChanges());
            out.println(" Total Carbon Bond Changes " + ms.getTotalCarbonBondChanges());
        }

        /*
        This is to skip reaction where the no change is detected.
        Example: R02996
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
            if (DEBUG) {
                out.println("CASE: Condition Default");
            }
            return true;
        } else if (ms.getBondEnergySum() == 0.
                && ms.getTotalFragmentChanges() == 0
                && ms.getTotalBondChanges() == 0
                && this.selectedMapping.getTotalStereoChanges() >= ms.getTotalStereoChanges()) {
            LOGGER.info("Condition 1 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 1");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                || this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum())) {
            LOGGER.info("Condition 2 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 2");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > 0
                && ms.getTotalFragmentChanges() > 0) {
            LOGGER.info("Condition 3 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 3");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() >= ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() >= ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() >= ms.getTotalCarbonBondChanges()) {
            /* Example reaction R05069*/
            LOGGER.info("Condition 4 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 4");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            LOGGER.info("Condition 5 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 5");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() == ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() >= ms.getTotalCarbonBondChanges()) {
            LOGGER.info("Condition 6 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 6");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            LOGGER.info("Condition 7 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 7");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()) {
            /* && this.selectedMapping.getBondEnergyChange() > ms.getBondEnergyChange()) {*/
            LOGGER.info("Condition 8 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 8");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()) {
            LOGGER.info("Condition 9 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 9");
            }
            return true;
        } else if (this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalStereoChanges() > ms.getTotalStereoChanges()) {
            LOGGER.info("Condition 10 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 10");
            }
            return true;
        } else if (this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()) {
            LOGGER.info("Condition 11 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 11");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() < ms.getTotalBondChanges()
                && this.selectedMapping.getBondEnergySum() < ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: R00652 N Vs O exchange*/
            LOGGER.info("Condition 12 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 12");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: Rhea 22881 N Vs O exchange*/
            LOGGER.info("Condition 13 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 13");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like:
            CC1=C2NC(C(CC(O)=O)C2(C)CCC(O)=O)C2(C)N=C(C(CCC(O)=O)C2(C)CC(O)=O)C(C)=C2N=C(C=C3N=C1C(CCC(O)=O)C3(C)C)C(CCC(O)=O)C2(C)CC(O)=O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O.NC(CCC(N)=O)C(O)=O.NC(CCC(N)=O)C(O)=O.O[H].O[H]>>CC1=C2NC(C(CC(O)=O)C2(C)CCC(O)=O)C2(C)N=C(C(CCC(O)=O)C2(C)CC(N)=O)C(C)=C2N=C(C=C3N=C1C(CCC(O)=O)C3(C)C)C(CCC(O)=O)C2(C)CC(N)=O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(O)=O)C(O)C1O.[H]OP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1O)n1cnc2c(N)ncnc12.NC(CCC(O)=O)C(O)=O.NC(CCC(O)=O)C(O)=O.[H]OP(O)(O)=O.OP(O)(O)=O
            **/
            LOGGER.info("Condition 14 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 14");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            /*This condition is for reactions like:
            R05421 (O-P over O-C)
            **/
            LOGGER.info("Condition 15 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                out.println("CASE: Condition 15");
            }
            return true;
        }

        if (DEBUG) {
            out.println("CASE: FAILED");
        }
        return false;
    }

    private synchronized double getTotalBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return total;
    }

    private synchronized int getTotalCarbonBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().filter((key) -> (key.getPattern().contains("C-C")
                || key.getPattern().contains("C=C")
                || key.getPattern().contains("C#C")
                || key.getPattern().contains("C%C")
                || key.getPattern().contains("C@C"))).map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return (int) total;
    }

    private synchronized int getTotalBondChangeEnergy(IPatternFingerprinter fingerprint, boolean skipHydrogen) {
        int total = 0;
        try {
            BondEnergies be = getInstance();
            for (IFeature feature : fingerprint.getFeatures()) {
                double val = feature.getWeight();
                String key = feature.getPattern();
                if (val > 0) {
//                    System.out.println("BOND BROKEN/FORMED: " + key + " : " + val);
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
                        int energy = be.getEnergies(temp[0], temp[1], SINGLE);
                        if (energy > 0) {
                            /*
                            * Ring energy correction factor example:R01081
                             */
                            if (key.contains("%")) {
                                total += val * (energy - 5.0);

                            } else {
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
                        String[] temp = key.split("$");
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
