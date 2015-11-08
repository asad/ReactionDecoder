/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator;
import org.openscience.smsd.tools.BondEnergies;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mapping.CallableAtomMappingTool;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
/**
 
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionMechanismTool implements Serializable {


    private final static boolean DEBUG = false;
    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(ReactionMechanismTool.class);
    private static final long serialVersionUID = 07342630505L;
    private static final Logger LOG = Logger.getLogger(ReactionMechanismTool.class.getName());
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
    public ReactionMechanismTool(
            IReaction reaction,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D) throws Exception {
        this(reaction, forcedMapping, generate2D, generate3D, new StandardizeReaction());
    }

    /**
     *
     * @param reaction
     * @param forcedMapping
     * @param generate2D
     * @param generate3D
     * @param standardizer
     * @throws CDKException
     * @throws AssertionError
     * @throws Exception
     */
    public ReactionMechanismTool(IReaction reaction,
            boolean forcedMapping,
            boolean generate2D,
            boolean generate3D,
            IStandardizer standardizer) throws CDKException, AssertionError, Exception {
        this.allSolutions = Collections.synchronizedList(new ArrayList<MappingSolution>());
        this.selectedMapping = null;

        if (!isBalanced(reaction)) {
            Logger.getLogger(ReactionMechanismTool.class.getName()).
                    log(Level.WARNING, "Atoms not balanced in the input reaction: {0}; "
                            + "unbalanced reaction may result in erroneous bond change assumptions!", reaction.getID());
            Logger.getLogger(ReactionMechanismTool.class.getName()).
                    log(Level.WARNING, " ");
            if (!forcedMapping) {
                return;
            }
        }

        if (!forcedMapping && reaction.getFlag(CDKConstants.MAPPED)
                && AtomContainerSetManipulator.getAtomCount(reaction.getReactants())
                == reaction.getMappingCount()) {
            try {
                logger.info("Using user defined mappings!");
                /*
                 Set Atom IDs
                 */
                for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                    for (IAtom a : ac.atoms()) {
                        a.setID(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "");
                    }
                }
                for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                    for (IAtom a : ac.atoms()) {
                        a.setID(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "");
                    }
                }
                /*
                 Set mapped flags
                 */
                for (IMapping map : reaction.mappings()) {
                    if (map.getChemObject(0) != null && map.getChemObject(1) != null) {
                        map.getChemObject(0).setFlag(CDKConstants.MAPPED, true);
                        map.getChemObject(1).setFlag(CDKConstants.MAPPED, true);
                    }
                }
                boolean selected = isMappingSolutionAcceptable(null, IMappingAlgorithm.USER_DEFINED, reaction, generate2D, generate3D);
                logger.info("is solution: " + IMappingAlgorithm.USER_DEFINED + " selected: " + selected);
            } catch (Exception e) {
                String ls = System.getProperty("line.separator");
                throw new CDKException(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
            }
        } else {
            try {
                if (DEBUG) {
                    SmilesGenerator withAtomClasses = SmilesGenerator.unique().aromatic().withAtomClasses();
                    System.err.println("Input reaction mapped " + withAtomClasses.createReactionSMILES(reaction));
                }

                boolean onlyCoreMappingByMCS = true;
                CallableAtomMappingTool amt
                        = new CallableAtomMappingTool(reaction, standardizer, onlyCoreMappingByMCS);
                Map<IMappingAlgorithm, Reactor> solutions = amt.getSolutions();
                logger.info("!!!!Calculating Best Mapping Model!!!!");
                boolean selected;
                for (IMappingAlgorithm algorithm : solutions.keySet()) {
                    Reactor reactor = solutions.get(algorithm);

                    if (DEBUG) {
                        SmilesGenerator withAtomClasses = SmilesGenerator.unique().aromatic().withAtomClasses();
                        System.out.println("reaction mapped " + withAtomClasses.createReactionSMILES(reactor.getReactionWithAtomAtomMapping()));
                    }
                    int atomCountR = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getReactants());
                    int atomCountP = getNonHydrogenMappingAtomCount(reactor.getReactionWithAtomAtomMapping().getProducts());

                    if (atomCountR != atomCountP) {
                        logger.warn("ERROR in Mapping " + reactor.toString());
                        String newline = System.getProperty("line.separator");
                        System.err.println("Unmapped atoms present in this reaction" + "(" + algorithm + ") algorithm.");
//                        throw new AssertionError(newline + "Unmapped atoms present in the reaction mapped by AAM "
//                                + "(" + algorithm + ") algorithm." + newline);
                    }
                    selected = isMappingSolutionAcceptable(solutions.get(algorithm), algorithm, reactor.getReactionWithAtomAtomMapping(), generate2D, generate3D);
                    logger.info("is solution: " + algorithm + " selected: " + selected);
                }
                System.gc();
            } catch (Exception e) {
                String ls = System.getProperty("line.separator");
                throw new CDKException(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
            }
//            System.out.println(this.getMappingDescription());
        }
    }

    private boolean isBalanced(
            IReaction r) {
        
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
                    System.out.println("Q=mol " + SmilesGenerator.generic().create(q));
                } catch (CDKException ex) {
                    Logger.getLogger(ReactionMechanismTool.class.getName()).log(Level.SEVERE, null, ex);
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
                    System.out.println("T=mol " + SmilesGenerator.generic().create(t));
                } catch (CDKException ex) {
                    Logger.getLogger(ReactionMechanismTool.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        if (DEBUG) {
            System.out.println("atomUniqueCounter1 " + leftHandAtomCount);
            System.out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }
        
        if (leftHandAtomCount != rightHandAtomCount) {
            System.err.println();
            System.err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                    + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
            System.err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            System.err.println();
            System.err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                    + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
            System.err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            return false;
        }
        
        if (DEBUG) {
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
        }
        return atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet());
    }

    private synchronized boolean isMappingSolutionAcceptable(Reactor reactor, IMappingAlgorithm ma, IReaction reaction, boolean generate2D, boolean generate3D) throws Exception {
        
        boolean chosen = false;
        try {
            BondChangeCalculator bcc;
            int fragmentDeltaChanges;
            if (reactor == null && ma.equals(IMappingAlgorithm.USER_DEFINED)) {
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

                logger.info(
                        "Score: " + fragmentDeltaChanges + " : " + bondChange);
                logger.info(
                        ", Energy Barrier: " + bondBreakingEnergy);
                logger.info(
                        ", Energy Delta: " + bcc.getEnergyDelta());

                int localScore = bondChange + fragmentDeltaChanges;
                bcc.getReaction().setFlag(CDKConstants.MAPPED, true);

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
                logger.info("MA: " + ma.description());
                if (isChangeFeasible(mappingSolution)) {
                    chosen = true;
                    mappingSolution.setChosen(chosen);
                    this.selectedMapping = mappingSolution;
                }
                this.allSolutions.add(mappingSolution);
            }
        } catch (Exception e) {
            String ls = System.getProperty("line.separator");
            throw new CDKException(ls + "ERROR: Unable to calculate bond changes: " + e.getMessage());
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
                System.out.println("\n selectedMapping.getAlgorithmID().description() " + selectedMapping.getAlgorithmID().description());
                System.out.println(" selectedMapping.getTotalBondChanges() " + selectedMapping.getTotalBondChanges());
                System.out.println(" selectedMapping.getSmallestFragmentCount() " + selectedMapping.getSmallestFragmentCount());
                System.out.println(" selectedMapping.getBondEnergyChange() " + selectedMapping.getBondEnergySum());
                System.out.println(" selectedMapping.getTotalFragmentChanges() " + selectedMapping.getTotalFragmentChanges());
                System.out.println(" Total Carbon Bond Changes " + selectedMapping.getTotalCarbonBondChanges());
            }
            System.out.println("\n ms.getAlgorithmID().description() " + ms.getAlgorithmID().description());
            System.out.println(" ms.getTotalBondChanges() " + ms.getTotalBondChanges());
            System.out.println(" ms.getSmallestFragmentCount() " + ms.getSmallestFragmentCount());
            System.out.println(" ms.getBondEnergyChange() " + ms.getBondEnergySum());
            System.out.println(" ms.getTotalFragmentChanges() " + ms.getTotalFragmentChanges());
            System.out.println(" Total Carbon Bond Changes " + ms.getTotalCarbonBondChanges());
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
            logger.info("Condition Default " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition Default");
            }
            return true;
        } else if (ms.getBondEnergySum() == 0.
                && ms.getTotalFragmentChanges() == 0
                && ms.getTotalBondChanges() == 0
                && this.selectedMapping.getTotalStereoChanges() >= ms.getTotalStereoChanges()) {
            logger.info("Condition 1 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 1");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                || this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum())) {
            logger.info("Condition 2 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 2");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > 0
                && ms.getTotalFragmentChanges() > 0) {
            logger.info("Condition 3 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 3");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() >= ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() >= ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            /* Example reaction R05069*/
            logger.info("Condition 4 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 4");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            logger.info("Condition 5 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 5");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getSmallestFragmentCount() == ms.getSmallestFragmentCount()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            logger.info("Condition 6 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 6");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            logger.info("Condition 7 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 7");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalFragmentChanges() > ms.getTotalFragmentChanges()) {
            /* && this.selectedMapping.getBondEnergyChange() > ms.getBondEnergyChange()) {*/
            logger.info("Condition 8 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 8");
            }
            return true;
        } else if (this.selectedMapping.getTotalFragmentChanges() == ms.getTotalFragmentChanges()
                && this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()) {
            logger.info("Condition 9 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 9");
            }
            return true;
        } else if (this.selectedMapping.getBondEnergySum() == ms.getBondEnergySum()
                && this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalStereoChanges() > ms.getTotalStereoChanges()) {
            logger.info("Condition 10 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 10");
            }
            return true;
        } else if (this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()) {
            logger.info("Condition 11 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 11");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() < ms.getTotalBondChanges()
                && this.selectedMapping.getBondEnergySum() < ms.getBondEnergySum()
                && this.selectedMapping.getTotalCarbonBondChanges() > 0
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: R00652 N Vs O exchange*/
            logger.info("Condition 12 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 12");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() > ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() > ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like: Rhea 22881 N Vs O exchange*/
            logger.info("Condition 13 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 13");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getSmallestFragmentCount() > ms.getSmallestFragmentCount()) {
            /*This condition is for reactions like:
            CC1=C2NC(C(CC(O)=O)C2(C)CCC(O)=O)C2(C)N=C(C(CCC(O)=O)C2(C)CC(O)=O)C(C)=C2N=C(C=C3N=C1C(CCC(O)=O)C3(C)C)C(CCC(O)=O)C2(C)CC(O)=O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O.NC(CCC(N)=O)C(O)=O.NC(CCC(N)=O)C(O)=O.O[H].O[H]>>CC1=C2NC(C(CC(O)=O)C2(C)CCC(O)=O)C2(C)N=C(C(CCC(O)=O)C2(C)CC(N)=O)C(C)=C2N=C(C=C3N=C1C(CCC(O)=O)C3(C)C)C(CCC(O)=O)C2(C)CC(N)=O.Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(O)=O)C(O)C1O.[H]OP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1O)n1cnc2c(N)ncnc12.NC(CCC(O)=O)C(O)=O.NC(CCC(O)=O)C(O)=O.[H]OP(O)(O)=O.OP(O)(O)=O
            **/
            logger.info("Condition 14 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 14");
            }
            return true;
        } else if (this.selectedMapping.getTotalBondChanges() == ms.getTotalBondChanges()
                && this.selectedMapping.getTotalCarbonBondChanges() == ms.getTotalCarbonBondChanges()
                && this.selectedMapping.getBondEnergySum() > ms.getBondEnergySum()) {
            /*This condition is for reactions like:
            R05421 (O-P over O-C)
            **/
            logger.info("Condition 15 " + ms.getAlgorithmID().description());
            if (DEBUG) {
                System.out.println("CASE: Condition 15");
            }
            return true;
        }

        if (DEBUG) {
            System.out.println("CASE: FAILED");
        }
        return false;
    }

    private synchronized double getTotalBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        for (IFeature key : fingerprint.getFeatures()) {
            double val = key.getWeight();
            if (val > 0.) {//&& !key.contains("H")
                total += val;
            }
        }
        return total;
    }

    private synchronized int getTotalCarbonBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        for (IFeature key : fingerprint.getFeatures()) {
            if (key.getPattern().contains("C-C")
                    || key.getPattern().contains("C=C")
                    || key.getPattern().contains("C#C")
                    || key.getPattern().contains("C%C")
                    || key.getPattern().contains("C@C")) {
                double val = key.getWeight();
                if (val > 0.) {//&& !key.contains("H")
                    total += val;
                }
            }
        }
        return (int) total;
    }

    private synchronized int getTotalBondChangeEnergy(IPatternFingerprinter fingerprint, boolean skipHydrogen) {
        int total = 0;
        try {
            BondEnergies be = BondEnergies.getInstance();
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
                        int energy = be.getEnergies(temp[0], temp[1], IBond.Order.SINGLE);
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
                        int energy = be.getEnergies(temp[0], temp[1], IBond.Order.DOUBLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("#")) {
                        String[] temp = key.split("#");
                        int energy = be.getEnergies(temp[0], temp[1], IBond.Order.TRIPLE);
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
                        int energy = be.getEnergies(temp[0], temp[1], IBond.Order.QUADRUPLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    }
                }
            }
        } catch (CDKException ex) {
            Logger.getLogger(ReactionMechanismTool.class.getName()).log(Level.SEVERE, null, ex);
        }
        return Math.abs(total);
    }

    /**
     *
     * @return
     */
    public String getMappingDescription() {
        return this.selectedMapping.toString();
    }

    public MappingSolution getSelectedSolution() {
        return this.selectedMapping;
    }

    public Collection<MappingSolution> getAllSolutions() {
        return Collections.unmodifiableCollection(this.allSolutions);
    }

    private int getNonHydrogenMappingAtomCount(IAtomContainerSet mol) {
        int count = Integer.MIN_VALUE;
        List<IAtomContainer> allAtomContainers = AtomContainerSetManipulator.getAllAtomContainers(mol);
        for (IAtomContainer ac : allAtomContainers) {
            IAtom[] atomArray = ExtAtomContainerManipulator.getAtomArray(ac);
            for (IAtom atom : atomArray) {
                if (atom.getSymbol().equalsIgnoreCase("H")) {
                    continue;
                }
                if (atom.getID() != null && Integer.parseInt(atom.getID()) > count) {
                    count = Integer.parseInt(atom.getID());
                }
            }
        }
        return count;
    }
}
