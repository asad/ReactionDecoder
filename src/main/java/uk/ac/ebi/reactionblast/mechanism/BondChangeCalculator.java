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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.abs;
import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.synchronizedList;
import static java.util.Collections.synchronizedMap;
import static java.util.Collections.unmodifiableCollection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import static org.openscience.cdk.aromaticity.Kekulization.kekulize;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import static org.openscience.cdk.graph.ConnectivityChecker.partitionIntoMolecules;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import static org.openscience.cdk.tools.CDKHydrogenAdder.getInstance;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getRelevantAtomContainer;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getRelevantAtomContainer;
import org.openscience.smsd.tools.BondEnergies;
import static org.openscience.smsd.tools.BondEnergies.getInstance;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomStereoChangeInformation;
import uk.ac.ebi.reactionblast.mechanism.helper.BondChange;
import uk.ac.ebi.reactionblast.mechanism.helper.MoleculeMoleculePair;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionCenterFragment;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCircularSMILES;
import uk.ac.ebi.reactionblast.mechanism.interfaces.AbstractChangeCalculator;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.BOND_CHANGE_INFORMATION;
import uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct.PRODUCT;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct.REACTANT;
import uk.ac.ebi.reactionblast.mechanism.interfaces.IChangeCalculator;
import static uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool.deepClone;

/**
 * This class marks the bond changes
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondChangeCalculator extends AbstractChangeCalculator implements IChangeCalculator {

    private final boolean DEBUG = false;
    private static final long serialVersionUID = 98698690880809981L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(BondChangeCalculator.class);
    private final BondChangeAnnotator bondChangeAnnotator;
    private final IPatternFingerprinter formedCleavedWFingerprint;
    private final IPatternFingerprinter orderChangesWFingerprint;
    private final IPatternFingerprinter stereoChangesWFingerprint;
    private final IPatternFingerprinter reactionCenterWFingerprint;
    private final Map<Integer, IPatternFingerprinter> reactionCenterFormedCleavedFingerprint;
    private final Map<Integer, IPatternFingerprinter> reactionCenterOrderChangeFingerprint;
    private final Map<Integer, IPatternFingerprinter> reactionCenterStereoChangeFingerprint;
    private final Map<IBond, String> bondFormedMap;
    private final Map<IBond, String> bondCleavedMap;
    private final Map<IBond, String> bondOrderRMap;
    private final Map<IBond, String> bondOrderPMap;
    private final Map<IAtom, String> AtomStereoRMap;
    private final Map<IAtom, String> AtomStereoPMap;
    private final List<ReactionCenterFragment> reactionCenterFragmentList;
    private final Set<MoleculeMoleculePair> reactionMoleculeMoleculePairList;
    private final IReaction mappedReaction;
    private int energySum;
    private int energyDelta;
    private int totalSmallestFragmentSize;

    /**
     *
     * @param reaction
     * @param generate2D
     * @param generate3D
     * @throws Exception
     */
    public BondChangeCalculator(IReaction reaction, boolean generate2D, boolean generate3D) throws Exception {
        if (DEBUG) {
            System.out.println("Bond Change Calculator START");
        }
        int rEnergy = 0;
        int pEnergy = 0;

        this.energySum = 0;
        this.energyDelta = 0;
        this.totalSmallestFragmentSize = 0;
        this.mappedReaction = reaction;
        if (DEBUG) {
            System.out.println("Bond Change Annotator START");
        }
        this.bondChangeAnnotator = new BondChangeAnnotator(this.mappedReaction, true, generate2D, generate3D);
        if (DEBUG) {
            System.out.println("Bond Change Annotator END");
        }

        BondEnergies be = getInstance();

        this.formedCleavedWFingerprint = new PatternFingerprinter();
        this.formedCleavedWFingerprint.setFingerprintID(reaction.getID() + ":" + "Bond Cleaved and Formed");
        this.orderChangesWFingerprint = new PatternFingerprinter();
        this.orderChangesWFingerprint.setFingerprintID(reaction.getID() + ":" + "Bond Order Change");
        this.stereoChangesWFingerprint = new PatternFingerprinter();
        this.stereoChangesWFingerprint.setFingerprintID(reaction.getID() + ":" + "Bond Stereo Change");
        this.reactionCenterWFingerprint = new PatternFingerprinter();
        this.reactionCenterWFingerprint.setFingerprintID(reaction.getID() + ":" + "Reaction Center");

        this.reactionCenterFormedCleavedFingerprint = new HashMap<>();
        this.reactionCenterOrderChangeFingerprint = new HashMap<>();
        this.reactionCenterStereoChangeFingerprint = new HashMap<>();

        this.reactionMoleculeMoleculePairList = new LinkedHashSet<>();

        this.bondFormedMap = synchronizedMap(new HashMap<>());
        this.bondCleavedMap = synchronizedMap(new HashMap<>());
        this.bondOrderRMap = synchronizedMap(new HashMap<>());
        this.bondOrderPMap = synchronizedMap(new HashMap<>());
        this.AtomStereoRMap = synchronizedMap(new HashMap<>());
        this.AtomStereoPMap = synchronizedMap(new HashMap<>());

        this.reactionCenterFragmentList = synchronizedList(new ArrayList<>());

        if (DEBUG) {
            System.out.println("Loop for stereo conformation changes count: " + bondChangeAnnotator.getConformationChangeList().size());
        }
        /*
         * Loop for stereo changes
         */
        for (AtomStereoChangeInformation atomConformation : bondChangeAnnotator.getConformationChangeList()) {

            /*
             * Stereo changes are marked only once in the Fingerprint
             */
            if (atomConformation.getReactantAtom() != null && atomConformation.getProductAtom() != null) {
                String keyE = atomConformation.getReactantAtom().getSymbol().concat("(E/Z)");
                IFeature eductFeature = new Feature(keyE, 1.0);
                stereoChangesWFingerprint.add(eductFeature);
            }

            /*
             * Stereo changes are marked on reactant and product for reaction center identification
             */
            if (atomConformation.getReactantAtom() != null) {
                atomConformation.getReactantAtom().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                AtomStereoRMap.put(atomConformation.getReactantAtom(), getMoleculeID(atomConformation.getReactantAtom(), reaction.getReactants()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomR1 = atomConformation.getReactantAtom();
                IAtomContainer moleculeR = getAtomContainer(atomConformation.getReactantAtom(), reaction.getReactants());

                if (moleculeR.getAtomCount() > 1) {
                    if (!atomR1.getSymbol().equals("H")) {
                        if (DEBUG) {
                            System.out.println("Educt CircularFingerprints START");
                        }
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeR, atomR1, REACTANT));
                        setCircularFingerprints(reaction.getID(), moleculeR, atomR1, reactionCenterStereoChangeFingerprint);
                        if (DEBUG) {
                            System.out.println("Educt CircularFingerprints END");
                        }
                    }
                }
            }
            if (atomConformation.getProductAtom() != null) {
                atomConformation.getProductAtom().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                AtomStereoPMap.put(atomConformation.getProductAtom(), getMoleculeID(atomConformation.getProductAtom(), reaction.getProducts()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomP1 = atomConformation.getProductAtom();
                IAtomContainer moleculeP = getAtomContainer(atomConformation.getProductAtom(), reaction.getProducts());

                if (moleculeP.getAtomCount() > 1) {
                    if (!atomP1.getSymbol().equals("H")) {
                        if (DEBUG) {
                            System.out.println("Product CircularFingerprints START");
                        }
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, PRODUCT));
                        setCircularFingerprints(reaction.getID(), moleculeP, atomP1, reactionCenterStereoChangeFingerprint);
                        if (DEBUG) {
                            System.out.println("Product CircularFingerprints END");
                        }
                    }
                }
            }
        }

        if (DEBUG) {
            System.out.println("Loop for stereo changes count: " + bondChangeAnnotator.getStereoChangeList().size());
        }

        /*
         * Loop for stereo changes
         */
        for (AtomStereoChangeInformation atomStereo : bondChangeAnnotator.getStereoChangeList()) {

            /*
             * Stereo changes are marked only once in the Fingerprint
             */
            if (atomStereo.getReactantAtom() != null && atomStereo.getProductAtom() != null) {
                String key = atomStereo.getReactantAtom().getSymbol().concat("(R/S)");
                IFeature eductFeature = new Feature(key, 1.0);
                stereoChangesWFingerprint.add(eductFeature);
            }

            /*
             * Stereo changes are marked on reactant and product for reaction center identification
             */
            if (atomStereo.getReactantAtom() != null) {
                atomStereo.getReactantAtom().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                AtomStereoRMap.put(atomStereo.getReactantAtom(), getMoleculeID(atomStereo.getReactantAtom(), reaction.getReactants()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomR1 = atomStereo.getReactantAtom();
                IAtomContainer moleculeR = getAtomContainer(atomStereo.getReactantAtom(), reaction.getReactants());

                if (moleculeR.getAtomCount() > 1) {

                    if (!atomR1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeR, atomR1, REACTANT));
                        setCircularFingerprints(reaction.getID(), moleculeR, atomR1, reactionCenterStereoChangeFingerprint);
                    }
                }
            }
            if (atomStereo.getProductAtom() != null) {
                atomStereo.getProductAtom().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                AtomStereoPMap.put(atomStereo.getProductAtom(), getMoleculeID(atomStereo.getProductAtom(), reaction.getProducts()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomP1 = atomStereo.getProductAtom();
                IAtomContainer moleculeP = getAtomContainer(atomStereo.getProductAtom(), reaction.getProducts());

                if (moleculeP.getAtomCount() > 1) {

                    if (!atomP1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, PRODUCT));
                        setCircularFingerprints(reaction.getID(), moleculeP, atomP1, reactionCenterStereoChangeFingerprint);
                    }
                }
            }
        }

        /*
         * Loop over atom order and generate unique list to atoms
         */
        Set<IAtom> reactantAtoms = new HashSet<>();
        Set<IAtom> productAtoms = new HashSet<>();

        if (DEBUG) {
            System.out.println("Bond Change List: " + bondChangeAnnotator.getBondChangeList().size());
        }

        for (BondChange bcinfo : bondChangeAnnotator.getBondChangeList()) {
            IBond bondR = bcinfo.getReactantBond();
            IBond bondP = bcinfo.getProductBond();

            // Mark Bond Order Changes
            if (bondR != null && bondP != null
                    && bondP.getProperties().get(BOND_CHANGE_INFORMATION).
                            equals(BOND_ORDER)
                    && bondR.getProperties().get(BOND_CHANGE_INFORMATION).
                            equals(BOND_ORDER)) {

                bondOrderRMap.put(bondR, getMoleculeID(bondR, reaction.getReactants()));
                bondR.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER);
                bondR.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER);

                bondOrderPMap.put(bondP, getMoleculeID(bondP, reaction.getProducts()));
                bondP.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER);
                bondP.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER);

                reactantAtoms.add(bondR.getAtom(0));
                reactantAtoms.add(bondR.getAtom(1));

                productAtoms.add(bondP.getAtom(0));
                productAtoms.add(bondP.getAtom(1));
                orderChangesWFingerprint.add(new Feature(getCanonisedBondChangePattern(bondR, bondP), 1.0));
            }
        }

        if (DEBUG) {
            System.out.println("Bond Order changes");
        }

        /*
         * Store changes in the bond order
         */
        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet products = reaction.getProducts();

        for (IAtom atom : reactantAtoms) {
            IAtomContainer relevantAtomContainer = getRelevantAtomContainer(reactants, atom);
            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, REACTANT));
            setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterOrderChangeFingerprint);
        }

        for (IAtom atom : productAtoms) {
            IAtomContainer relevantAtomContainer = getRelevantAtomContainer(products, atom);
            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, PRODUCT));
            setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterOrderChangeFingerprint);
        }

        if (DEBUG) {
            System.out.println("Bond formed, cleaved changes");
        }


        /*
         * Loop for formed, cleaved changes
         */
        for (BondChange bcinfo : bondChangeAnnotator.getBondChangeList()) {
            IBond bondR = bcinfo.getReactantBond();
            IBond bondP = bcinfo.getProductBond();

            if (DEBUG) {
                System.out.println("Bond formed, cleaved changes 1 ");
            }

            //Mark Formed Bonds in the Product
            if (bondP != null && (bondP.getProperties().get(BOND_CHANGE_INFORMATION).
                    equals(BOND_FORMED)
                    || bondP.getProperties().get(BOND_CHANGE_INFORMATION).
                            equals(PSEUDO_BOND))) {

                if (DEBUG) {
                    System.out.println("Bond formed, cleaved changes 1 - 1");
                }

                if (!bondP.getAtom(0).getSymbol().equals("PsH")
                        && !bondP.getAtom(1).getSymbol().equals("PsH")) {
                    this.energySum += be.getEnergies(bondP);
                    pEnergy += be.getEnergies(bondP);
                    bondFormedMap.put(bondP, getMoleculeID(bondP, reaction.getProducts()));
                    bondP.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, BOND_FORMED);
                    bondP.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, BOND_FORMED);

                    if (DEBUG) {
                        System.out.println("Bond formed, cleaved changes 1 - 1 - 1");
                    }

                    /*
                     * Update Reaction center FP
                     */
                    IAtomContainer moleculeP = getAtomContainer(bondP, reaction.getProducts());

                    if (DEBUG) {
                        System.out.println("Bond formed, cleaved changes FP");
                    }

                    if (moleculeP != null && moleculeP.getAtomCount() > 1) {

                        if (DEBUG) {
                            System.out.println("Bond formed, cleaved changes FP IN");
                        }
                        /*
                         * Mark reaction centers
                         */
                        IAtom atomP1 = bondP.getAtom(0);
                        IAtom atomP2 = bondP.getAtom(1);
                        if (DEBUG) {
                            System.out.println("Bond formed, cleaved changes 1 - 1 - 1 FP");
                        }
                        if (!atomP1.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, PRODUCT));
                            setCircularFingerprints(reaction.getID(), moleculeP, atomP1, reactionCenterFormedCleavedFingerprint);
                        }
                        if (!atomP2.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP2, PRODUCT));
                            setCircularFingerprints(reaction.getID(), moleculeP, atomP2, reactionCenterFormedCleavedFingerprint);
                        }

                        if (DEBUG) {
                            System.out.println("Bond formed, cleaved changes 1 - 1 - 1 FP Done");
                        }

                        IAtomContainer product = getAtomContainer(bondP, reaction.getProducts());
                        IAtomContainer cloneProduct = product.getBuilder().newInstance(IAtomContainer.class, product);
                        int chippedBondIndex = product.indexOf(bondP);
                        totalSmallestFragmentSize += chipTheBondCountSmallestFragmentSize(cloneProduct, chippedBondIndex);
                        formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bondP), 1.0));
                    }
                }
            }

            if (DEBUG) {
                System.out.println("Bond formed, cleaved changes 2");
            }

            //Mark Cleaved Bonds in Reactants
            if (bondR != null && (bondR.getProperties().get(BOND_CHANGE_INFORMATION).
                    equals(BOND_CLEAVED)
                    || bondR.getProperties().get(BOND_CHANGE_INFORMATION).
                            equals(PSEUDO_BOND))) {

                if (DEBUG) {
                    System.out.println("Bond formed, cleaved changes 1 - 2 - 1");
                }

                if (!bondR.getAtom(0).getSymbol().equals("PsH")
                        && !bondR.getAtom(1).getSymbol().equals("PsH")) {
                    this.energySum += be.getEnergies(bondR);
                    pEnergy += be.getEnergies(bondR);
                    bondCleavedMap.put(bondR, getMoleculeID(bondR, reaction.getReactants()));
                    bondR.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
                    bondR.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);

                    /*
                     * update reaction center product FP
                     */
                    IAtomContainer moleculeE = getAtomContainer(bondR, reaction.getReactants());

                    if (moleculeE != null && moleculeE.getAtomCount() > 1) {

                        /*
                         * Mark reaction centers
                         */
                        IAtom atomE1 = bondR.getAtom(0);
                        IAtom atomE2 = bondR.getAtom(1);
                        if (!atomE1.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeE, atomE1, REACTANT));
                            setCircularFingerprints(reaction.getID(), moleculeE, atomE1, reactionCenterFormedCleavedFingerprint);
                        }
                        if (!atomE2.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeE, atomE2, REACTANT));
                            setCircularFingerprints(reaction.getID(), moleculeE, atomE2, reactionCenterFormedCleavedFingerprint);
                        }

                        IAtomContainer reactant = getAtomContainer(bondR, reaction.getReactants());
                        IAtomContainer cloneReactant = reactant.getBuilder().newInstance(IAtomContainer.class, reactant);
                        int chippedBondIndex = reactant.indexOf(bondR);
                        totalSmallestFragmentSize += chipTheBondCountSmallestFragmentSize(cloneReactant, chippedBondIndex);
                        formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bondR), 1.0));
                    }
                }
            }
        }

        if (DEBUG) {
            System.out.println("RC Fingerprint");
        }

        /*
         * IMP for RC Fingerprint: compute all the unique reaction centers atoms
         */
        Map<IAtom, IAtom> reactionCenterMap = new HashMap<>();
        bondChangeAnnotator.getReactionCenterSet().stream().filter((atom) -> (!atom.getSymbol().equals("H"))).forEachOrdered((atom) -> {
            reactionCenterMap.put(atom, bondChangeAnnotator.getMappingMap().get(atom));
        });

        if (DEBUG) {
            System.out.println("RC Fingerprint charges like Mg2+ too Mg3+");
        }

        /*
         * Store changes in the charges like Mg2+ too Mg3+
         */
        for (IAtom atom : bondChangeAnnotator.getReactionCenterSet()) {
            if (!atom.getSymbol().equals("H")) {
                IAtomContainer relevantAtomContainer = getRelevantAtomContainer(reaction, atom);

                IAtomContainer relevantAtomContainer1 = getRelevantAtomContainer(reactants, atom);
                IAtomContainer relevantAtomContainer2 = getRelevantAtomContainer(products, atom);
                if (relevantAtomContainer != null && relevantAtomContainer.getAtomCount() == 1) {
                    EnumSubstrateProduct esp = null;

                    if (relevantAtomContainer1 != null) {
                        esp = REACTANT;
                    } else if (relevantAtomContainer2 != null) {
                        esp = PRODUCT;
                    }
                    if (!atom.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, esp));
                        setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterFormedCleavedFingerprint);
                    }
                }
            }

        }

        if (DEBUG) {
            System.out.println("RC Fingerprint ");
        }

        /*
         * Assign Reaction Center Fingerprints
         */
        for (Map.Entry<IAtom, IAtom> mapRC : reactionCenterMap.entrySet()) {

            IAtom sourceAtom = mapRC.getKey();
            IAtom sinkAtom = mapRC.getValue();

            IAtomContainer relevantAtomContainer1 = getRelevantAtomContainer(reaction.getReactants(), sourceAtom);
            IAtomContainer relevantAtomContainer2 = getRelevantAtomContainer(reaction.getProducts(), sinkAtom);

            if (relevantAtomContainer1 != null) {
                for (int i = 0; i < 3; i++) {
                    String circularSMILES = getCircularSMILES(relevantAtomContainer1, sourceAtom, i, true);
                    reactionCenterWFingerprint.add(new Feature(circularSMILES, 1.0));
                }
            }

            if (relevantAtomContainer2 != null) {
                for (int i = 0; i < 3; i++) {
                    String circularSMILES = getCircularSMILES(relevantAtomContainer2, sinkAtom, i, true);
                    reactionCenterWFingerprint.add(new Feature(circularSMILES, 1.0));
                }
            }

            if (relevantAtomContainer1 != null && relevantAtomContainer2 != null) {
                for (int i = 1; i < 4; i++) {
                    String circularSMILESSource = getCircularSMILES(relevantAtomContainer1, sourceAtom, i, true);
                    String circularSMILESSink = getCircularSMILES(relevantAtomContainer2, sinkAtom, i, true);
                    StringBuilder level = new StringBuilder();
                    level.append(circularSMILESSource).append(">>").append(circularSMILESSink);
                    reactionCenterWFingerprint.add(new Feature(level.toString(), 1.0));
                }
                MoleculeMoleculePair molMolPair = getMolMolPair(sourceAtom, sinkAtom, relevantAtomContainer1, relevantAtomContainer2);
                this.reactionMoleculeMoleculePairList.add(molMolPair);
            }
        }

        setEnergyDelta(rEnergy - pEnergy);

        if (DEBUG) {
            System.out.println("Bond Change Calculator END");
        }
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized BEMatrix getEductBEMatrix() {
        return bondChangeAnnotator.getEductBEMatrix();
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized BEMatrix getProductBEMatrix() {
        return bondChangeAnnotator.getProductBEMatrix();
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized RMatrix getRMatrix() {
        return bondChangeAnnotator.getRMatrix();
    }

    /**
     *
     */
    @Override
    public synchronized void printBMatrix() {
        bondChangeAnnotator.printBMatrix();
    }

    /**
     *
     */
    @Override
    public synchronized void printEMatrix() {
        bondChangeAnnotator.printEMatrix();
    }

    /**
     *
     */
    @Override
    public synchronized void printRMatrix() {
        bondChangeAnnotator.printRMatrix();
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeBMatrix(File outputFile) {
        bondChangeAnnotator.writeBMatrix(outputFile);
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeEMatrix(File outputFile) {
        bondChangeAnnotator.writeEMatrix(outputFile);
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeRMatrix(File outputFile) {
        bondChangeAnnotator.writeRMatrix(outputFile);
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized boolean hasRMatrix() {
        return bondChangeAnnotator.hasRMatrix();
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Map<IAtom, IAtom> getMappingMap() {
        return synchronizedMap(bondChangeAnnotator.getMappingMap());
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized List<BondChange> getBondChangeList() {
        return bondChangeAnnotator.getBondChangeList();
    }

    /**
     *
     * @return @throws CDKException
     */
    @Override
    public synchronized IPatternFingerprinter getOrderChangesWFingerprint() throws CDKException {
        return orderChangesWFingerprint;
    }

    @Override
    public synchronized IPatternFingerprinter getStereoChangesWFingerprint() throws CDKException {
        return stereoChangesWFingerprint;
    }

    /**
     *
     * @return @throws CDKException
     */
    @Override
    public synchronized IPatternFingerprinter getFormedCleavedWFingerprint() throws CDKException {
        return formedCleavedWFingerprint;
    }

    /**
     * Return Reaction center Fingerprint
     *
     * @return
     * @throws CDKException
     */
    @Override
    public synchronized IPatternFingerprinter getReactionCenterWFingerprint() throws CDKException {
        return reactionCenterWFingerprint;
    }

    /**
     * Print the atom changes
     *
     * @return
     */
    @Override
    public synchronized String toString() {

        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");

        result.append(NEW_LINE).append(getLicenseHeader());
        result.append(NEW_LINE).append(NEW_LINE).append("//DATA START//");
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Cleaved FingerPrint (Reactant)");
        bondCleavedMap.entrySet().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getID() != null) {
                id1 = bond.getAtom(0).getID();
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getID() != null) {
                id2 = bond.getAtom(1).getID();
                symbol2 = bond.getAtom(1).getSymbol();
            }
            if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Formed FingerPrint (Product)");
        bondFormedMap.entrySet().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getID() != null) {
                id1 = bond.getAtom(0).getID();
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getID() != null) {
                id2 = bond.getAtom(1).getID();
                symbol2 = bond.getAtom(1).getSymbol();
            }
            if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Order Change FingerPrint (Reactant)");
        bondOrderRMap.entrySet().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getID() != null) {
                id1 = bond.getAtom(0).getID();
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getID() != null) {
                id2 = bond.getAtom(1).getID();
                symbol2 = bond.getAtom(1).getSymbol();
            }
            if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Order Change FingerPrint (Product)");
        bondOrderPMap.entrySet().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getID() != null) {
                id1 = bond.getAtom(0).getID();
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getID() != null) {
                id2 = bond.getAtom(1).getID();
                symbol2 = bond.getAtom(1).getSymbol();
            }
            if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Stereo Change FingerPrint (Reactant)");
        AtomStereoRMap.entrySet().forEach((map) -> {
            IAtom atom = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (atom.getID() != null) {
                id1 = atom.getID();
                symbol1 = atom.getSymbol();
            }
            result.append(symbol1).append("(").append(id1).append(")" + "(R/S)" + "\t").append(molID);
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Stereo Change FingerPrint (Product)");
        AtomStereoPMap.entrySet().forEach((map) -> {
            IAtom atom = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (atom.getID() != null) {
                id1 = atom.getID();
                symbol1 = atom.getSymbol();
            }
            result.append(symbol1).append("(").append(id1).append(")" + "(R/S)" + "\t").append(molID);
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("//DATA END//").append(NEW_LINE);
        result.append(NEW_LINE).append(getLicenseFooter());
        result.append(NEW_LINE).append(NEW_LINE);
        return result.toString();
    }

    /**
     * Print the atom changes
     *
     * @param bondChangeInfoFile
     * @throws IOException
     */
    @Override
    public synchronized void writeBondChanges(File bondChangeInfoFile) throws IOException {
        FileWriter bcFW = new FileWriter(bondChangeInfoFile + ".txt");
        try (BufferedWriter bfw = new BufferedWriter(bcFW)) {
            bfw.newLine();

            bfw.write(getLicenseHeader());
            bfw.newLine();
            bfw.write("//DATA START//");
            bfw.newLine();
            bfw.newLine();

            bfw.write("Cleaved FingerPrint (Reactant)");
            for (Map.Entry<IBond, String> map : bondCleavedMap.entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getID() != null) {
                    id1 = bond.getAtom(0).getID();
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getID() != null) {
                    id2 = bond.getAtom(1).getID();
                    symbol2 = bond.getAtom(1).getSymbol();
                }
                if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol2 + "(" + id2 + ")" + "\t" + molID);
                }
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Formed FingerPrint (Product)");
            for (Map.Entry<IBond, String> map : bondFormedMap.entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getID() != null) {
                    id1 = bond.getAtom(0).getID();
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getID() != null) {
                    id2 = bond.getAtom(1).getID();
                    symbol2 = bond.getAtom(1).getSymbol();
                }
                if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol2 + "(" + id2 + ")" + "\t" + molID);
                }
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Order Change FingerPrint (Reactant)");
            for (Map.Entry<IBond, String> map : bondOrderRMap.entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getID() != null) {
                    id1 = bond.getAtom(0).getID();
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getID() != null) {
                    id2 = bond.getAtom(1).getID();
                    symbol2 = bond.getAtom(1).getSymbol();
                }
                if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol2 + "(" + id2 + ")" + "\t" + molID);
                }
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Order Change FingerPrint (Product)");
            for (Map.Entry<IBond, String> map : bondOrderPMap.entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getID() != null) {
                    id1 = bond.getAtom(0).getID();
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getID() != null) {
                    id2 = bond.getAtom(1).getID();
                    symbol2 = bond.getAtom(1).getSymbol();
                }
                if (!symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.isEmpty() && symbol2.isEmpty()) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.isEmpty() && !symbol2.isEmpty()) {
                    bfw.write(symbol2 + "(" + id2 + ")" + "\t" + molID);
                }
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Stereo Change FingerPrint (Reactant)");
            for (Map.Entry<IAtom, String> map : AtomStereoRMap.entrySet()) {
                IAtom atom = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (atom.getID() != null) {
                    id1 = atom.getID();
                    symbol1 = atom.getSymbol();
                }
                if (atom.getID() != null) {
                    id2 = atom.getID();
                    symbol2 = atom.getSymbol();
                }
                bfw.write(symbol1 + "(" + id1 + ")" + "(R/S)"
                        + "\t" + molID);
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Stereo Change FingerPrint (Product)");
            for (Map.Entry<IAtom, String> map : AtomStereoPMap.entrySet()) {
                IAtom atom = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (atom.getID() != null) {
                    id1 = atom.getID();
                    symbol1 = atom.getSymbol();
                }
                bfw.write(symbol1 + "(" + id1 + ")" + "(R/S)"
                        + "\t" + molID);
            }

            bfw.newLine();
            bfw.newLine();
            bfw.write("//DATA END//");
            bfw.flush();

            bfw.newLine();
            bfw.newLine();
            bfw.write(getLicenseFooter());
            bfw.newLine();
            bfw.newLine();
            bfw.flush();
        }
    }

    /**
     * @return atom(s) Formed at reactant side and with reactant ID
     */
    @Override
    public synchronized Map<IBond, String> getBondFormedProduct() {
        return synchronizedMap(bondFormedMap);
    }

    /**
     * @return atom(s) cleaved at reactant side and with reactant ID
     */
    @Override
    public synchronized Map<IBond, String> getBondCleavedReactant() {
        return synchronizedMap(bondCleavedMap);
    }

    /**
     * @return atom(s) order changed at reactant side and with reactant ID
     */
    @Override
    public synchronized Map<IBond, String> getBondOrderReactant() {
        return synchronizedMap(bondOrderRMap);
    }

    /**
     * @return atom(s) order changed at product side and with product ID
     */
    @Override
    public synchronized Map<IBond, String> getBondOrderProduct() {
        return synchronizedMap(bondOrderPMap);
    }

    /**
     * @return the atom effect by Stereo changes at reactant side and with
     * reactant ID
     */
    @Override
    public synchronized Map<IAtom, String> getStereoCenterAtomsReactant() {
        return synchronizedMap(AtomStereoRMap);
    }

    /**
     * @return the atom effect by Stereo changes at product side and with
     * product ID
     */
    @Override
    public synchronized Map<IAtom, String> getStereoCenterAtomsProduct() {
        return synchronizedMap(AtomStereoPMap);
    }

    /**
     *
     * @return (removed the unchanged H atoms)
     */
    @Override
    public synchronized IReaction getReactionWithCompressUnChangedHydrogens() {
        IReaction compressedReaction = null;

        try {
            //Clone reaction with bond change flags
            compressedReaction = deepClone(mappedReaction);
            compressedReaction.setProperties(mappedReaction.getProperties());
            //Add mapping to the clone
            Map<IAtom, IAtom> mappings = new HashMap<>();
            for (IMapping mapping : compressedReaction.mappings()) {
                mappings.put((IAtom) mapping.getChemObject(0), (IAtom) mapping.getChemObject(1));
            }

            for (IAtomContainer mol : compressedReaction.getReactants().atomContainers()) {
                List<IAtom> atoms = getAtoms(mol);
                if (atoms.size() > 1) {
                    atoms.stream().filter((atom) -> (atom.getSymbol().equalsIgnoreCase("H") && mappings.containsKey(atom))).filter((atom) -> (atom.getProperty(BOND_CHANGE_INFORMATION) == null)).forEachOrdered((atom) -> {
                        mol.removeAtom(atom);
                    });
                    /*
                    Do not remove radical changes Hydrogen changes p-sh
                     */
                    CDKHydrogenAdder hAdder = getInstance(mol.getBuilder());
                    try {
                        hAdder.addImplicitHydrogens(mol);
                        kekulize(mol);
                    } catch (CDKException ex) {
                    }
                }
            }

            for (IAtomContainer mol : compressedReaction.getProducts().atomContainers()) {
                List<IAtom> atoms = getAtoms(mol);
                if (atoms.size() > 1) {
                    atoms.stream().filter((atom) -> (atom.getSymbol().equalsIgnoreCase("H") && mappings.containsValue(atom))).filter((atom) -> (atom.getProperty(BOND_CHANGE_INFORMATION) == null)).forEachOrdered((atom) -> {
                        mol.removeAtom(atom);
                    });
                    /*
                    Do not remove radical changes Hydrogen changes p-sh
                     */
                    CDKHydrogenAdder cdkHAdder = getInstance(mol.getBuilder());
                    try {
                        cdkHAdder.addImplicitHydrogens(mol);
                        kekulize(mol);
                    } catch (CDKException ex) {
                    }
                }
            }

            /*
             * This is a VERY important step else mapping will point to NULL objects
             */
            cleanMapping(compressedReaction);

            for (Map.Entry<IAtom, IAtom> map : mappings.entrySet()) {
                if (map.getKey() != null && map.getValue() != null) {
                    map.getKey().setFlag(MAPPED, true);
                    map.getValue().setFlag(MAPPED, true);
                    compressedReaction.addMapping(new Mapping(map.getKey(), map.getValue()));
                }
            }
        } catch (CloneNotSupportedException ex) {
            LOGGER.error(SEVERE, null, ex);
        }

        try {
            KekulizeReaction(compressedReaction);
        } catch (CDKException ex) {
            LOGGER.error(Level.SEVERE, null, ex);
        }
        return compressedReaction;
    }

    private static void KekulizeReaction(IReaction r) throws CDKException {
        ElectronDonation model = ElectronDonation.daylight();
//        CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.all(6));
//        Aromaticity aromaticity = new Aromaticity(model, cycles);

        Aromaticity aromaticity = new Aromaticity(model,
                Cycles.or(Cycles.all(),
                        Cycles.or(Cycles.relevant(),
                                Cycles.essential())));
        // apply our configured model to each molecule
        for (IAtomContainer molecule : r.getReactants().atomContainers()) {
            aromaticity.apply(molecule);
        }
        for (IAtomContainer molecule : r.getProducts().atomContainers()) {
            aromaticity.apply(molecule);
        }
    }

    /**
     *
     * @param MappedReaction
     */
    private synchronized void cleanMapping(IReaction MappedReaction) {

        int count = MappedReaction.getMappingCount();
        for (int i = count; i > 0; i--) {
            MappedReaction.removeMapping(i);
        }

        for (int eMol = 0; eMol < MappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = MappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atomEMap = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                atomEMap.setFlag(MAPPED, false);
            }
        }

        for (int pMol = 0; pMol < MappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = MappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atomPMap = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                atomPMap.setFlag(MAPPED, false);
            }
        }
    }

    /**
     *
     * @return (removed the unchanged H atoms)
     * @throws Exception
     */
    @Override
    public synchronized IReaction getReaction() throws Exception {
        IReaction mappedReactionWithBondChanges = mappedReaction;
        return mappedReactionWithBondChanges;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Map<IAtom, IAtom> getAtomAtomMappings() {
        return bondChangeAnnotator.getMappingMap();
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMappingContainer getMappingContainer() {
        return bondChangeAnnotator.getMappingContainer();
    }

    /**
     *
     * @return @throws CDKException
     */
    public synchronized int getTotalBondBreakingEnergy() throws CDKException {
        return this.energySum;
    }

    /**
     * @return the energyDelta
     */
    public synchronized int getEnergyDelta() {
        return abs(energyDelta);
    }

    /**
     * @param energyDelta the energyDelta to set
     */
    private void setEnergyDelta(int energyDelta) {
        this.energyDelta = energyDelta;
    }

    /**
     *
     * @return
     */
    @Override
    public List<AtomStereoChangeInformation> getStereoChangeList() {
        return bondChangeAnnotator.getStereoChangeList();
    }

    /**
     *
     * @return
     */
    @Override
    public Iterable<AtomStereoChangeInformation> getConformationChangeList() {
        return bondChangeAnnotator.getConformationChangeList();
    }

    private String getLicenseHeader() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");

        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("ecBLAST (Enzymatic Reaction BLAST)");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("Contact: Syed Asad Rahman,");
        result.append(NEW_LINE).append("\t EMBL-EBI, Hinxton ");
        result.append(NEW_LINE).append("\t Cambridge CB10 1SD");
        result.append(NEW_LINE).append("\t United Kingdom ");
        result.append(NEW_LINE).append("e-mail: asad@ebi.ac.uk|s9asad@gmail.com, thornton@ebi.ac.uk");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("ecBLAST software can perform atom-atom mapping,");
        result.append(NEW_LINE).append("marks bond changes between reactions, calculate");
        result.append(NEW_LINE).append("similarity between small molecules, reactions");
        result.append(NEW_LINE).append("using in-house algorithm developed at EMBL-EBI.");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("Acknowledgment: Many thanks to Franz Fenninger,");
        result.append(NEW_LINE).append("Gilleain Torrance, Lorenzo Baldacci and Gemma L.");
        result.append(NEW_LINE).append("Holliday for their contributions/input.").append(NEW_LINE);
        result.append(NEW_LINE).append("The open source tools i.e. the SMSD, CDK and InChI");
        result.append(NEW_LINE).append("were used in this project. The licensed version of");
        result.append(NEW_LINE).append("the Chemaxon was used for 3D model building,");
        result.append(NEW_LINE).append("cleaning and hydrogen handling.");
        result.append(NEW_LINE).append("");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        return result.toString();
    }

    private String getLicenseFooter() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("NOTE: You can't distribute this tool or it's");
        result.append(NEW_LINE).append("components, including the output without prior");
        result.append(NEW_LINE).append("written approval from the authors. *Copyright");
        result.append(NEW_LINE).append("and limitations* of use except where indicated,");
        result.append(NEW_LINE).append("all content on this tool, including programs,");
        result.append(NEW_LINE).append("code, concept, text and collectively the \"content\"");
        result.append(NEW_LINE).append("is the property of authors and is protected by");
        result.append(NEW_LINE).append("copyright and other intellectual property laws.").append(NEW_LINE);
        result.append(NEW_LINE).append("Certain content might be licensed on a non-exclusive");
        result.append(NEW_LINE).append("basis from 3rd Parties. This is a beta version;");
        result.append(NEW_LINE).append("suggestions to improve this tool are welcome.");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        return result.toString();
    }

    /**
     * @return the reactionCenterFormedCleavedFingerprint
     */
    public Map<Integer, IPatternFingerprinter> getReactionCenterFormedCleavedFingerprint() {
        return synchronizedMap(reactionCenterFormedCleavedFingerprint);
    }

    /**
     * @return the reactionCenterOrderChangeFingerprint
     */
    public Map<Integer, IPatternFingerprinter> getReactionCenterOrderChangeFingerprint() {
        return synchronizedMap(reactionCenterOrderChangeFingerprint);
    }

    /**
     * @return the reactionCenterStereoChangeFingerprint
     */
    public Map<Integer, IPatternFingerprinter> getReactionCenterStereoChangeFingerprint() {
        return synchronizedMap(reactionCenterStereoChangeFingerprint);
    }

    /**
     *
     * @return
     */
    @Override
    public Collection<IAtom> getReactionCenterSet() {
        return bondChangeAnnotator.getReactionCenterSet();
    }

    /**
     * @return the Reaction Center Fragment List
     */
    @Override
    public Collection<ReactionCenterFragment> getReactionCenterFragmentList() {
        return unmodifiableCollection(reactionCenterFragmentList);
    }

    @Override
    public Collection<MoleculeMoleculePair> getReactionCentreTransformationPairs() {
        return unmodifiableCollection(reactionMoleculeMoleculePairList);
    }

    /**
     *
     * @return
     */
    @Override
    public Map<String, Collection<String>> getMoleculeMoleculeTransformationPairs() {
        Map<String, Collection<String>> uniqueRPAIRS = new TreeMap<>();
        this.getReactionCentreTransformationPairs().stream().map((m) -> {
            if (!uniqueRPAIRS.containsKey(m.getName().toString())) {
                LinkedList<String> l = new LinkedList<>();
                uniqueRPAIRS.put(m.getName().toString(), l);
            }
            return m;
        }).forEachOrdered((m) -> {
            uniqueRPAIRS.get(m.getName().toString()).add(m.getSmirks());
        });
        return uniqueRPAIRS;
    }

    private int chipTheBondCountSmallestFragmentSize(IAtomContainer cloneContainer, int chippedBondIndex) {
        int size = cloneContainer.getAtomCount();
        cloneContainer.removeBond(chippedBondIndex);
        boolean fragmentFlag = isConnected(cloneContainer);
        if (!fragmentFlag) {
            IAtomContainerSet partitionIntoMolecules = partitionIntoMolecules(cloneContainer);
            for (IAtomContainer ac : partitionIntoMolecules.atomContainers()) {
                if (size > ac.getAtomCount()) {
                    size = ac.getAtomCount();
                }
            }
        }
        return size;
    }

    /**
     * @return the totalSmallestFragmentSize
     */
    public int getTotalSmallestFragmentSize() {
        return totalSmallestFragmentSize;
    }
}
