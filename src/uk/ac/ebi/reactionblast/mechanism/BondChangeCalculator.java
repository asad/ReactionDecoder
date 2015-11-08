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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
import org.openscience.smsd.tools.BondEnergies;
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
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct;
import uk.ac.ebi.reactionblast.mechanism.interfaces.IChangeCalculator;
import uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool;

/**
 * This class marks the bond changes
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondChangeCalculator extends AbstractChangeCalculator implements IChangeCalculator {

    private static final long serialVersionUID = 98698690880809981L;
    private static final Logger LOG = Logger.getLogger(BondChangeCalculator.class.getName());
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
        int rEnergy = 0;
        int pEnergy = 0;

        this.energySum = 0;
        this.energyDelta = 0;
        this.totalSmallestFragmentSize = 0;
        this.mappedReaction = reaction;
        this.bondChangeAnnotator = new BondChangeAnnotator(this.mappedReaction, true, generate2D, generate3D);
        BondEnergies be = BondEnergies.getInstance();

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

        this.bondFormedMap = Collections.synchronizedMap(new HashMap<IBond, String>());
        this.bondCleavedMap = Collections.synchronizedMap(new HashMap<IBond, String>());
        this.bondOrderRMap = Collections.synchronizedMap(new HashMap<IBond, String>());
        this.bondOrderPMap = Collections.synchronizedMap(new HashMap<IBond, String>());
        this.AtomStereoRMap = Collections.synchronizedMap(new HashMap<IAtom, String>());
        this.AtomStereoPMap = Collections.synchronizedMap(new HashMap<IAtom, String>());

        this.reactionCenterFragmentList = Collections.synchronizedList(new ArrayList<ReactionCenterFragment>());

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
                atomConformation.getReactantAtom().setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO);
                AtomStereoRMap.put(atomConformation.getReactantAtom(), getMoleculeID(atomConformation.getReactantAtom(), reaction.getReactants()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomR1 = atomConformation.getReactantAtom();
                IAtomContainer moleculeR = getAtomContainer(atomConformation.getReactantAtom(), reaction.getReactants());

                if (moleculeR.getAtomCount() > 1) {
                    if (!atomR1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeR, atomR1, EnumSubstrateProduct.REACTANT));
                        setCircularFingerprints(reaction.getID(), moleculeR, atomR1, reactionCenterStereoChangeFingerprint);
                    }
                }
            }
            if (atomConformation.getProductAtom() != null) {
                atomConformation.getProductAtom().setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO);
                AtomStereoPMap.put(atomConformation.getProductAtom(), getMoleculeID(atomConformation.getProductAtom(), reaction.getProducts()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomP1 = atomConformation.getProductAtom();
                IAtomContainer moleculeP = getAtomContainer(atomConformation.getProductAtom(), reaction.getProducts());

                if (moleculeP.getAtomCount() > 1) {
                    if (!atomP1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, EnumSubstrateProduct.PRODUCT));
                        setCircularFingerprints(reaction.getID(), moleculeP, atomP1, reactionCenterStereoChangeFingerprint);
                    }
                }
            }
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
                atomStereo.getReactantAtom().setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO);
                AtomStereoRMap.put(atomStereo.getReactantAtom(), getMoleculeID(atomStereo.getReactantAtom(), reaction.getReactants()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomR1 = atomStereo.getReactantAtom();
                IAtomContainer moleculeR = getAtomContainer(atomStereo.getReactantAtom(), reaction.getReactants());

                if (moleculeR.getAtomCount() > 1) {

                    if (!atomR1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeR, atomR1, EnumSubstrateProduct.REACTANT));
                        setCircularFingerprints(reaction.getID(), moleculeR, atomR1, reactionCenterStereoChangeFingerprint);
                    }
                }
            }
            if (atomStereo.getProductAtom() != null) {
                atomStereo.getProductAtom().setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO);
                AtomStereoPMap.put(atomStereo.getProductAtom(), getMoleculeID(atomStereo.getProductAtom(), reaction.getProducts()));

                /*
                 * Update Reaction center FP
                 */
                IAtom atomP1 = atomStereo.getProductAtom();
                IAtomContainer moleculeP = getAtomContainer(atomStereo.getProductAtom(), reaction.getProducts());

                if (moleculeP.getAtomCount() > 1) {

                    if (!atomP1.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, EnumSubstrateProduct.PRODUCT));
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
        for (BondChange bcinfo : bondChangeAnnotator.getBondChangeList()) {
            IBond bondR = bcinfo.getReactantBond();
            IBond bondP = bcinfo.getProductBond();

            // Mark Bond Order Changes
            if (bondR != null && bondP != null
                    && bondP.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER)
                    && bondR.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER)) {

                bondOrderRMap.put(bondR, getMoleculeID(bondR, reaction.getReactants()));
                bondR.getAtom(0).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                bondR.getAtom(1).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);

                bondOrderPMap.put(bondP, getMoleculeID(bondP, reaction.getProducts()));
                bondP.getAtom(0).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                bondP.getAtom(1).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);

                reactantAtoms.add(bondR.getAtom(0));
                reactantAtoms.add(bondR.getAtom(1));

                productAtoms.add(bondP.getAtom(0));
                productAtoms.add(bondP.getAtom(1));
                orderChangesWFingerprint.add(new Feature(getCanonisedBondChangePattern(bondR, bondP), 1.0));
            }
        }


        /*
         * Store changes in the bond order
         */
        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet products = reaction.getProducts();

        for (IAtom atom : reactantAtoms) {
            IAtomContainer relevantAtomContainer = AtomContainerSetManipulator.getRelevantAtomContainer(reactants, atom);
            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, EnumSubstrateProduct.REACTANT));
            setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterOrderChangeFingerprint);
        }

        for (IAtom atom : productAtoms) {
            IAtomContainer relevantAtomContainer = AtomContainerSetManipulator.getRelevantAtomContainer(products, atom);
            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, EnumSubstrateProduct.PRODUCT));
            setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterOrderChangeFingerprint);
        }


        /*
         * Loop for formed, cleaved changes
         */
        for (BondChange bcinfo : bondChangeAnnotator.getBondChangeList()) {
            IBond bondR = bcinfo.getReactantBond();
            IBond bondP = bcinfo.getProductBond();

            //Mark Formed Bonds in the Product
            if (bondP != null && (bondP.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED)
                    || bondP.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND))) {

                if (!bondP.getAtom(0).getSymbol().equals("PsH")
                        && !bondP.getAtom(1).getSymbol().equals("PsH")) {
                    this.energySum += be.getEnergies(bondP);
                    pEnergy += be.getEnergies(bondP);
                    bondFormedMap.put(bondP, getMoleculeID(bondP, reaction.getProducts()));
                    bondP.getAtom(0).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    bondP.getAtom(1).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);

                    /*
                     * Update Reaction center FP
                     */
                    IAtomContainer moleculeP = getAtomContainer(bondP, reaction.getProducts());

                    if (moleculeP.getAtomCount() > 1) {

                        /*
                         * Mark reaction centers
                         */
                        IAtom atomP1 = bondP.getAtom(0);
                        IAtom atomP2 = bondP.getAtom(1);
                        if (!atomP1.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP1, EnumSubstrateProduct.PRODUCT));
                            setCircularFingerprints(reaction.getID(), moleculeP, atomP1, reactionCenterFormedCleavedFingerprint);
                        }
                        if (!atomP2.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeP, atomP2, EnumSubstrateProduct.PRODUCT));
                            setCircularFingerprints(reaction.getID(), moleculeP, atomP2, reactionCenterFormedCleavedFingerprint);
                        }
                    }
                    IAtomContainer product = getAtomContainer(bondP, reaction.getProducts());
                    IAtomContainer cloneProduct = product.getBuilder().newInstance(IAtomContainer.class, product);
                    int chippedBondIndex = product.getBondNumber(bondP);
                    totalSmallestFragmentSize += chipTheBondCountSmallestFragmentSize(cloneProduct, chippedBondIndex);
                    formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bondP), 1.0));
                }
            }

            //Mark Cleaved Bonds in Reactants
            if (bondR != null && (bondR.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED)
                    || bondR.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION).
                    equals(ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND))) {

                if (!bondR.getAtom(0).getSymbol().equals("PsH")
                        && !bondR.getAtom(1).getSymbol().equals("PsH")) {
                    this.energySum += be.getEnergies(bondR);
                    pEnergy += be.getEnergies(bondR);
                    bondCleavedMap.put(bondR, getMoleculeID(bondR, reaction.getReactants()));
                    bondR.getAtom(0).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                    bondR.getAtom(1).setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);

                    /*
                     * update reaction center product FP
                     */
                    IAtomContainer moleculeE = getAtomContainer(bondR, reaction.getReactants());

                    if (moleculeE.getAtomCount() > 1) {

                        /*
                         * Mark reaction centers
                         */
                        IAtom atomE1 = bondR.getAtom(0);
                        IAtom atomE2 = bondR.getAtom(1);
                        if (!atomE1.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeE, atomE1, EnumSubstrateProduct.REACTANT));
                            setCircularFingerprints(reaction.getID(), moleculeE, atomE1, reactionCenterFormedCleavedFingerprint);
                        }
                        if (!atomE2.getSymbol().equals("H")) {
                            reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(moleculeE, atomE2, EnumSubstrateProduct.REACTANT));
                            setCircularFingerprints(reaction.getID(), moleculeE, atomE2, reactionCenterFormedCleavedFingerprint);
                        }
                    }
                    IAtomContainer reactant = getAtomContainer(bondR, reaction.getReactants());
                    IAtomContainer cloneReactant = reactant.getBuilder().newInstance(IAtomContainer.class, reactant);
                    int chippedBondIndex = reactant.getBondNumber(bondR);
                    totalSmallestFragmentSize += chipTheBondCountSmallestFragmentSize(cloneReactant, chippedBondIndex);
                    formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bondR), 1.0));
                }
            }
        }


        /*
         * IMP for RC Fingerprint: Mine all the unique reaction centers atoms
         */
        Map<IAtom, IAtom> reactionCenterMap = new HashMap<>();
        for (IAtom atom : bondChangeAnnotator.getReactionCenterSet()) {
            if (!atom.getSymbol().equals("H")) {
                reactionCenterMap.put(atom, bondChangeAnnotator.getMappingMap().get(atom));
            }
        }

        /*
         * Store changes in the charges like Mg2+ too Mg3+
         */
        for (IAtom atom : bondChangeAnnotator.getReactionCenterSet()) {
            if (!atom.getSymbol().equals("H")) {
                IAtomContainer relevantAtomContainer = ReactionManipulator.getRelevantAtomContainer(reaction, atom);

                IAtomContainer relevantAtomContainer1 = AtomContainerSetManipulator.getRelevantAtomContainer(reactants, atom);
                IAtomContainer relevantAtomContainer2 = AtomContainerSetManipulator.getRelevantAtomContainer(products, atom);
                if (relevantAtomContainer != null && relevantAtomContainer.getAtomCount() == 1) {
                    EnumSubstrateProduct esp = null;

                    if (relevantAtomContainer1 != null) {
                        esp = EnumSubstrateProduct.REACTANT;
                    } else if (relevantAtomContainer2 != null) {
                        esp = EnumSubstrateProduct.PRODUCT;
                    }
                    if (!atom.getSymbol().equals("H")) {
                        reactionCenterFragmentList.addAll(getCircularReactionPatternFingerprints(relevantAtomContainer, atom, esp));
                        setCircularFingerprints(reaction.getID(), relevantAtomContainer, atom, reactionCenterFormedCleavedFingerprint);
                    }
                }
            }

        }

        /*
         * Assign Reaction Center Fingerprints
         */
        for (Map.Entry<IAtom, IAtom> mapRC : reactionCenterMap.entrySet()) {

            IAtom sourceAtom = mapRC.getKey();
            IAtom sinkAtom = mapRC.getValue();

            IAtomContainer relevantAtomContainer1 = AtomContainerSetManipulator.getRelevantAtomContainer(reaction.getReactants(), sourceAtom);
            IAtomContainer relevantAtomContainer2 = AtomContainerSetManipulator.getRelevantAtomContainer(reaction.getProducts(), sinkAtom);

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
    }

    @Override
    public synchronized BEMatrix getEductBEMatrix() {
        return bondChangeAnnotator.getEductBEMatrix();
    }

    @Override
    public synchronized BEMatrix getProductBEMatrix() {
        return bondChangeAnnotator.getProductBEMatrix();
    }

    @Override
    public synchronized RMatrix getRMatrix() {
        return bondChangeAnnotator.getRMatrix();
    }

    @Override
    public synchronized void printBMatrix() {
        bondChangeAnnotator.printBMatrix();
    }

    @Override
    public synchronized void printEMatrix() {
        bondChangeAnnotator.printEMatrix();
    }

    @Override
    public synchronized void printRMatrix() {
        bondChangeAnnotator.printRMatrix();
    }

    @Override
    public synchronized void writeBMatrix(File outputFile) {
        bondChangeAnnotator.writeBMatrix(outputFile);
    }

    @Override
    public synchronized void writeEMatrix(File outputFile) {
        bondChangeAnnotator.writeEMatrix(outputFile);
    }

    @Override
    public synchronized void writeRMatrix(File outputFile) {
        bondChangeAnnotator.writeRMatrix(outputFile);
    }

    @Override
    public synchronized boolean hasRMatrix() {
        return bondChangeAnnotator.hasRMatrix();
    }

    @Override
    public synchronized Map<IAtom, IAtom> getMappingMap() {
        return Collections.synchronizedMap(bondChangeAnnotator.getMappingMap());
    }

    @Override
    public synchronized List<BondChange> getBondChangeList() {
        return bondChangeAnnotator.getBondChangeList();
    }

    @Override
    public synchronized IPatternFingerprinter getOrderChangesWFingerprint() throws CDKException {
        return orderChangesWFingerprint;
    }

    @Override
    public synchronized IPatternFingerprinter getStereoChangesWFingerprint() throws CDKException {
        return stereoChangesWFingerprint;
    }

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
        String NEW_LINE = System.getProperty("line.separator");

        result.append(NEW_LINE).append(getLicenseHeader());
        result.append(NEW_LINE).append(NEW_LINE).append("//DATA START//");
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Cleaved FingerPrint (Reactant)");
        for (Map.Entry<IBond, String> map : bondCleavedMap.entrySet()) {
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
            if (!symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.equals("") && symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        }
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Formed FingerPrint (Product)");
        for (Map.Entry<IBond, String> map : bondFormedMap.entrySet()) {
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
            if (!symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.equals("") && symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        }
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Order Change FingerPrint (Reactant)");
        for (Map.Entry<IBond, String> map : bondOrderRMap.entrySet()) {
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
            if (!symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.equals("") && symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        }
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Order Change FingerPrint (Product)");
        for (Map.Entry<IBond, String> map : bondOrderPMap.entrySet()) {
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
            if (!symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")");
                result.append(getBondOrderSign(bond)).append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            } else if (!symbol1.equals("") && symbol2.equals("")) {
                result.append(symbol1);
                result.append("(").append(id1).append(")" + "\t").append(molID);
            } else if (symbol1.equals("") && !symbol2.equals("")) {
                result.append(symbol2);
                result.append("(").append(id2).append(")" + "\t").append(molID);
            }
        }
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Stereo Change FingerPrint (Reactant)");
        for (Map.Entry<IAtom, String> map : AtomStereoRMap.entrySet()) {
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
        }
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Stereo Change FingerPrint (Product)");
        for (Map.Entry<IAtom, String> map : AtomStereoPMap.entrySet()) {
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
        }
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
                if (!symbol1.equals("") && !symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.equals("") && symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.equals("") && !symbol2.equals("")) {
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
                if (!symbol1.equals("") && !symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.equals("") && symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.equals("") && !symbol2.equals("")) {
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
                if (!symbol1.equals("") && !symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.equals("") && symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.equals("") && !symbol2.equals("")) {
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
                if (!symbol1.equals("") && !symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + getBondOrderSign(bond)
                            + symbol2 + "(" + id2 + ")" + "\t" + molID);
                } else if (!symbol1.equals("") && symbol2.equals("")) {
                    bfw.write(symbol1 + "(" + id1 + ")" + "\t" + molID);
                } else if (symbol1.equals("") && !symbol2.equals("")) {
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
        return Collections.synchronizedMap(bondFormedMap);
    }

    /**
     * @return atom(s) cleaved at reactant side and with reactant ID
     */
    @Override
    public synchronized Map<IBond, String> getBondCleavedReactant() {
        return Collections.synchronizedMap(bondCleavedMap);
    }

    /**
     * @return atom(s) order changed at reactant side and with reactant ID
     */
    @Override
    public synchronized Map<IBond, String> getBondOrderReactant() {
        return Collections.synchronizedMap(bondOrderRMap);
    }

    /**
     * @return atom(s) order changed at product side and with product ID
     */
    @Override
    public synchronized Map<IBond, String> getBondOrderProduct() {
        return Collections.synchronizedMap(bondOrderPMap);
    }

    /**
     * @return the atom effect by Stereo changes at reactant side and with
     * reactant ID
     */
    @Override
    public synchronized Map<IAtom, String> getStereoCenterAtomsReactant() {
        return Collections.synchronizedMap(AtomStereoRMap);
    }

    /**
     * @return the atom effect by Stereo changes at product side and with
     * product ID
     */
    @Override
    public synchronized Map<IAtom, String> getStereoCenterAtomsProduct() {
        return Collections.synchronizedMap(AtomStereoPMap);
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
            compressedReaction = ExtReactionManipulatorTool.deepClone(mappedReaction);
            compressedReaction.setProperties(mappedReaction.getProperties());
            //Add mapping to the clone
            Map<IAtom, IAtom> mappings = new HashMap<>();
            for (IMapping mapping : compressedReaction.mappings()) {
                mappings.put((IAtom) mapping.getChemObject(0), (IAtom) mapping.getChemObject(1));
            }

            for (IAtomContainer mol : compressedReaction.getReactants().atomContainers()) {
                List<IAtom> atoms = getAtoms(mol);
                if (atoms.size() > 1) {
                    for (IAtom atom : atoms) {
                        /*
                         Do not remove radical changes Hydrogen changes p-sh
                         */
                        if (atom.getSymbol().equalsIgnoreCase("H") && mappings.containsKey(atom)) {
                            if (atom.getProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == null) {
                                mol.removeAtomAndConnectedElectronContainers(atom);
                            }
                        }
                    }
                    CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
                    try {
                        hAdder.addImplicitHydrogens(mol);
                        Kekulization.kekulize(mol);
                    } catch (CDKException ex) {
                    }
                }
            }

            for (IAtomContainer mol : compressedReaction.getProducts().atomContainers()) {
                List<IAtom> atoms = getAtoms(mol);
                if (atoms.size() > 1) {
                    for (IAtom atom : atoms) {
                        /*
                         Do not remove radical changes Hydrogen changes p-sh
                         */
                        if (atom.getSymbol().equalsIgnoreCase("H") && mappings.containsValue(atom)) {
                            if (atom.getProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == null) {
                                mol.removeAtomAndConnectedElectronContainers(atom);
                            }
                        }
                    }
                    CDKHydrogenAdder cdkHAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
                    try {
                        cdkHAdder.addImplicitHydrogens(mol);
                        Kekulization.kekulize(mol);
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
                    map.getKey().setFlag(CDKConstants.MAPPED, true);
                    map.getValue().setFlag(CDKConstants.MAPPED, true);
                    compressedReaction.addMapping(new Mapping(map.getKey(), map.getValue()));
                }
            }
        } catch (CloneNotSupportedException ex) {
            Logger.getLogger(BondChangeCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
        return compressedReaction;
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
                atomEMap.setFlag(CDKConstants.MAPPED, false);
            }
        }

        for (int pMol = 0; pMol < MappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = MappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atomPMap = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                atomPMap.setFlag(CDKConstants.MAPPED, false);
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

    @Override
    public synchronized Map<IAtom, IAtom> getAtomAtomMappings() {
        return bondChangeAnnotator.getMappingMap();
    }

    @Override
    public synchronized AtomAtomMappingContainer getMappingContainer() {
        return bondChangeAnnotator.getMappingContainer();
    }

    public synchronized int getTotalBondBreakingEnergy() throws CDKException {
        return this.energySum;
    }

    /**
     * @return the energyDelta
     */
    public synchronized int getEnergyDelta() {
        return Math.abs(energyDelta);
    }

    /**
     * @param energyDelta the energyDelta to set
     */
    private void setEnergyDelta(int energyDelta) {
        this.energyDelta = energyDelta;
    }

    @Override
    public List<AtomStereoChangeInformation> getStereoChangeList() {
        return bondChangeAnnotator.getStereoChangeList();
    }

    @Override
    public Iterable<AtomStereoChangeInformation> getConformationChangeList() {
        return bondChangeAnnotator.getConformationChangeList();
    }

    private String getLicenseHeader() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = System.getProperty("line.separator");

        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("ecBLAST (Enzymatic Reaction BLAST)");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("Contact: Syed Asad Rahman,");
        result.append(NEW_LINE).append("\t EMBL-EBI, Hinxton ");
        result.append(NEW_LINE).append("\t Cambridge CB10 1SD");
        result.append(NEW_LINE).append("\t United Kingdom ");
        result.append(NEW_LINE).append("e-mail: asad@ebi.ac.uk, thornton@ebi.ac.uk");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("ecBLAST software can perform atom-atom mapping,");
        result.append(NEW_LINE).append("marks bond changes between reactions, calculate");
        result.append(NEW_LINE).append("similarity between small molecules, reactions");
        result.append(NEW_LINE).append("using in-house algorithm developed at EMBL-EBI.");
        result.append(NEW_LINE).append("++++++++++++++++++++++++++++++++++++++++++++++").append(NEW_LINE);
        result.append(NEW_LINE).append("Acknowledgment: Many thanks to Franz Fenninger,");
        result.append(NEW_LINE).append("Gilleain Torrance, Lorenzo Baldacci and Gemma L.");
        result.append(NEW_LINE).append("Holliday for their contributions/inputs.").append(NEW_LINE);
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
        String NEW_LINE = System.getProperty("line.separator");
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
        return Collections.synchronizedMap(reactionCenterFormedCleavedFingerprint);
    }

    /**
     * @return the reactionCenterOrderChangeFingerprint
     */
    public Map<Integer, IPatternFingerprinter> getReactionCenterOrderChangeFingerprint() {
        return Collections.synchronizedMap(reactionCenterOrderChangeFingerprint);
    }

    /**
     * @return the reactionCenterStereoChangeFingerprint
     */
    public Map<Integer, IPatternFingerprinter> getReactionCenterStereoChangeFingerprint() {
        return Collections.synchronizedMap(reactionCenterStereoChangeFingerprint);
    }

    @Override
    public Collection<IAtom> getReactionCenterSet() {
        return bondChangeAnnotator.getReactionCenterSet();
    }

    /**
     * @return the Reaction Center Fragment List
     */
    @Override
    public Collection<ReactionCenterFragment> getReactionCenterFragmentList() {
        return Collections.unmodifiableCollection(reactionCenterFragmentList);
    }

    @Override
    public Collection<MoleculeMoleculePair> getReactionCentreTransformationPairs() {
        return Collections.unmodifiableCollection(reactionMoleculeMoleculePairList);
    }

    @Override
    public Map<String, Collection<String>> getMoleculeMoleculeTransformationPairs() {
        Map<String, Collection<String>> uniqueRPAIRS = new TreeMap<>();
        for (MoleculeMoleculePair m : this.getReactionCentreTransformationPairs()) {
            if (!uniqueRPAIRS.containsKey(m.getName().toString())) {
                LinkedList<String> l = new LinkedList<>();
                uniqueRPAIRS.put(m.getName().toString(), l);
            }
            uniqueRPAIRS.get(m.getName().toString()).add(m.getSmirks());
        }
        return uniqueRPAIRS;
    }

    private int chipTheBondCountSmallestFragmentSize(IAtomContainer cloneContainer, int chippedBondIndex) {
        int size = cloneContainer.getAtomCount();
        cloneContainer.removeBond(chippedBondIndex);
        boolean fragmentFlag = ConnectivityChecker.isConnected(cloneContainer);
        if (!fragmentFlag) {
            IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(cloneContainer);
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
