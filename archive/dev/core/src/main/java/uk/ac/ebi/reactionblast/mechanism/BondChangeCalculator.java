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

import uk.ac.ebi.reactionblast.mechanism.interfaces.IBondChangeCalculator;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.unmodifiableCollection;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import static org.openscience.cdk.graph.ConnectivityChecker.partitionIntoMolecules;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.smsd.tools.BondEnergies;
import static org.openscience.smsd.tools.BondEnergies.getInstance;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.helper.MoleculeMoleculePair;
import uk.ac.ebi.reactionblast.mechanism.interfaces.AbstractChangeCalculator;
import static java.util.Collections.synchronizedMap;
import java.util.HashMap;
import static org.openscience.cdk.CDKConstants.REACTIVE_CENTER;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;
import uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.MAPPED;
import static org.openscience.cdk.aromaticity.Kekulization.kekulize;
import org.openscience.cdk.interfaces.IMapping;
import static org.openscience.cdk.tools.CDKHydrogenAdder.getInstance;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.BOND_CHANGE_INFORMATION;
import static uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool.deepClone;
import static java.lang.Math.abs;
import static java.lang.System.getProperty;
import static java.util.logging.Logger.getLogger;

/**
 * This class marks the bond changes
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondChangeCalculator extends AbstractChangeCalculator implements IBondChangeCalculator {

    private final boolean DEBUG = false;
    private static final long serialVersionUID = 98698690880809981L;
    private static final Logger LOG = getLogger(BondChangeCalculator.class.getName());

    private final IPatternFingerprinter formedCleavedWFingerprint;
    private final IPatternFingerprinter orderChangesWFingerprint;
    private final IPatternFingerprinter stereoChangesWFingerprint;
    private final IPatternFingerprinter reactionCenterWFingerprint;
    private final Set<MoleculeMoleculePair> reactionMoleculeMoleculePairList;

    private final Map<IAtom, IAtom> mappings;
    private final IReaction mappedReaction;
    private int energySum;
    private int energyDelta;
    private int totalSmallestFragmentSize;

    private final Set<IBond> bondCleavedFormedChanges;
    private final Map<IBond, IBond> bondOrderChanges;
    private final Set<IAtom> atomStereoChanges;

    /**
     *
     * @param reaction
     * @throws Exception
     */
    public BondChangeCalculator(IReaction reaction) throws Exception {
        if (DEBUG) {
            System.out.println("Bond Change Calculator START");
        }

        BondEnergies be = getInstance();
        int rEnergy = 0;
        int pEnergy = 0;

        this.energySum = 0;
        this.energyDelta = 0;
        this.totalSmallestFragmentSize = 0;
        this.mappedReaction = reaction;
        if (DEBUG) {
            System.out.println("Bond Change Annotator START");
        }
        if (DEBUG) {
            System.out.println("stage 1");
        }
        /*
         * Mappings withou Hydrogens
         */
        this.mappings = ReactionMappingUtility.getMappings(mappedReaction);
        if (DEBUG) {
            System.out.println("stage 2");
        }
        this.bondCleavedFormedChanges = ReactionMappingUtility.getBondCleavedFormedChanges(mappedReaction, mappings);
        if (DEBUG) {
            System.out.println("stage 3");
        }
        this.bondOrderChanges = ReactionMappingUtility.getBondOrderChanges(mappedReaction, mappings);
        if (DEBUG) {
            System.out.println("stage 4");
        }
        this.atomStereoChanges = ReactionMappingUtility.getAtomStereoChanges(mappedReaction, mappings);
        if (DEBUG) {
            System.out.println("stage 5");
        }

        if (DEBUG) {
            System.out.println("Bond Change Annotator END");
        }

        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);

        if (DEBUG) {
            System.out.println("Bond Fingerprints Calculation BEGIN");
        }
        /*
         * Bond formed and cleaved fingerprints
         */

        if (DEBUG) {
            System.out.println("Bond C/F Fingerprints Calculation BEGIN");
        }
        this.formedCleavedWFingerprint = new PatternFingerprinter();
        this.formedCleavedWFingerprint.setFingerprintID(mappedReaction.getID() + ":" + "Bond Cleaved and Formed");

        for (IBond bond : bondCleavedFormedChanges) {
            this.formedCleavedWFingerprint.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond), 1.0));
            IAtomContainer container = getAtomContainer(bond, allReactants);
            if (container == null) {
                pEnergy += be.getEnergies(bond);
                container = getAtomContainer(bond, allProducts);
            } else {
                rEnergy += be.getEnergies(bond);
            }
            IAtomContainer cloneReactant = container.getBuilder().newInstance(IAtomContainer.class, container);
            int chippedBondIndex = container.getBondNumber(bond);
            this.totalSmallestFragmentSize += chipTheBondCountSmallestFragmentSize(cloneReactant, chippedBondIndex);
            this.energySum += be.getEnergies(bond);

        }
        if (DEBUG) {
            System.out.println("Bond C/F Fingerprints Calculation END");
        }

        if (DEBUG) {
            System.out.println("Bond Order Fingerprints Calculation BEGIN");
        }
        /*
         * Bond Order Change fingerprints
         */
        this.orderChangesWFingerprint = new PatternFingerprinter();
        this.orderChangesWFingerprint.setFingerprintID(mappedReaction.getID() + ":" + "Bond Order Change");

        for (Map.Entry<IBond, IBond> bond : bondOrderChanges.entrySet()) {
            this.orderChangesWFingerprint.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond.getKey(), bond.getValue()), 1.0));
        }

        if (DEBUG) {
            System.out.println("Bond Order Fingerprints Calculation END");
        }

        if (DEBUG) {
            System.out.println("Bond Stereo Fingerprints Calculation BEGIN");
        }
        /*
         * Stereo centre fingerprints (Chirality (R/S), (E/Z) configurations)
         */
        this.stereoChangesWFingerprint = new PatternFingerprinter();
        this.stereoChangesWFingerprint.setFingerprintID(mappedReaction.getID() + ":" + "Bond Stereo Change");

        for (IAtom atom : atomStereoChanges) {
            this.stereoChangesWFingerprint.add(new Feature(ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom), 1.0));
        }

        if (DEBUG) {
            System.out.println("Bond Stereo Fingerprints Calculation END");
        }

        if (DEBUG) {
            System.out.println("RC Fingerprints Calculation BEGIN");
        }

        /*
         * Generate Reaction centre Fingerprint and Mol-Mol Pairs (MMP)
         */
        this.reactionCenterWFingerprint = new PatternFingerprinter();
        this.reactionCenterWFingerprint.setFingerprintID(mappedReaction.getID() + ":" + "Reaction Center");

        this.reactionMoleculeMoleculePairList = new LinkedHashSet<>();

        int circleDiameter = 3;

        for (Map.Entry<IAtom, IAtom> map : mappings.entrySet()) {

            if (map.getKey().getSymbol().equals("H")) {
                continue;
            }

            IAtomContainer atomContainerR = getAtomContainer(map.getKey(), allReactants);
            IAtomContainer atomContainerP = getAtomContainer(map.getValue(), allProducts);

            List<String> rcR = new ArrayList<>();

            if (map.getKey().getFlag(REACTIVE_CENTER)) {
                for (int i = 0; i < circleDiameter; i++) {
                    String circularSMILES = getCircularSMILES(atomContainerR, map.getKey(), i, true);
                    reactionCenterWFingerprint.add(new Feature(circularSMILES, 1.0));
                    rcR.add(circularSMILES);
                }
            }

            List<String> rcP = new ArrayList<>();

            if (map.getValue().getFlag(REACTIVE_CENTER)) {
                for (int i = 0; i < circleDiameter; i++) {
                    String circularSMILES = getCircularSMILES(atomContainerP, map.getValue(), i, true);
                    reactionCenterWFingerprint.add(new Feature(circularSMILES, 1.0));
                    rcP.add(circularSMILES);
                }
            }

            if (map.getKey().getFlag(REACTIVE_CENTER) && map.getValue().getFlag(REACTIVE_CENTER)) {
                for (int i = 0; i < circleDiameter; i++) {
                    StringBuilder mmp = new StringBuilder();
                    mmp.append(rcR.get(i)).append(">>").append(rcP.get(i));
                    reactionCenterWFingerprint.add(new Feature(mmp.toString(), 1.0));
                }
                MoleculeMoleculePair molMolPair = getMolMolPair(map.getKey(), map.getValue(), atomContainerR, atomContainerP);
                this.reactionMoleculeMoleculePairList.add(molMolPair);
            }
        }

        if (DEBUG) {
            System.out.println("RC Fingerprints Calculation END");
        }

        if (DEBUG) {
            System.out.println("Bond Fingerprints Calculation END");
        }

        setEnergyDelta(rEnergy - pEnergy);
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
        getBondCleavedReactant().entrySet().stream().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
        getBondFormedProduct().entrySet().stream().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
        getBondOrderReactant().entrySet().stream().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
        getBondOrderProduct().entrySet().stream().forEach((map) -> {
            IBond bond = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id1 = "";
            String id2 = "";
            String symbol1 = "";
            String symbol2 = "";
            if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                symbol1 = bond.getAtom(0).getSymbol();
            }
            if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
        getStereoCenterAtomsReactant().entrySet().stream().forEach((map) -> {
            IAtom atom = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id = "";
            String symbol = "";
            if (atom.getProperty(ATOM_ATOM_MAPPING) != null) {
                id = atom.getProperty(ATOM_ATOM_MAPPING);
                symbol = atom.getSymbol();
            }
            result.append(symbol).append("(").append(id).append(")").append(ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom)).append("\t").append(molID);
        });
        result.append(NEW_LINE).append(NEW_LINE);
        result.append("Stereo Change FingerPrint (Product)");
        getStereoCenterAtomsProduct().entrySet().stream().forEach((map) -> {
            IAtom atom = map.getKey();
            String molID = map.getValue();
            result.append(NEW_LINE);
            String id = "";
            String symbol = "";
            if (atom.getProperty(ATOM_ATOM_MAPPING) != null) {
                id = atom.getProperty(ATOM_ATOM_MAPPING);
                symbol = atom.getSymbol();
            }
            result.append(symbol).append("(").append(id).append(")").append(ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom)).append("\t").append(molID);
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
            for (Map.Entry<IBond, String> map : getBondCleavedReactant().entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
            for (Map.Entry<IBond, String> map : getBondFormedProduct().entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
            for (Map.Entry<IBond, String> map : getBondOrderReactant().entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
            for (Map.Entry<IBond, String> map : getBondOrderProduct().entrySet()) {
                IBond bond = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id1 = "";
                String id2 = "";
                String symbol1 = "";
                String symbol2 = "";
                if (bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id1 = bond.getAtom(0).getProperty(ATOM_ATOM_MAPPING);
                    symbol1 = bond.getAtom(0).getSymbol();
                }
                if (bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING) != null) {
                    id2 = bond.getAtom(1).getProperty(ATOM_ATOM_MAPPING);
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
            for (Map.Entry<IAtom, String> map : getStereoCenterAtomsReactant().entrySet()) {
                IAtom atom = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id = "";
                String symbol = "";
                if (atom.getProperty(ATOM_ATOM_MAPPING) != null) {
                    id = atom.getProperty(ATOM_ATOM_MAPPING);
                    symbol = atom.getSymbol();
                }
                bfw.write(symbol + "(" + id + ")" + ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom)
                        + "\t" + molID);
            }
            bfw.newLine();
            bfw.newLine();
            bfw.write("Stereo Change FingerPrint (Product)");
            for (Map.Entry<IAtom, String> map : getStereoCenterAtomsProduct().entrySet()) {
                IAtom atom = map.getKey();
                String molID = map.getValue();
                bfw.newLine();
                String id = "";
                String symbol = "";
                if (atom.getProperty(ATOM_ATOM_MAPPING) != null) {
                    id = atom.getID();
                    symbol = atom.getSymbol();
                }
                bfw.write(symbol + "(" + id + ")" + ReactionMappingUtility.getCanonicalisedAtomChangePattern(atom)
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
            Map<IAtom, IAtom> mappingsClone = new HashMap<>();
            for (IMapping mapping : compressedReaction.mappings()) {
                mappingsClone.put((IAtom) mapping.getChemObject(0), (IAtom) mapping.getChemObject(1));
            }

            for (IAtomContainer mol : compressedReaction.getReactants().atomContainers()) {
                List<IAtom> atoms = getAtoms(mol);
                if (atoms.size() > 1) {
                    atoms.stream().filter((atom) -> (atom.getSymbol().equalsIgnoreCase("H") && mappingsClone.containsKey(atom))).filter((atom) -> (atom.getProperty(BOND_CHANGE_INFORMATION) == null)).forEach((atom) -> {
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
                    atoms.stream().filter((atom) -> (atom.getSymbol().equalsIgnoreCase("H") && mappingsClone.containsValue(atom))).filter((atom) -> (atom.getProperty(BOND_CHANGE_INFORMATION) == null)).forEach((atom) -> {
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

            for (Map.Entry<IAtom, IAtom> map : mappingsClone.entrySet()) {
                if (map.getKey() != null && map.getValue() != null) {
                    map.getKey().setFlag(MAPPED, true);
                    map.getValue().setFlag(MAPPED, true);
                    compressedReaction.addMapping(new Mapping(map.getKey(), map.getValue()));
                }
            }
        } catch (CloneNotSupportedException ex) {
            getLogger(BondChangeCalculator.class.getName()).log(SEVERE, null, ex);
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
        return this.mappedReaction;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Map<IAtom, IAtom> getAtomAtomMappings() {
        return mappings;
    }

    /**
     *
     * @return @throws CDKException
     */
    @Override
    public synchronized int getTotalBondBreakingEnergy() throws CDKException {
        return this.energySum;
    }

    /**
     * @return the energyDelta
     */
    @Override
    public synchronized int getEnergyDelta() {
        return abs(energyDelta);
    }

    /**
     * @param energyDelta the energyDelta to set
     */
    private void setEnergyDelta(int energyDelta) {
        this.energyDelta = energyDelta;
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
     *
     * @return
     */
    @Override
    public Collection<IAtom> getReactionCenterSet() {
        Set<IAtom> atoms = new LinkedHashSet<>();
        mappings.entrySet().stream().map((Map.Entry<IAtom, IAtom> map) -> {
            if (map.getKey().getFlag(REACTIVE_CENTER)) {
                atoms.add(map.getKey());
            }
            return map;
        }).filter((map) -> (map.getValue().getFlag(REACTIVE_CENTER)))
                .forEach((Map.Entry<IAtom, IAtom> map) -> {
                    atoms.add(map.getValue());
                });
        return atoms;
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
        this.getReactionCentreTransformationPairs().stream().map((MoleculeMoleculePair m) -> {
            if (!uniqueRPAIRS.containsKey(m.getName().toString())) {
                LinkedList<String> l = new LinkedList<>();
                uniqueRPAIRS.put(m.getName().toString(), l);
            }
            return m;
        }).forEach((MoleculeMoleculePair m) -> {
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
    @Override
    public int getTotalSmallestFragmentSize() {
        return totalSmallestFragmentSize;
    }

    @Override
    public Map<IBond, String> getBondCleavedReactant() {
        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        Map<IBond, String> bondCleavedMap = synchronizedMap(new HashMap<IBond, String>());

        this.bondCleavedFormedChanges.stream().forEach((bond) -> {
            IAtomContainer atomContainer = getAtomContainer(bond, allReactants);
            if (atomContainer != null) {
                bondCleavedMap.put(bond, atomContainer.getProperty(ATOM_ATOM_MAPPING));
            }
        });
        return bondCleavedMap;
    }

    @Override
    public Map<IBond, String> getBondFormedProduct() {
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        Map<IBond, String> bondFormedMap = synchronizedMap(new HashMap<IBond, String>());

        this.bondCleavedFormedChanges.stream().forEach((bond) -> {
            IAtomContainer atomContainer = getAtomContainer(bond, allProducts);
            if (atomContainer != null) {
                bondFormedMap.put(bond, atomContainer.getID());
            }
        });
        return bondFormedMap;
    }

    @Override
    public Map<IBond, String> getBondOrderProduct() {
        Map<IBond, String> bondOrderPMap = synchronizedMap(new HashMap<IBond, String>());
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        this.bondOrderChanges.entrySet().stream().forEach((map) -> {
            IAtomContainer atomContainer = getAtomContainer(map.getValue(), allProducts);
            bondOrderPMap.put(map.getValue(), atomContainer.getID());
        });

        return bondOrderPMap;
    }

    @Override
    public Map<IBond, String> getBondOrderReactant() {
        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        Map<IBond, String> bondOrderRMap = synchronizedMap(new HashMap<IBond, String>());
        this.bondOrderChanges.entrySet().stream().forEach((map) -> {
            IAtomContainer atomContainer = getAtomContainer(map.getKey(), allReactants);
            bondOrderRMap.put(map.getKey(), atomContainer.getID());
        });
        return bondOrderRMap;
    }

    @Override
    public Map<IAtom, String> getStereoCenterAtomsProduct() {
        Map<IAtom, String> atomStereoPMap = synchronizedMap(new HashMap<IAtom, String>());
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        this.atomStereoChanges.stream().forEach((atom) -> {
            IAtomContainer atomContainer = getAtomContainer(atom, allProducts);
            if (atomContainer != null) {
                atomStereoPMap.put(atom, atomContainer.getID());
            }
        });
        return atomStereoPMap;
    }

    @Override
    public Map<IAtom, String> getStereoCenterAtomsReactant() {
        Map<IAtom, String> atomStereoRMap = synchronizedMap(new HashMap<IAtom, String>());
        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        this.atomStereoChanges.stream().forEach((atom) -> {
            IAtomContainer atomContainer = getAtomContainer(atom, allReactants);
            if (atomContainer != null) {
                atomStereoRMap.put(atom, atomContainer.getID());
            }
        });

        return atomStereoRMap;
    }
}
