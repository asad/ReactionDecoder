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

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import static java.util.Collections.synchronizedMap;
import java.util.List;
import java.util.Map;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.Bond;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import static org.openscience.cdk.CDKConstants.REACTIVE_CENTER;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomStereoChangeInformation;
import uk.ac.ebi.reactionblast.mechanism.helper.BondChange;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.BOND_CHANGE_INFORMATION;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.E;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.R;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.S;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.Z;
import static java.lang.Math.abs;
import static java.lang.System.getProperty;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public final class BondChangeAnnotator extends DUModel {

    static final String NEW_LINE = getProperty("line.separator");
    private static final long serialVersionUID = 988987678877861L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(BondChangeAnnotator.class);
    private static final boolean DEBUG = false;

    /**
     *
     * @param reaction with Atom-Atom Mapping produced by Reactor class
     * @param withoutHydrogen
     * @param generate2D
     * @param generate3D
     * @throws Exception
     */
    protected BondChangeAnnotator(IReaction reaction,
            boolean withoutHydrogen,
            boolean generate2D,
            boolean generate3D) throws Exception {
        super(reaction, withoutHydrogen, generate2D, generate3D);
        if (DEBUG) {
            System.out.println("MARK Bond Change START");
        }
        markBondChanges();
        if (DEBUG) {
            System.out.println("MARK Bond Change END");
        }
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMappingContainer getMappingContainer() {
        return mapping;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized BEMatrix getEductBEMatrix() {
        return reactantBE;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized List<BondChange> getBondChangeList() {
        return bondChangeList;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Collection<IAtom> getReactionCenterSet() {
        return reactionCenterList;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized List<AtomStereoChangeInformation> getStereoChangeList() {
        return stereoChangeList;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized BEMatrix getProductBEMatrix() {
        return productBE;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Map<IAtom, IAtom> getMappingMap() {
        return synchronizedMap(mappingMap);
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized RMatrix getRMatrix() {
        return reactionMatrix;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized boolean hasRMatrix() {
        return reactionMatrix != null;
    }

    /**
     *
     */
    @Override
    public synchronized void printBMatrix() {
        printBEMatrix(reactantBE);
    }

    /**
     *
     * @return
     */
    @Override
    public List<AtomStereoChangeInformation> getConformationChangeList() {
        return conformationChangeList;
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeBMatrix(File outputFile) {
        try {
            writeBEMatrix(outputFile, reactantBE);
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     *
     */
    @Override
    public synchronized void printEMatrix() {
        printBEMatrix(productBE);
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeEMatrix(File outputFile) {
        try {
            writeBEMatrix(outputFile, productBE);
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     *
     */
    @Override
    public synchronized void printRMatrix() {
        printReactionMatrix(reactionMatrix);
    }

    /**
     *
     * @param outputFile
     */
    @Override
    public synchronized void writeRMatrix(File outputFile) {
        try {
            writeReactionMatrix(outputFile, reactionMatrix);
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     * set bond change flag
     *
     * @throws Exception
     */
    protected synchronized void markBondChanges() throws Exception {

        BEMatrix substrateBEMatrix = reactantBE;
        BEMatrix productBEMatrix = productBE;

        if (DEBUG) {
            System.out.println("markBondChanges method start");
//        System.out.println(reactantBE.toString());
//        System.out.println(productBE.toString());
//        System.out.println(reactionMatrix.toString
        }

        /*
         * Marking CDKConstants.ISINRING FLAGS
         */
        if (DEBUG) {
            System.out.println("Marking Rings");
        }
        for (IAtomContainer atomContainerQ : reactantSet.atomContainers()) {
            try {
                /*
                 * set Flag(CDKConstants.ISINRING)
                 */
                initializeMolecule(atomContainerQ);
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
//            IRingSet singleRingsQ = new SSSRFinder(atomContainerQ).findSSSR();
            //New Method
            CycleFinder cf = Cycles.mcb();
            Cycles cycles = cf.find(atomContainerQ); // ignore error - essential cycles do not check tractability
            IRingSet singleRingsQ = cycles.toRingSet();
            queryRingSet.add(singleRingsQ);
        }

        for (IAtomContainer atomContainerT : productSet.atomContainers()) {
            try {
                /*
                 * set Flag(CDKConstants.ISINRING)
                 */
                initializeMolecule(atomContainerT);
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
//            IRingSet singleRingsT = new SSSRFinder(atomContainerT).findSSSR();
            //New Method
            CycleFinder cf = Cycles.mcb();
            Cycles cycles = cf.find(atomContainerT); // ignore error - essential cycles do not check tractability
            IRingSet singleRingsT = cycles.toRingSet();
            targetRingSet.add(singleRingsT);
        }
        /*
         * Mining Stereo Atom Changes E/Z or R/S only
         */

        if (DEBUG) {
            System.out.println("Marking E/Z or R/S");
        }
        for (StereoChange sc : stereogenicCenters) {
            IAtom atomE = sc.getReactantAtom();
            IAtom atomP = sc.getProductAtom();

            IStereoAndConformation rsb = sc.getReactantAtomStereo();
            IStereoAndConformation psb = sc.getProductAtomStereo();

            if (atomE != null && atomP != null) {
                if (atomE.getSymbol().equals("P") || atomP.getSymbol().equals("P")) {
                    LOGGER.warn(NEW_LINE + "WARNING: The stereo change " + atomE.getSymbol()
                            + " not supported");
                    continue;
                }
                atomE.setProperty(ATOM_STEREO_CHANGE_INFORMATION, rsb);
                atomP.setProperty(ATOM_STEREO_CHANGE_INFORMATION, psb);
                atomE.setFlag(REACTIVE_CENTER, true);
                atomP.setFlag(REACTIVE_CENTER, true);
                getReactionCenterSet().add(atomE);
                getReactionCenterSet().add(atomP);

                if ((sc.getReactantAtomStereo().equals(E)
                        || sc.getProductAtomStereo().equals(Z))
                        || (sc.getReactantAtomStereo().equals(Z)
                        || sc.getProductAtomStereo().equals(E))) {
                    getConformationChangeList().add(new AtomStereoChangeInformation(atomE, atomP, sc.getReactantAtomStereo(), sc.getProductAtomStereo()));
                } else if ((sc.getReactantAtomStereo().equals(R)
                        || sc.getProductAtomStereo().equals(S))
                        || (sc.getReactantAtomStereo().equals(S)
                        || sc.getProductAtomStereo().equals(R))) {
                    getStereoChangeList().add(new AtomStereoChangeInformation(atomE, atomP, sc.getReactantAtomStereo(), sc.getProductAtomStereo()));
                }
            }
        }
        /*
         * Mining bond cleaved formed, order change except stereo information
         *
         */

        if (DEBUG) {
            System.out.println("Marking Bond Changes");
        }
        int sizeQ = reactionMatrix.getReactantsAtomArray().size();
        int sizeT = reactionMatrix.getProductsAtomArray().size();

        for (int i = 0; i < reactionMatrix.getRowDimension(); i++) {
            for (int j = i; j < reactionMatrix.getColumnDimension(); j++) {
                if (DEBUG) {
                    System.out.println("Marking Bond Changes-1");
                }
                if (i != j && reactionMatrix.getValue(i, j) == 0.) {
                    IBond affectedBondReactants = null;
                    IBond affectedBondProducts = null;
                    ECBLAST_BOND_CHANGE_FLAGS bondChangeInformation;
                    try {
                        if (i < sizeQ && j < sizeQ) {
                            affectedBondReactants = getBondOfReactantsByRMatrix(reactionMatrix.getReactantAtom(i), reactionMatrix.getReactantAtom(j));
                        }
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                    try {
                        if (i < sizeT && j < sizeT) {
                            affectedBondProducts = getBondOfProductsByRMatrix(reactionMatrix.getProductAtom(i), reactionMatrix.getProductAtom(j));
                        }
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                    if (affectedBondReactants == null && affectedBondProducts == null) {
                        continue;
                    }

                    int kekuleEffect = 0;
                    kekuleEffect = isAlternateKekuleChange(affectedBondReactants, affectedBondProducts);
                    if (kekuleEffect == 0) {
                        bondChangeInformation = BOND_ORDER;
                        if (affectedBondReactants != null) {
                            affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                            affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);
                            getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                            getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                            affectedBondReactants.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                        }
                        if (affectedBondProducts != null) {
                            affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                            affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);
                            getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                            getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                            affectedBondProducts.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                        }
                        getBondChangeList().add(new BondChange(affectedBondReactants, affectedBondProducts));
                    }
                }

                /*
                 * R-Matrix with changes
                 */
                if (DEBUG) {
                    System.out.println("Marking Bond Changes-2");
                }
                if (reactionMatrix.getValue(i, j) != 0.) {

                    /*
                     * DEBUG
                     */
                    if (DEBUG) {
                        System.out.println("Bond Change in R Matrix " + " i "
                                + (i + 1) + ", j " + (j + 1) + " " + reactionMatrix.getValue(i, j));
                    }

                    //Diagonal free valence electron changes 
                    if (i == j) {
                        IAtom reactantAtom;
                        IAtom productAtom;
                        try {
                            reactantAtom = reactionMatrix.getReactantAtom(i);
                            if (reactantAtom != null) {
                                reactantAtom.setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(reactantAtom);
                            }
                        } catch (CDKException ex) {
                            LOGGER.error(SEVERE, null, ex);
                        }
                        try {
                            productAtom = reactionMatrix.getProductAtom(j);
                            if (productAtom != null) {
                                productAtom.setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(productAtom);
                            }
                        } catch (CDKException ex) {
                            LOGGER.error(SEVERE, null, ex);
                        }
                    }

                    /*
                     * off diagonal changes
                     */
                    IBond affectedBondReactants;
                    IBond affectedBondProducts;
                    ECBLAST_BOND_CHANGE_FLAGS bondChangeInformation;
                    if (DEBUG) {
                        System.out.println("Marking Bond Changes-2");
                    }
                    try {
                        affectedBondReactants = getBondOfReactantsByRMatrix(reactionMatrix.getReactantAtom(i), reactionMatrix.getReactantAtom(j));
                        affectedBondProducts = getBondOfProductsByRMatrix(reactionMatrix.getProductAtom(i), reactionMatrix.getProductAtom(j));
                        if (affectedBondReactants == null && affectedBondProducts == null) {
                            continue;
                        }
                        if (affectedBondReactants != null
                                && affectedBondProducts != null
                                && affectedBondReactants.getProperties().containsKey(BOND_CHANGE_INFORMATION)
                                && affectedBondProducts.getProperties().containsKey(BOND_CHANGE_INFORMATION)) {
                            continue;
                        }
                        int kekuleEffect = isKekuleEffect(affectedBondReactants, affectedBondProducts);
                        if (kekuleEffect == 1) {
                            continue;
                        }

                        if (DEBUG) {
                            System.out.println(i + "," + j + " reactionMatrix.getValue(i, j) " + reactionMatrix.getValue(i, j));
                        }

                        /*
                         * Changes in the product
                         */
                        if (reactionMatrix.getValue(i, j) < 0.0d) {
                            if (DEBUG) {
                                System.out.println("Marking Bond Changes-2 product");
                            }

                            if (productBEMatrix.getValue(i, j) == 0.0d && affectedBondProducts == null) {
                                /*
                                 * Here the bond is cleaved (Reduced)
                                 */
                                bondChangeInformation = BOND_CLEAVED;
                            } else {
                                bondChangeInformation = BOND_ORDER;
                            }

                            if (affectedBondReactants != null) {

                                affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                                affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                                getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                                affectedBondReactants.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                            if (affectedBondProducts != null) {

                                affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                                affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                                getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                                affectedBondProducts.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                        } /*
                           * Changes in the educt
                         */ else if (reactionMatrix.getValue(i, j) > 0.d) {

                            if (DEBUG) {
                                System.out.println("Marking Bond Changes-2 educt");
                            }

                            if (substrateBEMatrix.getValue(i, j) == 0.0d && affectedBondReactants == null) {
                                /*
                                 * Here the bond is Formed (Gained)
                                 */
                                bondChangeInformation = BOND_FORMED;
                            } else {
                                bondChangeInformation = BOND_ORDER;
                            }

                            if (affectedBondReactants != null) {

                                affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                                affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                                getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                                affectedBondReactants.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                            if (affectedBondProducts != null) {

                                affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                                affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                                getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                                affectedBondProducts.setProperty(BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                        }
                        /*
                         * Store the bond changes
                         */
                        if (DEBUG) {
                            System.out.println("Marking Bond Changes-2 STORED ");
                        }

                        getBondChangeList().add(new BondChange(affectedBondReactants, affectedBondProducts));
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                }
            }
        }
        if (DEBUG) {
            System.out.println("Marking Bond Changes-DONE");
        }
        /*
         * Marking Missing Bond Changes
         */
        markHydrogenDisplacementBondChanges();
        /*
         * Marking Un-Mapped Atoms
         */
        markUnMappedAtoms();

        if (DEBUG) {
            System.out.println("markBondChanges method END");
        }
    }

    private synchronized void markHydrogenDisplacementBondChanges() {

        if (DEBUG) {
            System.out.println("markHydrogenDisplacementBondChanges method START");
        }
        /*
         * Mark Hydrogen bond broken/Formed in the reaction
         *
         */
        for (IAtom eductAtom : mappingMap.keySet()) {
            IAtom productAtom = mappingMap.get(eductAtom);

            IBond affectedBondReactants = null;
            IBond affectedBondProducts = null;

            IAtomContainer rMol = getAtomContainer(eductAtom, reactantSet);
            IAtomContainer pMol = getAtomContainer(productAtom, productSet);

            if (eductAtom.getSymbol().equalsIgnoreCase("H")
                    && productAtom.getSymbol().equalsIgnoreCase("H")) {
                List<IBond> connectedEductBondsList = rMol.getConnectedBondsList(eductAtom);
                List<IBond> connectedProductBondsList = pMol.getConnectedBondsList(productAtom);

                /*
                 * Look for cases like O + CH4 -> H-C-H + H2O
                 */
                if (rMol.getConnectedBondsCount(eductAtom) > 0
                        && pMol.getConnectedBondsCount(productAtom) > 0) {

                    for (IBond eBond : connectedEductBondsList) {
                        boolean isBondChange = true;
                        String attachedEAtomID = eBond.getOther(eductAtom).getID();

                        for (IBond pBond : connectedProductBondsList) {
                            String attachedPAtomID = pBond.getOther(productAtom).getID();

                            if (attachedEAtomID.equalsIgnoreCase(attachedPAtomID)) {
                                isBondChange = false;
                                break;
                            }
                        }
                        if (isBondChange) {
                            eBond.setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
                            affectedBondReactants = eBond;
                        }
                    }

                    for (IBond pBond : connectedProductBondsList) {
                        boolean isBondChange = true;
                        String attachedPAtomID = pBond.getOther(productAtom).getID();
                        for (IBond eBond : connectedEductBondsList) {
                            String attachedEAtomID = eBond.getOther(eductAtom).getID();
                            if (attachedPAtomID.equalsIgnoreCase(attachedEAtomID)) {
                                isBondChange = false;
                                break;
                            }
                        }
                        if (isBondChange) {
                            pBond.setProperty(BOND_CHANGE_INFORMATION, BOND_FORMED);
                            affectedBondProducts = pBond;
                        }
                    }
                    if (affectedBondReactants != null && affectedBondProducts != null) {
                        affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                        affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);
                        affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                        affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);

                        getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                        getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                        getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                        getReactionCenterSet().add(affectedBondProducts.getAtom(1));

                        BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
                        getBondChangeList().add(bondChange);
                    }
                } else if (rMol.getConnectedBondsCount(eductAtom) == 0
                        && pMol.getConnectedBondsCount(productAtom) > 0) {

                    IAtom psudoAtom = new PseudoAtom("PsH");
                    IBond eBond = new Bond(psudoAtom, eductAtom, SINGLE);
                    IBond pBond = connectedProductBondsList.iterator().next();

                    eBond.setProperty(BOND_CHANGE_INFORMATION, PSEUDO_BOND);
                    pBond.setProperty(BOND_CHANGE_INFORMATION, BOND_FORMED);
                    affectedBondReactants = eBond;
                    affectedBondProducts = pBond;

                    affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);

                    getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                    getReactionCenterSet().add(affectedBondProducts.getAtom(1));

                    BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
//                        BondChange bondChange = new BondChange(null, affectedBondProducts);
                    getBondChangeList().add(bondChange);

                } else if (rMol.getConnectedBondsCount(eductAtom) > 0
                        && pMol.getConnectedBondsCount(productAtom) == 0) {

                    IBond eBond = connectedEductBondsList.iterator().next();

                    IAtom psudoAtom = new PseudoAtom("PsH");
                    IBond pBond = new Bond(psudoAtom, productAtom, SINGLE);

                    eBond.setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
                    pBond.setProperty(BOND_CHANGE_INFORMATION, PSEUDO_BOND);

                    affectedBondReactants = eBond;
                    affectedBondProducts = pBond;

                    affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);

                    getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                    getReactionCenterSet().add(affectedBondReactants.getAtom(1));

                    BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
//                        BondChange bondChange = new BondChange(affectedBondReactants, null);
                    getBondChangeList().add(bondChange);

                }
            }
        }
        if (DEBUG) {
            System.out.println("markHydrogenDisplacementBondChanges method END");
        }

    }

    private synchronized void markUnMappedAtoms() {
        if (DEBUG) {
            System.out.println("markUnMappedAtoms method START");
        }
        for (IAtomContainer acE : reactantSet.atomContainers()) {
            for (IBond affectedBondReactants : acE.bonds()) {
                boolean isNotMapped = false;
                if ((!mappingMap.containsKey(affectedBondReactants.getAtom(0))
                        && mappingMap.containsKey(affectedBondReactants.getAtom(1)))
                        || (mappingMap.containsKey(affectedBondReactants.getAtom(0))
                        && !mappingMap.containsKey(affectedBondReactants.getAtom(1)))) {
                    isNotMapped = true;
                }
                IBond pBond = null;
                if (isNotMapped) {

                    if (DEBUG) {
                        System.out.println("affectedBondReactants-0 " + affectedBondReactants.getAtom(0).getID());
                        System.out.println("affectedBondReactants-1 " + affectedBondReactants.getAtom(1).getID());
                    }
                    affectedBondReactants.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    affectedBondReactants.getAtom(1).setFlag(REACTIVE_CENTER, true);
                    getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                    getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                    affectedBondReactants.setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
                    BondChange bondChange = new BondChange(affectedBondReactants, pBond);
                    getBondChangeList().add(bondChange);
                }
            }
        }

        for (IAtomContainer acP : productSet.atomContainers()) {
            for (IBond affectedBondProducts : acP.bonds()) {
                boolean isNotMapped = false;
                if ((!mappingMap.containsValue(affectedBondProducts.getAtom(0))
                        && mappingMap.containsValue(affectedBondProducts.getAtom(1)))
                        || (mappingMap.containsValue(affectedBondProducts.getAtom(0))
                        && !mappingMap.containsValue(affectedBondProducts.getAtom(1)))) {
                    isNotMapped = true;
                }
                IBond eBond = null;
                if (isNotMapped) {
                    affectedBondProducts.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    affectedBondProducts.getAtom(1).setFlag(REACTIVE_CENTER, true);
                    getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                    getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                    affectedBondProducts.setProperty(BOND_CHANGE_INFORMATION, BOND_FORMED);
                    BondChange bondChange = new BondChange(eBond, affectedBondProducts);
                    getBondChangeList().add(bondChange);
                }
            }
        }
        if (DEBUG) {
            System.out.println("markUnMappedAtoms method END");
        }
    }

    private synchronized IBond getBondOfReactantsByRMatrix(IAtom atom1, IAtom atom2) {
        for (int i = 0; i < reactantSet.getAtomContainerCount(); i++) {
            if (reactantSet.getAtomContainer(i).getBond(atom1, atom2) != null) {
                return reactantSet.getAtomContainer(i).getBond(atom1, atom2);
            }
        }
        return null;
    }

    private synchronized IBond getBondOfProductsByRMatrix(IAtom atom1, IAtom atom2) {
        for (int i = 0; i < productSet.getAtomContainerCount(); i++) {
            if (productSet.getAtomContainer(i).getBond(atom1, atom2) != null) {
                return productSet.getAtomContainer(i).getBond(atom1, atom2);
            }
        }
        return null;
    }

    /**
     *
     * @param affectedBondReactants
     * @param affectedBondProducts
     * @return
     */
    public int isKekuleEffect(IBond affectedBondReactants, IBond affectedBondProducts) {
        if (affectedBondReactants != null && affectedBondProducts != null) {
            if (affectedBondReactants.getFlag(ISINRING)
                    == affectedBondProducts.getFlag(ISINRING)) {

                if ((!affectedBondReactants.getFlag(ISAROMATIC)
                        && affectedBondProducts.getFlag(ISAROMATIC))
                        || (affectedBondReactants.getFlag(ISAROMATIC)
                        && !affectedBondProducts.getFlag(ISAROMATIC))) {
                    IRingSet smallestRingSetR = getSmallestRingSet(affectedBondReactants, queryRingSet);
                    IRingSet smallestRingSetP = getSmallestRingSet(affectedBondProducts, targetRingSet);

                    int changeCounter = 0;
                    for (IAtomContainer rRing : smallestRingSetR.atomContainers()) {
                        for (IAtomContainer pRing : smallestRingSetP.atomContainers()) {
                            if (rRing.getAtomCount() == pRing.getAtomCount()) {
                                for (IBond rBond : rRing.bonds()) {
                                    for (IBond pBond : pRing.bonds()) {
                                        if (rBond.getAtom(0).getID().equals(pBond.getAtom(0).getID())
                                                && rBond.getAtom(1).getID().equals(pBond.getAtom(1).getID())) {
                                            if (!rBond.getOrder().equals(pBond.getOrder())) {
                                                changeCounter++;
                                            }
                                        }
                                        if (rBond.getAtom(0).getID().equals(pBond.getAtom(1).getID())
                                                && rBond.getAtom(1).getID().equals(pBond.getAtom(0).getID())) {
                                            if (!rBond.getOrder().equals(pBond.getOrder())) {
                                                changeCounter++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (changeCounter == 1) {
                        return 0;
                    } else {
                        return 1;
                    }
                }
            }
        }

        return -1;
    }

    /**
     *
     * @param affectedBondReactants
     * @param affectedBondProducts
     * @return
     */
    public int isAlternateKekuleChange(IBond affectedBondReactants, IBond affectedBondProducts) {
        if (affectedBondReactants != null && affectedBondProducts != null) {
            if (affectedBondReactants.getFlag(ISINRING)
                    == affectedBondProducts.getFlag(ISINRING)) {

                if ((!affectedBondReactants.getFlag(ISAROMATIC)
                        && affectedBondProducts.getFlag(ISAROMATIC))
                        || (affectedBondReactants.getFlag(ISAROMATIC)
                        && !affectedBondProducts.getFlag(ISAROMATIC))) {
                    IRingSet smallestRingSetR = getSmallestRingSet(affectedBondReactants, queryRingSet);
                    IRingSet smallestRingSetP = getSmallestRingSet(affectedBondProducts, targetRingSet);
                    int countR = getNeighbourBondOrderCountFromRing(affectedBondReactants, smallestRingSetR);
                    int countP = getNeighbourBondOrderCountFromRing(affectedBondProducts, smallestRingSetP);
                    int sizeR = smallestRingSetR.getAtomContainerCount() > 0 ? smallestRingSetR.atomContainers().iterator().next().getBondCount() : 0;
                    int sizeP = smallestRingSetP.getAtomContainerCount() > 0 ? smallestRingSetP.atomContainers().iterator().next().getBondCount() : 0;

                    IAtom atomR1 = affectedBondReactants.getAtom(0);
                    IAtom atomR2 = affectedBondReactants.getAtom(1);

                    IAtom atomP1 = affectedBondProducts.getAtom(0);
                    IAtom atomP2 = affectedBondProducts.getAtom(1);

                    if (atomR1.getID().equals(atomP1.getID())
                            && atomR2.getID().equals(atomP2.getID())) {
                        if (atomR1.getHybridization() != null
                                && atomR2.getHybridization() != null
                                && atomP1.getHybridization() != null
                                && atomP2.getHybridization() != null) {

                            if (!atomR1.getHybridization().equals(atomP1.getHybridization())
                                    && !atomR2.getHybridization().equals(atomP2.getHybridization())) {
                                if (countR == countP && sizeR == sizeP) {
                                    return 1;
                                } else {
                                    if (!affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                        if (abs(countR - countP) <= 1 && sizeR == sizeP) {
                                            return 1;
                                        }
                                    }
                                    if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                        if (abs(countR - countP) == 0 && sizeR == sizeP) {
                                            return 0;
                                        }
                                    }
                                }
                            } else if (atomR1.getHybridization().equals(atomP1.getHybridization())
                                    && atomR2.getHybridization().equals(atomP2.getHybridization())) {

                                if (!affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (abs(countR - countP) <= 1
                                            && sizeR == sizeP) {
                                        return 0;
                                    }
                                }
                                if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (abs(countR - countP) > 1 && sizeR == sizeP) {
                                        return 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return -1;
    }
}
