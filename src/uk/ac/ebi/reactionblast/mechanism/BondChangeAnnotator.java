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

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.ringsearch.SSSRFinder;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomStereoChangeInformation;
import uk.ac.ebi.reactionblast.mechanism.helper.BondChange;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public final class BondChangeAnnotator extends DUModel {

    private static final long serialVersionUID = 988987678877861L;

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
        markBondChanges();
    }

    @Override
    public synchronized AtomAtomMappingContainer getMappingContainer() {
        return mapping;
    }

    @Override
    public synchronized BEMatrix getEductBEMatrix() {
        return reactantBE;
    }

    @Override
    public synchronized List<BondChange> getBondChangeList() {
        return bondChangeList;
    }

    @Override
    public synchronized Collection<IAtom> getReactionCenterSet() {
        return reactionCenterList;
    }

    @Override
    public synchronized List<AtomStereoChangeInformation> getStereoChangeList() {
        return stereoChangeList;
    }

    @Override
    public synchronized BEMatrix getProductBEMatrix() {
        return productBE;
    }

    @Override
    public synchronized Map<IAtom, IAtom> getMappingMap() {
        return Collections.synchronizedMap(mappingMap);
    }

    @Override
    public synchronized RMatrix getRMatrix() {
        return reactionMatrix;
    }

    @Override
    public synchronized boolean hasRMatrix() {
        return reactionMatrix != null;
    }

    @Override
    public synchronized void printBMatrix() {
        printBEMatrix(reactantBE);
    }

    @Override
    public List<AtomStereoChangeInformation> getConformationChangeList() {
        return conformationChangeList;
    }

    @Override
    public synchronized void writeBMatrix(File outputFile) {
        try {
            writeBEMatrix(outputFile, reactantBE);
        } catch (IOException ex) {
            Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public synchronized void printEMatrix() {
        printBEMatrix(productBE);
    }

    @Override
    public synchronized void writeEMatrix(File outputFile) {
        try {
            writeBEMatrix(outputFile, productBE);
        } catch (IOException ex) {
            Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public synchronized void printRMatrix() {
        printReactionMatrix(reactionMatrix);
    }

    @Override
    public synchronized void writeRMatrix(File outputFile) {
        try {
            writeReactionMatrix(outputFile, reactionMatrix);
        } catch (IOException ex) {
            Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
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

//        System.out.println(reactantBE.toString());
//        System.out.println(productBE.toString());
//        System.out.println(this.reactionMatrix.toString());
        for (IAtomContainer atomContainerQ : reactantSet.atomContainers()) {
            try {
                /*
                 * set Flag(CDKConstants.ISINRING)
                 */
                initializeMolecule(atomContainerQ);
            } catch (CDKException ex) {
                Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
            }
            IRingSet singleRingsQ = new SSSRFinder(atomContainerQ).findSSSR();
            queryRingSet.add(singleRingsQ);
        }

        for (IAtomContainer atomContainerT : productSet.atomContainers()) {
            try {
                /*
                 * set Flag(CDKConstants.ISINRING)
                 */
                initializeMolecule(atomContainerT);
            } catch (CDKException ex) {
                Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
            }
            IRingSet singleRingsT = new SSSRFinder(atomContainerT).findSSSR();
            targetRingSet.add(singleRingsT);
        }
        /*
         * Mining Stereo Atom Changes E/Z or R/S only
         */

        for (StereoChange sc : stereogenicCenters) {
            IAtom atomE = sc.getReactantAtom();
            IAtom atomP = sc.getProductAtom();

            IStereoAndConformation rsb = sc.getReactantAtomStereo();
            IStereoAndConformation psb = sc.getProductAtomStereo();

            if (atomE != null && atomP != null) {
                if (atomE.getSymbol().equals("P") || atomP.getSymbol().equals("P")) {
                    System.out.println("\nWARNING: The stereo change " + atomE.getSymbol()
                            + " not supported");
                    continue;
                }
                atomE.setProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION, rsb);
                atomP.setProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION, psb);
                atomE.setFlag(CDKConstants.REACTIVE_CENTER, true);
                atomP.setFlag(CDKConstants.REACTIVE_CENTER, true);
                getReactionCenterSet().add(atomE);
                getReactionCenterSet().add(atomP);

                if ((sc.getReactantAtomStereo().equals(IStereoAndConformation.E)
                        || sc.getProductAtomStereo().equals(IStereoAndConformation.Z))
                        || (sc.getReactantAtomStereo().equals(IStereoAndConformation.Z)
                        || sc.getProductAtomStereo().equals(IStereoAndConformation.E))) {
                    getConformationChangeList().add(new AtomStereoChangeInformation(atomE, atomP, sc.getReactantAtomStereo(), sc.getProductAtomStereo()));
                } else if ((sc.getReactantAtomStereo().equals(IStereoAndConformation.R)
                        || sc.getProductAtomStereo().equals(IStereoAndConformation.S))
                        || (sc.getReactantAtomStereo().equals(IStereoAndConformation.S)
                        || sc.getProductAtomStereo().equals(IStereoAndConformation.R))) {
                    getStereoChangeList().add(new AtomStereoChangeInformation(atomE, atomP, sc.getReactantAtomStereo(), sc.getProductAtomStereo()));
                }
            }
        }
        /*
         * Mining bond cleaved formed, order change except stereo information
         *
         */
        int sizeQ = reactionMatrix.getReactantsAtomArray().size();
        int sizeT = reactionMatrix.getProductsAtomArray().size();

        for (int i = 0; i < reactionMatrix.getRowDimension(); i++) {
            for (int j = i; j < reactionMatrix.getColumnDimension(); j++) {
                if (i != j && reactionMatrix.getValue(i, j) == 0.) {
                    IBond affectedBondReactants = null;
                    IBond affectedBondProducts = null;
                    ECBLAST_BOND_CHANGE_FLAGS bondChangeInformation;
                    try {
                        if (i < sizeQ && j < sizeQ) {
                            affectedBondReactants = getBondOfReactantsByRMatrix(reactionMatrix.getReactantAtom(i), reactionMatrix.getReactantAtom(j));
                        }
                    } catch (CDKException ex) {
                        Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    try {
                        if (i < sizeT && j < sizeT) {
                            affectedBondProducts = getBondOfProductsByRMatrix(reactionMatrix.getProductAtom(i), reactionMatrix.getProductAtom(j));
                        }
                    } catch (CDKException ex) {
                        Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    if (affectedBondReactants == null && affectedBondProducts == null) {
                        continue;
                    }

                    int kekuleEffect = isAlternateKekuleChange(affectedBondReactants, affectedBondProducts);
                    if (kekuleEffect == 0) {
                        bondChangeInformation = ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER;
                        if (affectedBondReactants != null) {
                            affectedBondReactants.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                            affectedBondReactants.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                            getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                            getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                            affectedBondReactants.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                        }
                        if (affectedBondProducts != null) {
                            affectedBondProducts.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                            affectedBondProducts.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                            getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                            getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                            affectedBondProducts.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                        }
                        getBondChangeList().add(new BondChange(affectedBondReactants, affectedBondProducts));
                    }
                }

                /*
                 * R-Matrix with changes
                 */
                if (reactionMatrix.getValue(i, j) != 0.) {

                    //Diagonal free valence electron changes 
                    if (i == j) {
                        IAtom reactantAtom;
                        IAtom productAtom;
                        try {
                            reactantAtom = reactionMatrix.getReactantAtom(i);
                            if (reactantAtom != null) {
                                reactantAtom.setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(reactantAtom);
                            }
                        } catch (CDKException ex) {
                            Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        try {
                            productAtom = reactionMatrix.getProductAtom(j);
                            if (productAtom != null) {
                                productAtom.setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(productAtom);
                            }
                        } catch (CDKException ex) {
                            Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }

                    /*
                     * off diagonal changes
                     */
                    IBond affectedBondReactants;
                    IBond affectedBondProducts;
                    ECBLAST_BOND_CHANGE_FLAGS bondChangeInformation;

                    try {
                        affectedBondReactants = getBondOfReactantsByRMatrix(reactionMatrix.getReactantAtom(i), reactionMatrix.getReactantAtom(j));
                        affectedBondProducts = getBondOfProductsByRMatrix(reactionMatrix.getProductAtom(i), reactionMatrix.getProductAtom(j));
                        if (affectedBondReactants == null && affectedBondProducts == null) {
                            continue;
                        }
                        if (affectedBondReactants != null
                                && affectedBondProducts != null
                                && affectedBondReactants.getProperties().containsKey(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION)
                                && affectedBondProducts.getProperties().containsKey(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION)) {
                            continue;
                        }
                        int kekuleEffect = isKekuleEffect(affectedBondReactants, affectedBondProducts);
                        if (kekuleEffect == 1) {
                            continue;
                        }

                        /*
                         * Changes in the product
                         */
                        if (reactionMatrix.getValue(i, j) < 0.0d) {
                            if (productBEMatrix.getValue(i, j) == 0.0d
                                    && affectedBondProducts == null) {
                                /*
                                 * Here the bond is cleaved (Reduced)
                                 */
                                bondChangeInformation = ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED;
                            } else {
                                bondChangeInformation = ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER;
                            }

                            if (affectedBondReactants != null) {

                                affectedBondReactants.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                affectedBondReactants.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                                getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                                affectedBondReactants.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                            if (affectedBondProducts != null) {

                                affectedBondProducts.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                affectedBondProducts.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                                getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                                affectedBondProducts.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                        } /*
                         * Changes in the educt
                         */ else if (reactionMatrix.getValue(i, j) > 0.d) {

                            if (substrateBEMatrix.getValue(i, j) == 0.0d
                                    && affectedBondReactants == null) {
                                /*
                                 * Here the bond is Formed (Gained)
                                 */
                                bondChangeInformation = ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED;
                            } else {

                                bondChangeInformation = ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER;
                            }

                            if (affectedBondReactants != null) {

                                affectedBondReactants.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                affectedBondReactants.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                                getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                                affectedBondReactants.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                            if (affectedBondProducts != null) {

                                affectedBondProducts.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                affectedBondProducts.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                                getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                                getReactionCenterSet().add(affectedBondProducts.getAtom(1));
                                affectedBondProducts.getProperties().put(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION, bondChangeInformation);
                            }
                        }
                        /*
                         * Store the bond changes
                         */
                        getBondChangeList().add(new BondChange(affectedBondReactants, affectedBondProducts));
                    } catch (CDKException ex) {
                        Logger.getLogger(BondChangeAnnotator.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        }
        markHydrogenDisplacementBondChanges();
        markUnMappedAtoms();
    }

    private synchronized void markHydrogenDisplacementBondChanges() {

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
                        String attachedEAtomID = eBond.getConnectedAtom(eductAtom).getID();

                        for (IBond pBond : connectedProductBondsList) {
                            String attachedPAtomID = pBond.getConnectedAtom(productAtom).getID();

                            if (attachedEAtomID.equalsIgnoreCase(attachedPAtomID)) {
                                isBondChange = false;
                                break;
                            }
                        }
                        if (isBondChange) {
                            eBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                                    ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                            affectedBondReactants = eBond;
                        }
                    }

                    for (IBond pBond : connectedProductBondsList) {
                        boolean isBondChange = true;
                        String attachedPAtomID = pBond.getConnectedAtom(productAtom).getID();
                        for (IBond eBond : connectedEductBondsList) {
                            String attachedEAtomID = eBond.getConnectedAtom(eductAtom).getID();
                            if (attachedPAtomID.equalsIgnoreCase(attachedEAtomID)) {
                                isBondChange = false;
                                break;
                            }
                        }
                        if (isBondChange) {
                            pBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                                    ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                            affectedBondProducts = pBond;
                        }
                    }
                    if (affectedBondReactants != null && affectedBondProducts != null) {
                        affectedBondReactants.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                        affectedBondReactants.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);
                        affectedBondProducts.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                        affectedBondProducts.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);

                        getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                        getReactionCenterSet().add(affectedBondReactants.getAtom(1));
                        getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                        getReactionCenterSet().add(affectedBondProducts.getAtom(1));

                        BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
                        getBondChangeList().add(bondChange);
                    }
                } else if (rMol.getConnectedBondsCount(eductAtom) == 0
                        && pMol.getConnectedBondsCount(productAtom) > 0) {

                    IAtom psudoAtom = new Atom("PsH");
                    IBond eBond = new Bond(psudoAtom, eductAtom, IBond.Order.SINGLE);
                    IBond pBond = connectedProductBondsList.iterator().next();

                    eBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND);
                    pBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    affectedBondReactants = eBond;
                    affectedBondProducts = pBond;

                    affectedBondProducts.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                    affectedBondProducts.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);

                    getReactionCenterSet().add(affectedBondProducts.getAtom(0));
                    getReactionCenterSet().add(affectedBondProducts.getAtom(1));

                    BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
//                        BondChange bondChange = new BondChange(null, affectedBondProducts);
                    getBondChangeList().add(bondChange);

                } else if (rMol.getConnectedBondsCount(eductAtom) > 0
                        && pMol.getConnectedBondsCount(productAtom) == 0) {

                    IBond eBond = connectedEductBondsList.iterator().next();
                    IAtom psudoAtom = new Atom("PsH");
                    IBond pBond = new Bond(psudoAtom, productAtom, IBond.Order.SINGLE);

                    eBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
                    pBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.PSEUDO_BOND);

                    affectedBondReactants = eBond;
                    affectedBondProducts = pBond;

                    affectedBondReactants.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER, true);
                    affectedBondReactants.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER, true);

                    getReactionCenterSet().add(affectedBondReactants.getAtom(0));
                    getReactionCenterSet().add(affectedBondReactants.getAtom(1));

                    BondChange bondChange = new BondChange(affectedBondReactants, affectedBondProducts);
//                        BondChange bondChange = new BondChange(affectedBondReactants, null);
                    getBondChangeList().add(bondChange);

                }
            }
        }
    }

    private synchronized void markUnMappedAtoms() {
        for (IAtomContainer acE : reactantSet.atomContainers()) {
            for (IBond eBond : acE.bonds()) {
                boolean isNotMapped = false;
                if ((!mappingMap.containsKey(eBond.getAtom(0))
                        && mappingMap.containsKey(eBond.getAtom(1)))
                        || (mappingMap.containsKey(eBond.getAtom(0))
                        && !mappingMap.containsKey(eBond.getAtom(1)))) {
                    isNotMapped = true;
                }
                IBond pBond = null;
                if (isNotMapped) {
                    eBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    BondChange bondChange = new BondChange(eBond, pBond);
                    getBondChangeList().add(bondChange);
                }
            }
        }

        for (IAtomContainer acP : productSet.atomContainers()) {
            for (IBond pBond : acP.bonds()) {
                boolean isNotMapped = false;
                if ((!mappingMap.containsValue(pBond.getAtom(0))
                        && mappingMap.containsValue(pBond.getAtom(1)))
                        || (mappingMap.containsValue(pBond.getAtom(0))
                        && !mappingMap.containsValue(pBond.getAtom(1)))) {
                    isNotMapped = true;
                }
                IBond eBond = null;
                if (isNotMapped) {
                    pBond.setProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION,
                            ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
                    BondChange bondChange = new BondChange(eBond, pBond);
                    getBondChangeList().add(bondChange);
                }
            }
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
            if (affectedBondReactants.getFlag(CDKConstants.ISINRING)
                    == affectedBondProducts.getFlag(CDKConstants.ISINRING)) {

                if ((!affectedBondReactants.getFlag(CDKConstants.ISAROMATIC)
                        && affectedBondProducts.getFlag(CDKConstants.ISAROMATIC))
                        || (affectedBondReactants.getFlag(CDKConstants.ISAROMATIC)
                        && !affectedBondProducts.getFlag(CDKConstants.ISAROMATIC))) {
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
            if (affectedBondReactants.getFlag(CDKConstants.ISINRING)
                    == affectedBondProducts.getFlag(CDKConstants.ISINRING)) {

                if ((!affectedBondReactants.getFlag(CDKConstants.ISAROMATIC)
                        && affectedBondProducts.getFlag(CDKConstants.ISAROMATIC))
                        || (affectedBondReactants.getFlag(CDKConstants.ISAROMATIC)
                        && !affectedBondProducts.getFlag(CDKConstants.ISAROMATIC))) {
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
                                        if (Math.abs(countR - countP) <= 1 && sizeR == sizeP) {
                                            return 1;
                                        }
                                    }
                                    if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                        if (Math.abs(countR - countP) == 0 && sizeR == sizeP) {
                                            return 0;
                                        }
                                    }
                                }
                            } else if (atomR1.getHybridization().equals(atomP1.getHybridization())
                                    && atomR2.getHybridization().equals(atomP2.getHybridization())) {

                                if (!affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (Math.abs(countR - countP) <= 1
                                            && sizeR == sizeP) {
                                        return 0;
                                    }
                                }
                                if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (Math.abs(countR - countP) > 1 && sizeR == sizeP) {
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
