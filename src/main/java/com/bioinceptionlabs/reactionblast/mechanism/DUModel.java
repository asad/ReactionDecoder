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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.MoleculeInitializer;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
abstract class DUModel extends MechanismHelpers.Utility implements IChangeCalculator, Serializable {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(DUModel.class);
    private static final ILoggingTool STEREO_LOGGER
            = LoggingToolFactory.createLoggingTool(DUModel.class);

    // ---- Methods merged from StereoCenteralityTool.java ----

    private static IAtom getAtomByID(String id, IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            if (a.getID() != null && a.getID().equals(id)) {
                return a;
            }
        }
        return null;
    }

    @SuppressWarnings("unchecked")
    static Map<IAtom, BondChangeCalculator.IStereoAndConformation> getChirality2D(IReaction reaction)
            throws CDKException, CloneNotSupportedException {
        Map<IAtom, BondChangeCalculator.IStereoAndConformation> chiralityMap = new HashMap<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
            Map<IAtom, BondChangeCalculator.IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH);
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
            Map<IAtom, BondChangeCalculator.IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH);
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        return chiralityMap;
    }

    @SuppressWarnings("unchecked")
    static Map<IAtom, BondChangeCalculator.IStereoAndConformation> getChirality2D(IAtomContainer ac) {
        Map<IAtom, BondChangeCalculator.IStereoAndConformation> chiralityMap = new HashMap<>();
        for (IAtom atom : ac.atoms()) {
            chiralityMap.put(atom, BondChangeCalculator.IStereoAndConformation.NONE);
        }
        try {
            for (IStereoElement<?, ?> element : ac.stereoElements()) {
                if (element instanceof ITetrahedralChirality) {
                    ITetrahedralChirality tc = (ITetrahedralChirality) element;
                    IAtom focus = tc.getChiralAtom();
                    ITetrahedralChirality.Stereo stereo = tc.getStereo();
                    if (stereo == ITetrahedralChirality.Stereo.CLOCKWISE) {
                        chiralityMap.put(focus, BondChangeCalculator.IStereoAndConformation.R);
                    } else if (stereo == ITetrahedralChirality.Stereo.ANTI_CLOCKWISE) {
                        chiralityMap.put(focus, BondChangeCalculator.IStereoAndConformation.S);
                    }
                    if (focus != null) {
                        focus.setProperty("Stereo", chiralityMap.get(focus));
                    }
                } else if (element instanceof IDoubleBondStereochemistry) {
                    IDoubleBondStereochemistry dbs = (IDoubleBondStereochemistry) element;
                    IDoubleBondStereochemistry.Conformation conf = dbs.getStereo();
                    IAtom a0 = dbs.getStereoBond().getBegin();
                    IAtom a1 = dbs.getStereoBond().getEnd();
                    BondChangeCalculator.IStereoAndConformation sc;
                    if (conf == IDoubleBondStereochemistry.Conformation.OPPOSITE) {
                        sc = BondChangeCalculator.IStereoAndConformation.E;
                    } else if (conf == IDoubleBondStereochemistry.Conformation.TOGETHER) {
                        sc = BondChangeCalculator.IStereoAndConformation.Z;
                    } else {
                        continue;
                    }
                    chiralityMap.put(a0, sc);
                    chiralityMap.put(a1, sc);
                    if (a0 != null) a0.setProperty("Stereo", sc);
                    if (a1 != null) a1.setProperty("Stereo", sc);
                }
            }
        } catch (Exception e) {
            STEREO_LOGGER.debug("Stereo perception error: " + e.getMessage());
        }
        chiralityMap.forEach((atom, sc) -> {
            if (atom != null) atom.setProperty("Stereo", sc);
        });
        return chiralityMap;
    }

    // ---- End of StereoCenteralityTool methods ----

    private static final long serialVersionUID = 179876660968690L;
    final IAtomContainerSet reactantSet;
    final IAtomContainerSet productSet;
    final Map<IAtom, IAtom> mappingMap;
    final List<MechanismHelpers.BondChange> bondChangeList;
    final Set<IAtom> reactionCenterList;
    final List<MechanismHelpers.AtomStereoChangeInformation> stereoChangeList;
    final List<MechanismHelpers.AtomStereoChangeInformation> conformationChangeList;
    final List<DUModel.StereoChange> stereogenicCenters;
    protected final boolean generate3DCoordinates;
    protected final boolean generate2DCoordinates;
    protected final MechanismHelpers.AtomAtomMappingContainer mapping;
    protected final BEMatrix reactantBE;
    protected final BEMatrix productBE;
    protected final RMatrix reactionMatrix;
    protected final IRingSet queryRingSet;
    protected final IRingSet targetRingSet;

    /**
     *
     * @param reaction
     * @param withoutHydrogen
     * @param generate2D
     * @param generate3D
     * @throws CDKException
     * @throws Exception
     */
    DUModel(IReaction reaction,
            boolean withoutHydrogen,
            boolean generate2D,
            boolean generate3D) throws CDKException, Exception {

        this.reactantSet = reaction.getReactants();
        this.productSet = reaction.getProducts();
        this.bondChangeList = new ArrayList<>();
        this.reactionCenterList = new LinkedHashSet<>();
        this.stereoChangeList = new ArrayList<>();
        this.conformationChangeList = new ArrayList<>();
        this.mappingMap = new HashMap<>();

        this.generate3DCoordinates = generate3D;
        this.generate2DCoordinates = generate2D;
        this.queryRingSet = new RingSet();
        this.targetRingSet = new RingSet();

        /*
         Set Atom-Atom Mapping
         */
        LOGGER.debug("setMappingMap");
        setMappingMap(reaction.mappings());
        LOGGER.debug("Done setMappingMap");

        LOGGER.debug("Mark Aromatic Bonds");
        List<IBond> rBonds = new ArrayList<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            /*
             * Aromatise and mark rings (imp for detecting keukal changes)
             */
            LOGGER.debug("MoleculeInitializer");
            MoleculeInitializer.initializeMolecule(ac);
            LOGGER.debug("MoleculeInitializer Done");
            for (IBond bond : ac.bonds()) {
                rBonds.add(bond);
            }
        }
        List<IBond> pBonds = new ArrayList<>();
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            /*
             Aromatise and mark rings (imp for detecting keukal changes)
             */
            LOGGER.debug("MoleculeInitializer");
            MoleculeInitializer.initializeMolecule(ac);
            LOGGER.debug("Done");
            for (IBond bond : ac.bonds()) {
                pBonds.add(bond);
            }
        }

        LOGGER.debug("Done Marking Aromatic Bonds");

        try {

            LOGGER.debug("=====Educt createBEMatrix=====");
            this.reactantBE = createBEMatrix(reactantSet, rBonds, withoutHydrogen, mappingMap);
            LOGGER.debug("=====Product createBEMatrix=====");
            this.productBE = createBEMatrix(productSet, pBonds, withoutHydrogen, mappingMap);
            LOGGER.debug("=====AAM Container=====");
            this.mapping = new MechanismHelpers.AtomAtomMappingContainer(reaction, withoutHydrogen);
            LOGGER.debug("=====createRMatrix=====");
            this.reactionMatrix = createRMatrix(reactantBE, productBE, mapping);
        } catch (Exception e) {
            throw new Exception("WARNING: Unable to compute reaction matrix", e);
        }
        /*
         * Stereo mapping
         */
        LOGGER.debug("Assign Stereo");
        Map<IAtom, BondChangeCalculator.IStereoAndConformation> chiralityCDK2D = new HashMap<>();
        try {
            chiralityCDK2D = getChirality2D(reaction);
        } catch (CDKException | CloneNotSupportedException ex) {
            throw new Exception("WARNING: 2D CDK based stereo perception failed", ex);
        }
        LOGGER.debug("Done Assign Stereo");
        /*
         * Generate stereo information
         */
        LOGGER.debug("Assign Stereo Center");
        try {
            this.stereogenicCenters = new StereogenicCenterCalculator().compare(reaction, chiralityCDK2D);
        } catch (Exception e) {
            throw new Exception("WARNING: 2D CDK based stereo centers perception failed", e);
        }
        LOGGER.debug("Done Assign Stereo Center");
    }

    /**
     * @param mappings to be set
     */
    private void setMappingMap(Iterable<IMapping> mappings) {
        Iterator<IMapping> mappingIterator = mappings.iterator();
        while (mappingIterator.hasNext()) {
            IMapping mappingObject = mappingIterator.next();
            IAtom atomEduct = (IAtom) mappingObject.getChemObject(0);
            IAtom atomProduct = (IAtom) mappingObject.getChemObject(1);
            mappingMap.put(atomEduct, atomProduct);
        }
    }

    private BEMatrix createBEMatrix(IAtomContainerSet molset, List<IBond> bonds, boolean withoutH, Map<IAtom, IAtom> mappings) throws Exception {
        BEMatrix res = new BEMatrix(withoutH, molset, bonds, mappings);
        res.setMatrixAtoms();
        return res;
    }

    /**
     *
     * @param reactantBE
     * @param productBE
     * @param mapping
     * @return
     * @throws CDKException
     */
    private RMatrix createRMatrix(BEMatrix reactantBE, BEMatrix productBE, MechanismHelpers.AtomAtomMappingContainer mapping) throws Exception {
        return new RMatrix(reactantBE, productBE, mapping);
    }

    @Override
    public String toString() {
        return "DUModel{" + "reactantSet=" + reactantSet
                + ", productSet=" + productSet
                + ", mappingMap=" + mappingMap
                + ", bondChangeList=" + bondChangeList
                + ", reactionCenterList=" + reactionCenterList
                + ", stereoChangeList=" + stereoChangeList
                + ", conformationChangeList=" + conformationChangeList + '}';
    }

    /**
     * Stereo change data holder (merged from StereoChange.java).
     */
    static class StereoChange implements Serializable {

        private static final long serialVersionUID = 6778787889667901L;
        private final BondChangeCalculator.IStereoAndConformation rAtomStereo;
        private final BondChangeCalculator.IStereoAndConformation pAtomStereo;
        private final IAtom rAtom;
        private final IAtom pAtom;

        StereoChange(BondChangeCalculator.IStereoAndConformation rAtomStereo,
                BondChangeCalculator.IStereoAndConformation pAtomStereo,
                IAtom rAtom,
                IAtom pAtom) {
            this.rAtomStereo = rAtomStereo;
            this.pAtomStereo = pAtomStereo;
            this.rAtom = rAtom;
            this.pAtom = pAtom;
        }

        @Override
        public String toString() {
            return "StereoChange{" + "rAtomStereo=" + rAtomStereo + ", pAtomStereo=" + pAtomStereo + ", rAtom="
                    + rAtom.getSymbol() + rAtom.getID() + ", pAtom=" + pAtom.getSymbol() + pAtom.getID() + '}';
        }

        public BondChangeCalculator.IStereoAndConformation getReactantAtomStereo() {
            return rAtomStereo;
        }

        public BondChangeCalculator.IStereoAndConformation getProductAtomStereo() {
            return pAtomStereo;
        }

        public IAtom getReactantAtom() {
            return rAtom;
        }

        public IAtom getProductAtom() {
            return pAtom;
        }
    }
}
