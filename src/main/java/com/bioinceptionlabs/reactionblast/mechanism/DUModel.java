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
import static java.util.Collections.synchronizedList;
import static java.util.Collections.synchronizedMap;
import static java.util.Collections.synchronizedSet;
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
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.helper.MoleculeInitializer;
import com.bioinceptionlabs.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import com.bioinceptionlabs.reactionblast.mechanism.helper.AtomStereoChangeInformation;
import com.bioinceptionlabs.reactionblast.mechanism.helper.BondChange;
import com.bioinceptionlabs.reactionblast.mechanism.interfaces.IChangeCalculator;
import com.bioinceptionlabs.reactionblast.stereo.IStereoAndConformation;
import com.bioinceptionlabs.reactionblast.stereo.ebi.StereoCenteralityTool;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
abstract class DUModel extends StereoCenteralityTool implements IChangeCalculator, Serializable {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(DUModel.class);

    private static final long serialVersionUID = 179876660968690L;
    final IAtomContainerSet reactantSet;
    final IAtomContainerSet productSet;
    final Map<IAtom, IAtom> mappingMap;
    final List<BondChange> bondChangeList;
    final Set<IAtom> reactionCenterList;
    final List<AtomStereoChangeInformation> stereoChangeList;
    final List<AtomStereoChangeInformation> conformationChangeList;
    final List<StereoChange> stereogenicCenters;
    protected final boolean generate3DCoordinates;
    protected final boolean generate2DCoordinates;
    protected final AtomAtomMappingContainer mapping;
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
        this.bondChangeList = synchronizedList(new ArrayList<>());
        this.reactionCenterList = synchronizedSet(new LinkedHashSet<>());
        this.stereoChangeList = synchronizedList(new ArrayList<>());
        this.conformationChangeList = synchronizedList(new ArrayList<>());
        this.mappingMap = synchronizedMap(new HashMap<>());

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
            this.mapping = new AtomAtomMappingContainer(reaction, withoutHydrogen);
            LOGGER.debug("=====createRMatrix=====");
            this.reactionMatrix = createRMatrix(reactantBE, productBE, mapping);
        } catch (Exception e) {
            throw new Exception("WARNING: Unable to compute reaction matrix", e);
        }
        /*
         * Stereo mapping
         */
        LOGGER.debug("Assign Stereo");
        Map<IAtom, IStereoAndConformation> chiralityCDK2D = new HashMap<>();
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
    private RMatrix createRMatrix(BEMatrix reactantBE, BEMatrix productBE, AtomAtomMappingContainer mapping) throws Exception {
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
}
