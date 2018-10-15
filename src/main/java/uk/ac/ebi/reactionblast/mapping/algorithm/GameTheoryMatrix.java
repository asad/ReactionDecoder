/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.algorithm;


/*
 * GameTheoryMatrix.java
 *
 * Created on 13 February 2006, 11:58
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
import java.io.IOException;
import java.util.BitSet;
import static java.util.Collections.sort;
import static java.util.Collections.synchronizedList;
import static java.util.Collections.synchronizedSortedMap;
import static java.util.Collections.unmodifiableList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import static uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator.getFingerprinterSize;
import uk.ac.ebi.reactionblast.mapping.container.BestMatchContainer;
import uk.ac.ebi.reactionblast.mapping.container.HydrogenFreeFingerPrintContainer;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.mapping.interfaces.IGraphTheoryMatrix;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.tools.AtomContainerSetComparator;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

//~--- classes ----------------------------------------------------------------
/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class GameTheoryMatrix extends BaseGameTheory implements IGraphTheoryMatrix {

    private static final long serialVersionUID = 0x2c36427fd2L;
    //~--- constructors -------------------------------------------------------
    private static final ILoggingTool LOGGER = createLoggingTool(GameTheoryMatrix.class);
//    private void initializeMappingFLAGS(Holder mh) throws Exception {
//        ReactionContainer reactionStructureInformation = mh.getReactionContainer();
    /*Reset all the flags*/
//        for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
//            for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
//                reactionStructureInformation.setEductModified(substrateIndex, true);
//                reactionStructureInformation.setProductModified(productIndex, true);
//            }
//        }
//    }
    private Holder matrixHolder;
    private MoleculeMoleculeMapping reactionBlastMolMapping;
    private final List<String> eductCounter;
    private final List<String> productCounter;
    private final IReaction reaction;
    private final Map<Integer, BitSet> substrateductFPMap;
    private final Map<Integer, BitSet> productFPMap;
    private final FingerprintGenerator fpr;
    private final HydrogenFreeFingerPrintContainer hydFreeFPContainer;
    private final boolean removeHydrogen;
    private final String reactionID;
    private final ReactionContainer structureMapObj;
    private final BestMatch bestMatchContainer;
    private final IMappingAlgorithm theory;

    /**
     * Creates a new instance of GameTheoryMatrix
     *
     * @param theory
     * @param reaction
     * @param removeHydrogen
     * @throws Exception
     */
    public GameTheoryMatrix(
            IMappingAlgorithm theory,
            IReaction reaction,
            boolean removeHydrogen) throws Exception {
        this.theory = theory;
        this.removeHydrogen = removeHydrogen;
        this.reaction = reaction;
        this.reactionID = reaction.getID();

        this.substrateductFPMap = synchronizedSortedMap(new TreeMap<>());
        this.productFPMap = synchronizedSortedMap(new TreeMap<>());
        this.fpr = new FingerprintGenerator();
        this.eductCounter = synchronizedList(new LinkedList<>());
        this.productCounter = synchronizedList(new LinkedList<>());
        this.structureMapObj = new ReactionContainer();
        this.bestMatchContainer = new BestMatchContainer();
        this.hydFreeFPContainer = new HydrogenFreeFingerPrintContainer();
        this.reactionBlastMolMapping = new MoleculeMoleculeMapping();

        try {
            StoichiometricCoefficientReplicator_Structure_FingerPrint_MapGenerator();
            BuildScoringMatrix();
        } catch (Exception e) {
            LOGGER.error(e);
        }
    }
    //~--- methods ------------------------------------------------------------

    private synchronized void BuildScoringMatrix() throws Exception {
        try {
            matrixHolder = new Holder(
                    theory,
                    reactionID,
                    eductCounter,
                    productCounter,
                    structureMapObj,
                    bestMatchContainer,
                    hydFreeFPContainer);
            this.reactionBlastMolMapping.setMolMappings(reactionID, matrixHolder.getMappingMolPair());
            /*
             * Set FLAGS to True, to allow MCS calculation
             * This step requires ALL-BY-ALL mappings, hence
             * all the falgs are true
             * 
             */
//            initializeMappingFLAGS(matrixHolder);
            UpdateMatrix(matrixHolder, removeHydrogen);
        } catch (Exception e) {
            LOGGER.error(SEVERE, null, e);
        }
    }

    /**
     * clears all the containers in this object
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        structureMapObj.Clear();
        hydFreeFPContainer.Clear();
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @return reactant molecule and their index (key) in a Map
     */
    @Override
    public synchronized List<String> getEductCounter() {
        return unmodifiableList(eductCounter);
    }

    /**
     *
     * @return product molecule and their index (key) in a Map
     */
    @Override
    public synchronized List<String> getProductCounter() {
        return unmodifiableList(productCounter);
    }

    private synchronized void StoichiometricCoefficientReplicator_Structure_FingerPrint_MapGenerator() {

        List<IAtomContainer> ac = new LinkedList<>();
        List<IAtomContainer> pd = new LinkedList<>();
        sortAtomContainer(ac, pd);
        for (int key = 0; key < ac.size(); key++) {
            try {
                IAtomContainer mol = ac.get(key).clone();
                String eductID = ac.get(key).getID().trim();
                mol.setID(eductID);

                BitSet FP;
                if (hydFreeFPContainer.isKeyPresent(eductID)) {
                    FP = hydFreeFPContainer.getFingerPrint(eductID);
                } else if (mol.getAtomCount() > 0) {
                    IAtomContainer tempMol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
                    FP = fpr.getFingerprint(tempMol);
                } else {
                    FP = new BitSet(getFingerprinterSize());
                }
                hydFreeFPContainer.setValue(eductID, FP);
                eductCounter.add(key, eductID);
                structureMapObj.putEduct(key, mol);
                structureMapObj.setEductModified(key, true);
                substrateductFPMap.put(key, FP);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }

        }
        for (int key = 0; key < pd.size(); key++) {
            try {
                IAtomContainer mol = pd.get(key).clone();
                String productID = pd.get(key).getID().trim();
                mol.setID(productID);

                BitSet fingerPrint;

                if (hydFreeFPContainer.isKeyPresent(productID)) {
                    fingerPrint = hydFreeFPContainer.getFingerPrint(productID);
                } else if (mol.getAtomCount() > 0) {
                    IAtomContainer tempMol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
                    fingerPrint = fpr.getFingerprint(tempMol);
                } else {
                    fingerPrint = new BitSet(getFingerprinterSize());
                }
                hydFreeFPContainer.setValue(productID, fingerPrint);
                productCounter.add(key, productID);
                structureMapObj.putProduct(key, mol);
                structureMapObj.setProductModified(key, true);
                productFPMap.put(key, fingerPrint);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
    }

    /**
     * @return the matrixHolder
     */
    @Override
    public synchronized Holder getMatrixHolder() {
        return matrixHolder;
    }

    @Override
    public synchronized MoleculeMoleculeMapping getReactionMolMapping() {
        return reactionBlastMolMapping;
    }

    private synchronized void sortAtomContainer(List<IAtomContainer> ac, List<IAtomContainer> pd) {
        for (IAtomContainer e : reaction.getReactants().atomContainers()) {
//            System.out.println("ID: e " + e.getID());
            ac.add(e);
        }
        for (IAtomContainer p : reaction.getProducts().atomContainers()) {
//            System.out.println("ID: P " + p.getID());
            pd.add(p);
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            sort(ac, comparator);
            sort(pd, comparator);
        } catch (Exception e) {
            LOGGER.debug("ERROR: in AtomMappingTool: " + e.getMessage());
        }
    }

    @Override
    public synchronized int getDelta() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public synchronized void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) {
        this.reactionBlastMolMapping = reactionMolMapping;
    }

}
