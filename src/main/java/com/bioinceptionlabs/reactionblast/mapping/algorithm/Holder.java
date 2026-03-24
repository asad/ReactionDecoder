/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping.algorithm;

import java.io.IOException;
import java.io.Serializable;
import static java.lang.System.arraycopy;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import static java.util.logging.Level.SEVERE;
import static com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.Similarity.getTanimotoSimilarity;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.HydrogenFreeFingerPrintContainer;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.MolMapping;
import com.bioinceptionlabs.reactionblast.mapping.Reactor.Debugger;
import com.bioinceptionlabs.reactionblast.mapping.BestMatch;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import com.bioinceptionlabs.reactionblast.tools.EBIMatrix;

import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class Holder extends Debugger implements Cloneable, Serializable {

    private static final long serialVersionUID = 18989786786L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(Holder.class);
    /*
     Final methods
     */
    private final List<MolMapping> mappingMolPair;
    private final EBIMatrix stereoMatrix;
    private final EBIMatrix cliqueMatrix;
    private final EBIMatrix graphSimilarityMatrix;
    private final EBIMatrix fragmentMatrix;
    private final EBIMatrix energyMatrix;
    private final EBIMatrix carbonOverlapMatrix;
    private final EBIMatrix fpSimMatrixWithoutHydrogen;
    private final int row;
    private final int coloumn;

    /*
     Local methods
     */
    private ReactionContainer structureInformation;
    private BestMatch bestMatchContainer;
    private List<String> eductCounter;
    private List<String> productCounter;
    private String reactionID;
    private HydrogenFreeFingerPrintContainer hydFPFree;
    private IMappingAlgorithm theory;

    /**
     *
     * @param theory
     * @param reactionID
     * @param eductCounter
     * @param productCounter
     * @param reactionContainer
     * @param bestMatchContainer
     * @param hydFPFree
     * @throws IOException
     */
    public Holder(
            IMappingAlgorithm theory,
            String reactionID,
            List<String> eductCounter,
            List<String> productCounter,
            ReactionContainer reactionContainer,
            BestMatch bestMatchContainer,
            HydrogenFreeFingerPrintContainer hydFPFree) throws IOException {

        this(eductCounter.size(), productCounter.size());
        this.theory = theory;
        this.reactionID = reactionID;
        this.eductCounter = eductCounter;
        this.productCounter = productCounter;
        this.structureInformation = reactionContainer;
        this.bestMatchContainer = bestMatchContainer;
        this.hydFPFree = hydFPFree;
        LOGGER.debug("setFingerprintAndMolMapping");
        setFingerprintAndMolMapping();
    }

    /**
     *
     * @param row
     * @param column
     */
    public Holder(int row, int column) {
        this.row = row;
        this.coloumn = column;
        this.graphSimilarityMatrix = new EBIMatrix(row, column);
        this.stereoMatrix = new EBIMatrix(row, column);
        this.cliqueMatrix = new EBIMatrix(row, column);
        this.fragmentMatrix = new EBIMatrix(row, column);
        this.carbonOverlapMatrix = new EBIMatrix(row, column);
        this.fpSimMatrixWithoutHydrogen = new EBIMatrix(row, column);
        this.energyMatrix = new EBIMatrix(row, column);
        this.mappingMolPair = new ArrayList<>();
        LOGGER.debug("initialize the Matrix");
        initialize();
    }

    /**
     * @return the stereoMatrix
     */
    public EBIMatrix getStereoMatrix() {
        return stereoMatrix;
    }

    /**
     * @return the cliqueMatrix
     */
    public EBIMatrix getCliqueMatrix() {
        return cliqueMatrix;
    }

    /**
     * @return the graphSimilarityMatrix
     */
    public EBIMatrix getGraphSimilarityMatrix() {
        return graphSimilarityMatrix;
    }

    /**
     * @return the fragmentMatrix
     */
    public EBIMatrix getFragmentMatrix() {
        return fragmentMatrix;
    }

    /**
     * @return the energyMatrix
     */
    public EBIMatrix getEnergyMatrix() {
        return energyMatrix;
    }

    private void initialize() {
        graphSimilarityMatrix.initMatrix(0.0);
        stereoMatrix.initMatrix(0.0);
        cliqueMatrix.initMatrix(0.0);
        fragmentMatrix.initMatrix(0.0);
        fpSimMatrixWithoutHydrogen.initMatrix(0.0);
        energyMatrix.initMatrix(0.0);
        carbonOverlapMatrix.initMatrix(0.0);
    }

    private void setFingerprintAndMolMapping() {
        // Pre-cache trimmed names and fingerprints to avoid redundant lookups
        int eSize = eductCounter.size();
        int pSize = productCounter.size();
        String[] eNames = new String[eSize];
        BitSet[] eFPs = new BitSet[eSize];
        for (int i = 0; i < eSize; i++) {
            eNames[i] = eductCounter.get(i).trim();
            try {
                eFPs[i] = hydFPFree.getFingerPrint(eNames[i]);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        String[] pNames = new String[pSize];
        BitSet[] pFPs = new BitSet[pSize];
        for (int j = 0; j < pSize; j++) {
            pNames[j] = productCounter.get(j).trim();
            try {
                pFPs[j] = hydFPFree.getFingerPrint(pNames[j]);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        // Single combined loop for fingerprint similarity + mol mapping
        for (int i = 0; i < eSize; i++) {
            for (int j = 0; j < pSize; j++) {
                try {
                    if (eFPs[i] != null && pFPs[j] != null) {
                        float hydrogenSimVal = getTanimotoSimilarity(eFPs[i], pFPs[j]);
                        LOGGER.debug("FP " + hydrogenSimVal);
                        fpSimMatrixWithoutHydrogen.setValue(i, j, hydrogenSimVal);
                    }
                    getMappingMolPair().add(new MolMapping(eNames[i], pNames[j], i, j));
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
    }

    /**
     * @return the fpSimMatrixWithoutHydrogen
     */
    public EBIMatrix getFPSimilarityMatrix() {
        return fpSimMatrixWithoutHydrogen;
    }

    /**
     * @return the structureInformation
     */
    public ReactionContainer getReactionContainer() {
        return structureInformation;
    }

    /**
     * @return the mappingMolPair
     */
    public List<MolMapping> getMappingMolPair() {
        return mappingMolPair;
    }

    /**
     * Cloned EBIMatrix Objects
     *
     * @return
     * @throws CloneNotSupportedException
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        Holder mhClone = new Holder(this.row, this.coloumn);
        mhClone.setTheory(this.getTheory());

        double[][] arrayCopy = this.getGraphSimilarityMatrix().getArrayCopy();
        EBIMatrix matrix = mhClone.getGraphSimilarityMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getFragmentMatrix().getArrayCopy();
        matrix = mhClone.getFragmentMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getStereoMatrix().getArrayCopy();
        matrix = mhClone.getStereoMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getCliqueMatrix().getArrayCopy();
        matrix = mhClone.getCliqueMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getEnergyMatrix().getArrayCopy();
        matrix = mhClone.getEnergyMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getCarbonOverlapMatrix().getArrayCopy();
        matrix = mhClone.getCarbonOverlapMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        arrayCopy = this.getFPSimilarityMatrix().getArrayCopy();
        matrix = mhClone.getFPSimilarityMatrix();
        setData(arrayCopy, matrix, row, coloumn);

        mhClone.structureInformation = this.getReactionContainer();
        mhClone.bestMatchContainer = this.getBestMatchContainer();
        return mhClone;
    }

    /**
     * Make a deep duplicate of a matrix
     *
     * @param sourceData
     * @param sinkMatrix
     * @param rows
     * @param coloumns
     */
    private void setData(double[][] sourceData, EBIMatrix sinkMatrix, int rows, int coloumns) {
        double[][] newDataMatrix = sinkMatrix.getArray();
        for (int i = 0; i < rows; i++) {
            arraycopy(sourceData[i], 0, newDataMatrix[i], 0, coloumns);
        }
    }

    /**
     * @return the bestMatchContainer
     */
    public BestMatch getBestMatchContainer() {
        return bestMatchContainer;
    }

    /**
     * @return the theory
     */
    public IMappingAlgorithm getTheory() {
        return theory;
    }

    /**
     * @param theory the theory to set
     */
    public void setTheory(IMappingAlgorithm theory) {
        this.theory = theory;
    }

    /**
     * @return the carbonOverlapMatrix
     */
    public EBIMatrix getCarbonOverlapMatrix() {
        return carbonOverlapMatrix;
    }
}
