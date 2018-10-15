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

import java.io.IOException;
import java.io.Serializable;
import static java.lang.System.arraycopy;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import static java.util.logging.Level.SEVERE;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.mapping.container.HydrogenFreeFingerPrintContainer;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.container.helper.MolMapping;
import uk.ac.ebi.reactionblast.mapping.helper.Debugger;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;
import static java.util.Collections.synchronizedList;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Holder extends Debugger implements Cloneable, Serializable {

    private static final boolean DEBUG = false;
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
        if (DEBUG) {
            out.println("setFingerprint");
        }
        setFingerprint();
        if (DEBUG) {
            out.println("setMolMapping");
        }
        setMolMapping();
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
        this.mappingMolPair = synchronizedList(new ArrayList<>());
        if (DEBUG) {
            out.println("initialize the Matrix");
        }
        initialize();
    }

    /**
     * @return the stereoMatrix
     */
    public synchronized EBIMatrix getStereoMatrix() {
        return stereoMatrix;
    }

    /**
     * @return the cliqueMatrix
     */
    public synchronized EBIMatrix getCliqueMatrix() {
        return cliqueMatrix;
    }

    /**
     * @return the graphSimilarityMatrix
     */
    public synchronized EBIMatrix getGraphSimilarityMatrix() {
        return graphSimilarityMatrix;
    }

    /**
     * @return the fragmentMatrix
     */
    public synchronized EBIMatrix getFragmentMatrix() {
        return fragmentMatrix;
    }

    /**
     * @return the energyMatrix
     */
    public synchronized EBIMatrix getEnergyMatrix() {
        return energyMatrix;
    }

    private void initialize() {
//        System.out.println("\nInitialize Matrix with Zero\n");
        graphSimilarityMatrix.initMatrix(0.0);
        stereoMatrix.initMatrix(0.0);
        cliqueMatrix.initMatrix(0.0);
        fragmentMatrix.initMatrix(0.0);
        fpSimMatrixWithoutHydrogen.initMatrix(0.0);
        energyMatrix.initMatrix(0.0);
        carbonOverlapMatrix.initMatrix(0.0);
    }

    private void setFingerprint() {
        for (int i = 0; i < eductCounter.size(); i++) {
            for (int j = 0; j < productCounter.size(); j++) {
                try {
                    String eductName = eductCounter.get(i).trim();
                    String productName = productCounter.get(j).trim();
                    BitSet hydrogenEductFP = hydFPFree.getFingerPrint(eductName);
                    BitSet hydrogenProductFP = hydFPFree.getFingerPrint(productName);
                    float hydrogenSimVal = getTanimotoSimilarity(hydrogenEductFP, hydrogenProductFP);
                    if (DEBUG) {
                        out.println("FP " + hydrogenSimVal);
                    }
                    fpSimMatrixWithoutHydrogen.setValue(i, j, hydrogenSimVal);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
    }

    private void setMolMapping() {
        for (int i = 0; i < eductCounter.size(); i++) {
            for (int j = 0; j < productCounter.size(); j++) {
                try {
                    String eductName = eductCounter.get(i).trim();
                    String productName = productCounter.get(j).trim();
                    MolMapping m = new MolMapping(eductName, productName, i, j);
                    getMappingMolPair().add(m);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
    }

    /**
     * @return the fpSimMatrixWithoutHydrogen
     */
    public synchronized EBIMatrix getFPSimilarityMatrix() {
        return fpSimMatrixWithoutHydrogen;
    }

    /**
     * @return the structureInformation
     */
    public synchronized ReactionContainer getReactionContainer() {
        return structureInformation;
    }

    /**
     * @return the mappingMolPair
     */
    public synchronized List<MolMapping> getMappingMolPair() {
        return mappingMolPair;
    }

    /**
     * Cloned EBIMatrix Objects
     *
     * @return
     * @throws CloneNotSupportedException
     */
    @Override
    public synchronized Object clone() throws CloneNotSupportedException {
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
    public synchronized BestMatch getBestMatchContainer() {
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
