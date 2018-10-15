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

 /*
 * Maximize.java
 *
 * Created on February 3, 2006, 12:06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package uk.ac.ebi.reactionblast.mapping.algorithm.checks;

import java.io.Serializable;
import static java.lang.Double.MIN_VALUE;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;
import static java.util.Collections.synchronizedList;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class ChooseWinner extends Selector implements Serializable {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 0x296558709L;
    private EBIMatrix stereoMatrix;
    private EBIMatrix energyMatrix;
    private EBIMatrix carbonOverlapMatrix;

    private EBIMatrix similarityMatrix = null;
    private List<Cell> crossMappingTracer = null;
    private Map<Integer, IAtomContainer> educts = null;
    private Map<Integer, IAtomContainer> products = null;

    //~--- constructors -------------------------------------------------------
    /**
     * Creates a new instance of Maximize
     *
     * @param eductNameList
     * @param productNameList
     */
    public ChooseWinner(List<String> eductNameList,
            List<String> productNameList) {

        this.rowSize = eductNameList.size();
        this.colSize = productNameList.size();
        this.flagMatrix = new boolean[rowSize][colSize];
    }

    /**
     * @return the stereoMatrix
     */
    public EBIMatrix getStereoMatrix() {
        return stereoMatrix;
    }

    /**
     * @param stereoMatrix the stereoMatrix to set
     */
    public void setStereoMatrix(EBIMatrix stereoMatrix) {
        this.stereoMatrix = stereoMatrix;
    }

    /**
     * @return the energyMatrix
     */
    public EBIMatrix getEnergyMatrix() {
        return energyMatrix;
    }

    /**
     * @return the energyMatrix
     */
    public EBIMatrix getCarbonOverlapMatrix() {
        return carbonOverlapMatrix;
    }

    /**
     * @param energyMatrix the energyMatrix to set
     */
    public void setEnergyMatrix(EBIMatrix energyMatrix) {
        this.energyMatrix = energyMatrix;
    }

    /**
     *
     * @param eductMap
     * @param productMap
     * @param mHolder
     */
    public synchronized void searchWinners(Map<Integer, IAtomContainer> eductMap, Map<Integer, IAtomContainer> productMap, Holder mHolder) {
        initFlagMatrix();
        this.educts = eductMap;
        this.products = productMap;
        this.similarityMatrix = mHolder.getGraphSimilarityMatrix();
        this.setStereoMatrix(mHolder.getStereoMatrix());
        this.setEnergyMatrix(mHolder.getEnergyMatrix());
        this.setCarbonOverlapMatrix(mHolder.getCarbonOverlapMatrix());

        this.crossMappingTracer = synchronizedList(new ArrayList<Cell>());

        boolean isMappingFesiable = checkStatusFlag();
        List<Double> scores = new ArrayList<>();

//        System.out.println("isMappingFesiable " + isMappingFesiable);
        if (isMappingFesiable) {
            double maximumSimilarity = 0.0;

            boolean maxValueI;
            boolean maxValueJ;

            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {

                    double similarity = similarityMatrix.getValue(i, j);
                    //System.out.println("similarity "+similarity);
                    //matrix.
                    if (similarity > MIN_VALUE) {
                        maxValueI = isMajorSubgraphRow(similarityMatrix, i, j);
                        maxValueJ = isMajorSubgraphColumn(similarityMatrix, i, j);

//                        System.out.println("similarity " + similarity);
//                        System.out.println("i= " + i + ", j=" + j);
//                        System.out.println(" maxValue I " + maxValueI);
//                        System.out.println(" maxValue J " + maxValueJ);
                        if (maxValueI && maxValueJ) {

                            if (similarity > maximumSimilarity) {
                                maximumSimilarity = similarity;
                                initFlagMatrix();
                                scores.clear();
                            }
                            if (similarity == maximumSimilarity) {
                                this.flagMatrix[i][j] = true;
                                scores.add(similarity);
                            }
                        }
                    }
                }
            }
        }
        if (DEBUG) {
            System.out.println("resolveDeadLocks");
        }
        resolveDeadLocks(scores);
        if (DEBUG) {
            System.out.println("setWinOverFlags");
        }
        setWinOverFlags();
        if (DEBUG) {
            System.out.println("Done");
        }
    }

    /*
    * @return true if a cell in the this.flagMatrix
    * was set true else false
     */
    /**
     *
     * @return
     */
    public synchronized boolean getFlag() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                if (this.flagMatrix[i][j]) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     *
     * @return this.flagMatrix containing the selected cells
     *
     */
    public synchronized boolean[][] getFlagMatrix() {
        return this.flagMatrix;
    }

    private synchronized void initFlagMatrix() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                this.flagMatrix[i][j] = false;
            }
        }
    }

    private synchronized void initFlagMatrix(boolean[][] boolMatrix, int rowSize, int colSize) {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                boolMatrix[i][j] = false;
            }
        }
    }

    private synchronized boolean checkStatusFlag() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                if (similarityMatrix.getValue(i, j) > MIN_VALUE) {
                    return true;
                }
            }
        }
        return false;
    }

    private synchronized void setWinOverFlags() {
        educts.keySet().stream().forEach((indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {
                Cell cell = new Cell();
                cell.indexI = indexI;
                cell.indexJ = indexJ;
                cell.eductName = educts.get(indexI).getID();
                cell.productName = products.get(indexJ).getID();
                if (this.flagMatrix[indexI][indexJ]) {
                    if (cell.eductName.equalsIgnoreCase(cell.productName)) {
                    } else {
                        this.flagMatrix[indexI][indexJ] = !checkTwinMapping(cell);
                    }
                }
            });
        });
    }

    private synchronized boolean checkTwinMapping(Cell refCell) {
        boolean _statusFlag = false;
        for (Cell cell : crossMappingTracer) {

            //System.out.println("Visitor AtomContainer " + moleculeName + " i:" + _value.get(1) + " j:" + j);
            if (cell.eductName.equals(refCell.eductName)
                    && cell.productName.equals(refCell.productName)) {
                if (cell.indexI == refCell.indexI
                        || cell.indexJ == refCell.indexJ) {
                    _statusFlag = true;
                    break;
                }
            }
        }

        if (!_statusFlag) {
            crossMappingTracer.add(refCell);
        }
        return _statusFlag;
    }

    private synchronized void resolveDeadLocks(List<Double> scores) {
        boolean[][] deadlockFreeFlagMatrix = new boolean[rowSize][colSize];
        initFlagMatrix(deadlockFreeFlagMatrix, rowSize, colSize);
        scores.stream().map((score) -> new DeadLockResolver().resolver(score)).forEach((choosenCell) -> {
            deadlockFreeFlagMatrix[choosenCell.indexI][choosenCell.indexJ] = true;
        });
        this.flagMatrix = deadlockFreeFlagMatrix;
    }

    /**
     * @param carbonOverlapMatrix the carbonOverlapMatrix to set
     */
    public void setCarbonOverlapMatrix(EBIMatrix carbonOverlapMatrix) {
        this.carbonOverlapMatrix = carbonOverlapMatrix;
    }

    /**
     * Chosen cell of the matrix is stored here
     */
    class Cell {

        @Override
        public String toString() {
            return "Cells{" + "eductName=" + eductName
                    + ", productName=" + productName
                    + ", indexI=" + indexI
                    + ", indexJ=" + indexJ + '}';
        }

        String eductName;
        String productName;
        int indexI;
        int indexJ;

    }

    /**
     * Resolves deadlocks if more than one cell clashes with same scores. The
     * decision is then made on the max. stereo score and max. energy score
     */
    class DeadLockResolver {

        private synchronized double getMaxStereo(List<ChooseWinner.Cell> choosenCells) {
            double max = 0.0;
            for (ChooseWinner.Cell cell : choosenCells) {
                double val = getStereoMatrix().getValue(cell.indexI, cell.indexJ);
                if (val > max) {
                    max = val;
                }
            }
            return max;
        }

        private synchronized double getMinEnergy(List<ChooseWinner.Cell> choosenCells) {
            double min = Double.MAX_VALUE;
            for (ChooseWinner.Cell cell : choosenCells) {
                double val = getEnergyMatrix().getValue(cell.indexI, cell.indexJ);
                if (val < min) {
                    min = val;
                }
            }
            return min;
        }

        private synchronized double getMaxCarbonOverlap(List<ChooseWinner.Cell> choosenCells) {
            double max = Double.MIN_VALUE;
            for (ChooseWinner.Cell cell : choosenCells) {
                double val = getCarbonOverlapMatrix().getValue(cell.indexI, cell.indexJ);
                if (val > max) {
                    max = val;
                }
            }
            return max;
        }

        /**
         * Returns best match for the given score
         *
         * @param choosenScore
         * @return
         */
        public synchronized ChooseWinner.Cell resolver(double choosenScore) {
            List<ChooseWinner.Cell> choosenCells = new ArrayList<>();
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    if (flagMatrix[i][j] && choosenScore > MIN_VALUE
                            && similarityMatrix.getValue(i, j) == choosenScore) {
                        ChooseWinner.Cell cell = new ChooseWinner.Cell();
                        cell.indexI = i;
                        cell.indexJ = j;
                        cell.eductName = educts.get(i).getID();
                        cell.productName = products.get(j).getID();
                        choosenCells.add(cell);
                    }
                }
            }

            Cell winner = choosenCells.listIterator().next();

            double maxStereo = getMaxStereo(choosenCells);
            double minEnergy = getMinEnergy(choosenCells);
            double maxCarbonOverlap = getMaxCarbonOverlap(choosenCells);
            if (DEBUG) {
                System.out.println("maxStereo " + maxStereo);
            }
            if (DEBUG) {
                System.out.println("minEnergy " + minEnergy);
            }
            if (DEBUG) {
                System.out.println("maxCarbon " + maxCarbonOverlap);
            }
            for (ChooseWinner.Cell cell : choosenCells) {
                double stereoVal = stereoMatrix.getValue(cell.indexI, cell.indexJ);
                double energyVal = energyMatrix.getValue(cell.indexI, cell.indexJ);
                double carbonOverlap = carbonOverlapMatrix.getValue(cell.indexI, cell.indexJ);

                if (DEBUG) {
                    System.out.println("Comparing " + cell.indexI + "," + cell.indexJ);
                }

                if (DEBUG) {
                    System.out.print("maxStereo " + maxStereo);
                    System.out.println(", stereoVal " + stereoVal);
                }
                if (DEBUG) {
                    System.out.print("minEnergy " + minEnergy);
                    System.out.println(", energyVal " + energyVal);
                }

                if (DEBUG) {
                    System.out.print("maxCarbon " + maxCarbonOverlap);
                    System.out.println(", Carbon " + carbonOverlap);
                }

                if (stereoVal >= maxStereo
                        && energyVal <= minEnergy
                        && (maxCarbonOverlap == carbonOverlap
                        || maxCarbonOverlap == 0)) {
                    winner = cell;
                    break;
                }

            }
            return winner;
        }
    }
}
