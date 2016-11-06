/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.util.Collections.synchronizedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class ChooseWinner extends Selector implements Serializable {

    private static final long serialVersionUID = 0x296558709L;
    private static final Logger LOG = getLogger(ChooseWinner.class.getName());
    private EBIMatrix stereoMatrix;
    private EBIMatrix energyMatrix;

    private EBIMatrix similarityMatrix = null;
    private List<Cells> crossMappingTracer = null;
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
    public EBIMatrix getStereoMatrix(
            ) {
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
        this.crossMappingTracer = synchronizedList(new ArrayList<Cells>());

        boolean isMappingFesiable = checkStatusFlag();
        List<Double> scores = new ArrayList<>();

        //System.out.println("isMappingFesiable " + isMappingFesiable);
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
//        System.out.println("resolveDeadLocks");
        resolveDeadLocks(scores);
//        System.out.println("setWinOverFlags");
        setWinOverFlags();
//        System.out.println("Done");
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
        for (Integer indexI : educts.keySet()) {
            for (Integer indexJ : products.keySet()) {
                Cells cell = new Cells();
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
            }
        }
    }

    private synchronized boolean checkTwinMapping(Cells refCell) {
        boolean _statusFlag = false;
        for (Cells cell : crossMappingTracer) {

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
        for (Double score : scores) {
            Cells choosenCell = new DeadLockResolver().resolver(score);
            deadlockFreeFlagMatrix[choosenCell.indexI][choosenCell.indexJ] = true;
        }
        this.flagMatrix = deadlockFreeFlagMatrix;
    }

    /**
     * Chosen cell of the matrix is stored here
     */
    class Cells {
        String eductName;
        String productName;
        int indexI;
        int indexJ;

    }

    /**
     * Resolves deadlocks if more than one cell clashes with same scores. The
     * decision is then made on the max. stereo score and min. energy score
     */
    class DeadLockResolver {

        private synchronized double getMaxStereo(List<ChooseWinner.Cells> choosenCells) {
            double max = -999;
            for (ChooseWinner.Cells cell : choosenCells) {
                double val = getStereoMatrix().getValue(cell.indexI, cell.indexJ);
                if (val > max) {
                    max = val;
                }
            }
            return max;
        }

        private synchronized double getMinEnergy(List<ChooseWinner.Cells> choosenCells) {
            double min = 999999;
            for (ChooseWinner.Cells cell : choosenCells) {
                double val = getEnergyMatrix().getValue(cell.indexI, cell.indexJ);
                if (val < min) {
                    min = val;
                }
            }
            return min;
        }

        /**
         * Returns best match for the given score
         *
         * @param choosenScore
         * @return
         */
        public synchronized ChooseWinner.Cells resolver(double choosenScore) {
            List<ChooseWinner.Cells> choosenCells = new ArrayList<>();
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    if (flagMatrix[i][j] && choosenScore > MIN_VALUE
                            && similarityMatrix.getValue(i, j) == choosenScore) {
                        ChooseWinner.Cells cells = new ChooseWinner.Cells();
                        cells.indexI = i;
                        cells.indexJ = j;
                        choosenCells.add(cells);
                    }
                }
            }
            double maxStereo = getMaxStereo(choosenCells);
            double minEnergy = getMinEnergy(choosenCells);
            //            System.out.println("maxStereo " + maxStereo);
            //            System.out.println("minEnergy " + minEnergy);
            for (ChooseWinner.Cells cell : choosenCells) {
                double stereoVal = stereoMatrix.getValue(cell.indexI, cell.indexJ);
                double energyVal = energyMatrix.getValue(cell.indexI, cell.indexJ);
                //                System.out.println("stereoVal " + stereoVal);
                //                System.out.println("energyVal " + energyVal);
                if (stereoVal <= maxStereo && energyVal <= minEnergy) {
                    return cell;
                }
            }
            return choosenCells.listIterator().next();
        }
    }
}
