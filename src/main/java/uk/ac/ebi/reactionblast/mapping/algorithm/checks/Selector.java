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
package uk.ac.ebi.reactionblast.mapping.algorithm.checks;

import java.io.IOException;
import java.io.Serializable;
import org.openscience.cdk.exception.CDKException;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class Selector implements Serializable {

    /**
     * Refills Matrix with true similarity scores
     *
     * @param orignal
     * @return
     * @throws IOException
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static synchronized Holder modifyMatrix(Holder orignal)
            throws IOException, CDKException, CloneNotSupportedException {
        ReactionContainer reactionStructureInformationContainer = orignal.getReactionContainer();
        Holder localHolder = (Holder) orignal.clone();

        int inputRowSize = orignal.getCliqueMatrix().getRowDimension();
        int inputColSize = orignal.getCliqueMatrix().getColumnDimension();
        /*
         * Flag matrix is assigned
         */
        for (int i = 0; i < inputRowSize; i++) {
            for (int j = 0; j < inputColSize; j++) {
                double totalAtomCount = (double) (reactionStructureInformationContainer.getProduct(j).getAtomCount() + reactionStructureInformationContainer.getEduct(i).getAtomCount());
                double cliqueValue = orignal.getCliqueMatrix().getValue(i, j);
                double simValue = cliqueValue / totalAtomCount;
                /*
                Modify the cloned matrix with real similarity values
                 */
                if (cliqueValue >= 1) {
                    localHolder.getGraphSimilarityMatrix().set(i, j, simValue);
                }
            }
        }

//        System.out.println("**********Modified Min Input Matrix**************");
//        printSimMatrix(localHolder);
        return localHolder;
    }

    int rowSize;
    int colSize;
    boolean[][] flagMatrix;

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     * @throws IOException
     * @throws CDKException
     */
    protected synchronized boolean isMajorSubgraphColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {

        double queryColScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
        if (queryColScore > 0) {
            for (int col = 0; col < colSize; col++) {
                if (flagMatrix[IndexI][col] && col != IndexJ) {
//                    int pSize = mh.getReactionContainer().getProduct(col).getAtomCount();
                    double colCSize = mh.getCliqueMatrix().getValue(IndexI, col);
                    if (queryColScore < colCSize) { //&& colCSize < pSize) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     * @throws IOException
     * @throws CDKException
     */
    protected synchronized boolean isMajorSubgraphRow(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
        double queryRowScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
        if (queryRowScore > 0) {
            for (int row = 0; row < rowSize; row++) {
                if (flagMatrix[row][IndexJ] && row != IndexI) {
                    int eSize = mh.getReactionContainer().getEduct(row).getAtomCount();
                    double rowRSize = mh.getCliqueMatrix().getValue(row, IndexJ);
                    if (queryRowScore < rowRSize) {//&& rowRSize == eSize) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     * @throws IOException
     * @throws CDKException
     */
    protected synchronized boolean isMinEnergyColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
        //        System.out.println("CHECKING ENERGY\nI " + EdMap.get(IndexI) + " J " + PdMap.get(IndexJ));
        double refEnergy = mh.getEnergyMatrix().getValue(IndexI, IndexJ);
        if (mh.getCliqueMatrix().getValue(IndexI, IndexJ) > 0) {
            for (int col = 0; col < colSize; col++) {
                if (flagMatrix[IndexI][col] && col != IndexJ) {
                    double colEnergy = mh.getEnergyMatrix().getValue(IndexI, col);
                    if (refEnergy > 0 && colEnergy > 0 && refEnergy > colEnergy) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     * @throws IOException
     * @throws CDKException
     */
    protected synchronized boolean isMinEnergyRow(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
        double refEnergy = mh.getEnergyMatrix().getValue(IndexI, IndexJ);
        if (mh.getCliqueMatrix().getValue(IndexI, IndexJ) > 0) {
            for (int row = 0; row < rowSize; row++) {
                if (flagMatrix[row][IndexJ] && row != IndexI) {
                    double rowEnergy = mh.getEnergyMatrix().getValue(row, IndexJ);
                    if (rowEnergy > 0 && refEnergy > 0 && refEnergy > rowEnergy) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     * @throws IOException
     * @throws CDKException
     */
    protected synchronized boolean isMinorSubgraphColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
        boolean flag = true;
        double queryColScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
        if (queryColScore > 0.) {
            for (int col = 0; col < colSize; col++) {
                if (flagMatrix[IndexI][col] && col != IndexJ) {
                    double targetColScore = mh.getCliqueMatrix().getValue(IndexI, col);
                    if (queryColScore > targetColScore && targetColScore > 0.) {
                        return true;
                    }
                }
            }
        }
        return flag;
    }

    /**
     *
     * @param mh
     * @param IndexI
     * @param IndexJ
     * @return
     */
    protected synchronized boolean isMinorSubgraphRow(Holder mh, int IndexI, int IndexJ) {
        boolean flag = true;
        double queryRowScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
        if (queryRowScore > 0.) {
            for (int row = 0; row < rowSize; row++) {
                if (flagMatrix[row][IndexJ] && row != IndexI) {
                    double targetRowScore = mh.getCliqueMatrix().getValue(row, IndexJ);
                    if (queryRowScore > targetRowScore && targetRowScore > 0.) {
                        return false;
                    }
                }
            }
        }
        return flag;
    }

    /**
     *
     * @param similarityMatrix
     * @param IndexI
     * @param IndexJ
     * @return
     */
    protected synchronized boolean isMinorSubgraphColumn(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
        boolean flag = true;
        double queryColScore = similarityMatrix.getValue(IndexI, IndexJ);
        if (queryColScore > 0.) {
            for (int col = 0; col < colSize; col++) {
                if (flagMatrix[IndexI][col] && col != IndexJ) {
                    double targetColScore = similarityMatrix.getValue(IndexI, col);
                    if (queryColScore > targetColScore && targetColScore > 0.) {
                        return true;
                    }
                }
            }
        }
        return flag;
    }

    /**
     *
     * @param similarityMatrix
     * @param IndexI
     * @param IndexJ
     * @return
     */
    protected synchronized boolean isMinorSubgraphRow(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
        boolean flag = true;
        double queryRowScore = similarityMatrix.getValue(IndexI, IndexJ);
        if (queryRowScore > 0.) {
            for (int row = 0; row < rowSize; row++) {
                if (flagMatrix[row][IndexJ] && row != IndexI) {
                    double targetRowScore = similarityMatrix.getValue(row, IndexJ);
                    if (queryRowScore > targetRowScore && targetRowScore > 0.) {
                        return false;
                    }
                }
            }
        }
        return flag;
    }

    /**
     *
     * @param similarityMatrix
     * @param IndexI
     * @param IndexJ
     * @return
     */
    protected synchronized boolean isMajorSubgraphColumn(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {

        double queryColumnscore = similarityMatrix.getValue(IndexI, IndexJ);
        if (queryColumnscore > 0) {
            for (int col = 0; col < colSize; col++) {
                if (flagMatrix[IndexI][col] && col != IndexJ) {
                    double colCSize = similarityMatrix.getValue(IndexI, col);
                    if (queryColumnscore < colCSize) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     *
     * @param similarityMatrix
     * @param IndexI
     * @param IndexJ
     * @return
     */
    protected synchronized boolean isMajorSubgraphRow(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
        double queryRowScore = similarityMatrix.getValue(IndexI, IndexJ);
        if (queryRowScore > 0) {
            for (int row = 0; row < rowSize; row++) {
                if (flagMatrix[row][IndexJ] && row != IndexI) {
                    double rowRSize = similarityMatrix.getValue(row, IndexJ);
                    if (queryRowScore < rowRSize) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}
