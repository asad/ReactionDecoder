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

import java.io.Serializable;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import java.util.List;

import org.openscience.cdk.interfaces.IAtomContainer;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getTotalFormalCharge;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ReactionIsomorphismHandler implements Serializable {

    private static final long serialVersionUID = 0x1bfce07abac99fL;
    private int rowSize = -1;
    private int colSize = -1;
    private boolean[][] flagSimilarityMatrix = null;
    private boolean[][] flagStereoMatrix = null;
    private boolean isomorphismFlag;
    private Holder matrixHolder;
    private final Holder matrixHolderWithSimilarityCheck;
    private final Holder matrixHolderWithStereoCheck;

    /**
     *
     * @param mHolder
     * @param EdMapOrignal
     * @param PdMapOrignal
     * @throws Exception
     */
    public ReactionIsomorphismHandler(Holder mHolder,
            List<String> EdMapOrignal,
            List<String> PdMapOrignal) throws Exception {
        this.matrixHolder = mHolder;
        this.matrixHolderWithSimilarityCheck = (Holder) mHolder.clone();
        this.matrixHolderWithStereoCheck = (Holder) mHolder.clone();
        this.isomorphismFlag = false;

        rowSize = matrixHolder.getCliqueMatrix().getRowDimension();
        colSize = matrixHolder.getCliqueMatrix().getColumnDimension();

        if (rowSize > 1 && rowSize == colSize) {

            setReactionMappingFlags();

            boolean flag1 = checkSimilarityWithStereo();
            boolean flag2 = checkSimilarityWithFingerprint();

            if (flag1 && !flag2) {
                isomorphismFlag = true;
                this.matrixHolder = this.matrixHolderWithStereoCheck;
            } else if (flag1 && flag2) {
                isomorphismFlag = true;
                this.matrixHolder = this.matrixHolderWithSimilarityCheck;
            }
        }
    }

    private boolean checkSimilarityWithStereo() {

        boolean stSimilarity = false;

        boolean RowT = true;
        boolean ColT = true;

        /*
         * check diagonal elements
         */
        for (int i = 0; i < rowSize; i++) {
            if (!flagStereoMatrix[i][i]) {
                RowT = false;
                break;
            }
        }

        /*
         * check across diagonal elements
         */
        for (int i = rowSize - 1; i >= 0; i--) {
            if (!flagStereoMatrix[i][i]) {
                ColT = false;
                break;
            }
        }

        if (RowT) {
            stSimilarity = true;
            for (int i = 0; i < rowSize; i++) {

                matrixHolderWithStereoCheck.getGraphSimilarityMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCliqueMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getStereoMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCarbonOverlapMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getFragmentMatrix().set(i, i, MAX_VALUE);
                matrixHolderWithStereoCheck.getEnergyMatrix().set(i, i, MAX_VALUE);
            }
        } else if (ColT) {
            stSimilarity = true;
            for (int i = rowSize - 1; i >= 0; i--) {
                matrixHolderWithStereoCheck.getGraphSimilarityMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCliqueMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getStereoMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCarbonOverlapMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getFragmentMatrix().set(i, i, MAX_VALUE);
                matrixHolderWithStereoCheck.getEnergyMatrix().set(i, i, MAX_VALUE);
            }
        }

        return stSimilarity;
    }

    private boolean checkSimilarityWithFingerprint() {
        boolean fpFlag = false;
        boolean RowT = true;
        boolean ColT = true;
        /*
         * check diagonal elements
         */
        for (int i = 0; i < rowSize; i++) {
            if (!flagSimilarityMatrix[i][i]) {
                RowT = false;
                break;
            }
        }

        /*
         * check across diagonal elements
         */
        for (int i = rowSize - 1; i >= 0; i--) {
            if (!flagSimilarityMatrix[i][i]) {
                ColT = false;
                break;
            }
        }

        if (RowT) {
            fpFlag = true;
            for (int i = 0; i < rowSize; i++) {
                matrixHolderWithSimilarityCheck.getGraphSimilarityMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getCliqueMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getStereoMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCarbonOverlapMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getFragmentMatrix().set(i, i, MAX_VALUE);
                matrixHolderWithSimilarityCheck.getEnergyMatrix().set(i, i, MAX_VALUE);
            }
        } else if (ColT) {
            fpFlag = true;
            for (int i = rowSize - 1; i >= 0; i--) {
                matrixHolderWithSimilarityCheck.getGraphSimilarityMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getCliqueMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getStereoMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithStereoCheck.getCarbonOverlapMatrix().setValue(i, i, MIN_VALUE);
                matrixHolderWithSimilarityCheck.getFragmentMatrix().set(i, i, MAX_VALUE);
                matrixHolderWithSimilarityCheck.getEnergyMatrix().set(i, i, MAX_VALUE);
            }
        }
        return fpFlag;
    }

    private void setReactionMappingFlags() throws Exception {
        flagSimilarityMatrix = new boolean[rowSize][colSize];
        flagStereoMatrix = new boolean[rowSize][colSize];

        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                flagSimilarityMatrix[i][j] = false;
                flagStereoMatrix[i][j] = false;
            }
        }

        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {

                IAtomContainer ac1 = matrixHolder.getReactionContainer().getEduct(i);
                IAtomContainer ac2 = matrixHolder.getReactionContainer().getProduct(j);
                //matrix.
                if (matrixHolder.getFPSimilarityMatrix().getValue(i, j) == 1.
                        && getTotalFormalCharge(ac1)
                        == getTotalFormalCharge(ac2)) {

                    try {

                        AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
                        BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);

                        Substructure isomorphism = new Substructure(ac1, ac2, atomMatcher, bondMatcher, false);
                        if (isomorphism.isSubgraph()) {
                            isomorphism.setChemFilters(true, true, true);

                            if (isomorphism.getTanimotoSimilarity() == 1.0) {

                                if (!isomorphism.isStereoMisMatch()) {
                                    flagStereoMatrix[i][j] = true;
                                }
                            }
                        }
                    } catch (Exception ex) {
                        flagStereoMatrix[i][j] = false;
                    }
                    flagSimilarityMatrix[i][j] = true;
                }
            }
        }
    }

    /**
     *
     * @return
     */
    public boolean getIsomorphismFlag() {
        return isomorphismFlag;
    }

    /**
     * @return the matrixHolder
     */
    public Holder getMatrixHolder() {
        return matrixHolder;
    }

    /**
     * @param matrixHolder the matrixHolder to set
     */
    public void setMatrixHolder(Holder matrixHolder) {
        this.matrixHolder = matrixHolder;
    }
}
