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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;
import java.util.logging.Logger;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getTotalFormalCharge;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogens;
import static org.openscience.smsd.ExtAtomContainerManipulator.Utility.isMatch;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import com.bioinceptionlabs.reactionblast.mapping.ReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.SmsdReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.legacy.EBIMatrix;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
interface IResult {
    abstract Holder getUpdatedHolder();
    abstract boolean isSubAndCompleteMatchFlag();
}

/**
 * Consolidated mapping checks for the algorithm package.
 * Merges: Selector, ChooseWinner, MaxSelection, MinSelection,
 *         ReactionIsomorphismHandler, RuleBasedMappingHandler
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class MappingChecks {

    private static final ReactionMappingEngine MAPPING_ENGINE
            = SmsdReactionMappingEngine.getInstance();

    private MappingChecks() { /* utility class */ }

    // ========== Selector (abstract base for ChooseWinner/MaxSelection/MinSelection) ==========

    public static abstract class Selector implements Serializable {

        public static Holder modifyMatrix(Holder orignal)
                throws IOException, CDKException, CloneNotSupportedException {
            ReactionContainer reactionStructureInformationContainer = orignal.getReactionContainer();
            Holder localHolder = (Holder) orignal.clone();
            int inputRowSize = orignal.getCliqueMatrix().getRowDimension();
            int inputColSize = orignal.getCliqueMatrix().getColumnDimension();
            for (int i = 0; i < inputRowSize; i++) {
                for (int j = 0; j < inputColSize; j++) {
                    double totalAtomCount = (double) (reactionStructureInformationContainer.getProduct(j).getAtomCount()
                            + reactionStructureInformationContainer.getEduct(i).getAtomCount());
                    double cliqueValue = orignal.getCliqueMatrix().getValue(i, j);
                    double simValue = cliqueValue / totalAtomCount;
                    if (cliqueValue >= 1) {
                        localHolder.getGraphSimilarityMatrix().set(i, j, simValue);
                    }
                }
            }
            return localHolder;
        }

        int rowSize;
        int colSize;
        boolean[][] flagMatrix;

        protected boolean isMajorSubgraphColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
            double queryColScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
            if (queryColScore > 0) {
                for (int col = 0; col < colSize; col++) {
                    if (flagMatrix[IndexI][col] && col != IndexJ) {
                        double colCSize = mh.getCliqueMatrix().getValue(IndexI, col);
                        if (queryColScore < colCSize) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        protected boolean isMajorSubgraphRow(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
            double queryRowScore = mh.getCliqueMatrix().getValue(IndexI, IndexJ);
            if (queryRowScore > 0) {
                for (int row = 0; row < rowSize; row++) {
                    if (flagMatrix[row][IndexJ] && row != IndexI) {
                        int eSize = mh.getReactionContainer().getEduct(row).getAtomCount();
                        double rowRSize = mh.getCliqueMatrix().getValue(row, IndexJ);
                        if (queryRowScore < rowRSize) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        protected boolean isMinEnergyColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
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

        protected boolean isMinEnergyRow(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
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

        protected boolean isMinorSubgraphColumn(Holder mh, int IndexI, int IndexJ) throws IOException, CDKException {
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

        protected boolean isMinorSubgraphRow(Holder mh, int IndexI, int IndexJ) {
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

        protected boolean isMajorSubgraphColumn(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
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

        protected boolean isMajorSubgraphRow(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
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

        protected boolean isMinorSubgraphColumn(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
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

        protected boolean isMinorSubgraphRow(EBIMatrix similarityMatrix, int IndexI, int IndexJ) {
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
    }

    // ========== ChooseWinner ==========

    public static class ChooseWinner extends Selector implements Serializable {

        private static final long serialVersionUID = 0x296558709L;
        private static final ILoggingTool LOGGER = createLoggingTool(ChooseWinner.class);
        private EBIMatrix stereoMatrix;
        private EBIMatrix energyMatrix;
        private EBIMatrix carbonOverlapMatrix;
        private EBIMatrix similarityMatrix = null;
        private List<Cell> crossMappingTracer = null;
        private Map<Integer, IAtomContainer> educts = null;
        private Map<Integer, IAtomContainer> products = null;

        public ChooseWinner(List<String> eductNameList, List<String> productNameList) {
            this.rowSize = eductNameList.size();
            this.colSize = productNameList.size();
            this.flagMatrix = new boolean[rowSize][colSize];
        }

        public EBIMatrix getStereoMatrix() { return stereoMatrix; }
        public void setStereoMatrix(EBIMatrix stereoMatrix) { this.stereoMatrix = stereoMatrix; }
        public EBIMatrix getEnergyMatrix() { return energyMatrix; }
        public EBIMatrix getCarbonOverlapMatrix() { return carbonOverlapMatrix; }
        public void setEnergyMatrix(EBIMatrix energyMatrix) { this.energyMatrix = energyMatrix; }
        public void setCarbonOverlapMatrix(EBIMatrix carbonOverlapMatrix) { this.carbonOverlapMatrix = carbonOverlapMatrix; }

        public void searchWinners(Map<Integer, IAtomContainer> eductMap,
                Map<Integer, IAtomContainer> productMap, Holder mHolder) {
            initFlagMatrix();
            this.educts = eductMap;
            this.products = productMap;
            this.similarityMatrix = mHolder.getGraphSimilarityMatrix();
            this.setStereoMatrix(mHolder.getStereoMatrix());
            this.setEnergyMatrix(mHolder.getEnergyMatrix());
            this.setCarbonOverlapMatrix(mHolder.getCarbonOverlapMatrix());
            this.crossMappingTracer = new ArrayList<>();
            boolean isMappingFesiable = checkStatusFlag();
            List<Double> scores = new ArrayList<>();
            if (isMappingFesiable) {
                double maximumSimilarity = 0.0;
                boolean maxValueI;
                boolean maxValueJ;
                for (int i = 0; i < rowSize; i++) {
                    for (int j = 0; j < colSize; j++) {
                        double similarity = similarityMatrix.getValue(i, j);
                        if (similarity > MIN_VALUE) {
                            maxValueI = isMajorSubgraphRow(similarityMatrix, i, j);
                            maxValueJ = isMajorSubgraphColumn(similarityMatrix, i, j);
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
            resolveDeadLocks(scores);
            setWinOverFlags();
        }

        public boolean getFlag() {
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    if (this.flagMatrix[i][j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        public boolean[][] getFlagMatrix() { return this.flagMatrix; }

        private void initFlagMatrix() {
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    this.flagMatrix[i][j] = false;
                }
            }
        }

        private void initFlagMatrix(boolean[][] boolMatrix, int rowSize, int colSize) {
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    boolMatrix[i][j] = false;
                }
            }
        }

        private boolean checkStatusFlag() {
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    if (similarityMatrix.getValue(i, j) > MIN_VALUE) {
                        return true;
                    }
                }
            }
            return false;
        }

        private void setWinOverFlags() {
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

        private boolean checkTwinMapping(Cell refCell) {
            boolean _statusFlag = false;
            for (Cell cell : crossMappingTracer) {
                if (cell.eductName.equals(refCell.eductName)
                        && cell.productName.equals(refCell.productName)) {
                    if (cell.indexI == refCell.indexI || cell.indexJ == refCell.indexJ) {
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

        private void resolveDeadLocks(List<Double> scores) {
            boolean[][] deadlockFreeFlagMatrix = new boolean[rowSize][colSize];
            initFlagMatrix(deadlockFreeFlagMatrix, rowSize, colSize);
            scores.stream().map((score) -> new DeadLockResolver().resolver(score)).forEach((choosenCell) -> {
                deadlockFreeFlagMatrix[choosenCell.indexI][choosenCell.indexJ] = true;
            });
            this.flagMatrix = deadlockFreeFlagMatrix;
        }

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

        class DeadLockResolver {
            private double getMaxStereo(List<ChooseWinner.Cell> choosenCells) {
                double max = 0.0;
                for (ChooseWinner.Cell cell : choosenCells) {
                    double val = getStereoMatrix().getValue(cell.indexI, cell.indexJ);
                    if (val > max) { max = val; }
                }
                return max;
            }

            private double getMinEnergy(List<ChooseWinner.Cell> choosenCells) {
                double min = Double.MAX_VALUE;
                for (ChooseWinner.Cell cell : choosenCells) {
                    double val = getEnergyMatrix().getValue(cell.indexI, cell.indexJ);
                    if (val < min) { min = val; }
                }
                return min;
            }

            private double getMaxCarbonOverlap(List<ChooseWinner.Cell> choosenCells) {
                double max = Double.MIN_VALUE;
                for (ChooseWinner.Cell cell : choosenCells) {
                    double val = getCarbonOverlapMatrix().getValue(cell.indexI, cell.indexJ);
                    if (val > max) { max = val; }
                }
                return max;
            }

            public ChooseWinner.Cell resolver(double choosenScore) {
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
                for (var cell : choosenCells) {
                    double stereoVal = stereoMatrix.getValue(cell.indexI, cell.indexJ);
                    double energyVal = energyMatrix.getValue(cell.indexI, cell.indexJ);
                    double carbonOverlap = carbonOverlapMatrix.getValue(cell.indexI, cell.indexJ);
                    if (stereoVal >= maxStereo
                            && energyVal <= minEnergy
                            && (maxCarbonOverlap == carbonOverlap || maxCarbonOverlap == 0)) {
                        winner = cell;
                        break;
                    }
                }
                return winner;
            }
        }
    }

    // ========== MaxSelection ==========

    public static class MaxSelection extends Selector implements IResult {

        private final static ILoggingTool LOGGER = createLoggingTool(MaxSelection.class);
        private static final long serialVersionUID = 0x192aa60a59L;
        private final Holder mHolder;
        private final Holder updatedHolder;
        private boolean SubAndCompleteFlag;

        public MaxSelection(Holder mHolder, List<String> EdMapOrignal, List<String> PdMapOrignal)
                throws IOException, Exception {
            this.mHolder = mHolder;
            this.updatedHolder = (Holder) mHolder.clone();
            rowSize = mHolder.getCliqueMatrix().getRowDimension();
            colSize = mHolder.getCliqueMatrix().getColumnDimension();
            this.flagMatrix = new boolean[rowSize][colSize];
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    this.flagMatrix[i][j] = false;
                }
            }
            try {
                PhaseOneMatcher();
                SubAndCompleteFlag = PhaseTwoMatcher();
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        private void PhaseOneMatcher() throws IOException, CDKException {
            boolean maxValueI;
            boolean maxValueJ;
            for (int i = 0; i < rowSize; i++) {
                IAtomContainer ac1 = this.mHolder.getReactionContainer().getEduct(i);
                for (int j = 0; j < colSize; j++) {
                    IAtomContainer ac2 = this.mHolder.getReactionContainer().getProduct(j);
                    double cliqueVal = mHolder.getCliqueMatrix().getValue(i, j);
                    double simVal = mHolder.getGraphSimilarityMatrix().getValue(i, j);
                    if (cliqueVal > 0.0) {
                        maxValueI = isMajorSubgraphRow(mHolder, i, j);
                        maxValueJ = isMajorSubgraphColumn(mHolder, i, j);
                        if (maxValueI || maxValueJ) {
                            if (simVal > 0 && ac1 != null && ac2 != null) {
                                int eMolSize = ac1.getAtomCount();
                                int pMolSize = ac2.getAtomCount();
                                if (eMolSize == 1 && pMolSize == 1
                                        && (ac1.atoms().iterator().next().getSymbol()
                                                .equals(ac2.atoms().iterator().next().getSymbol()))) {
                                    this.flagMatrix[i][j] = true;
                                } else if (eMolSize > 1 && pMolSize > 1
                                        && ac1.getAtomCount() <= ac2.getAtomCount()
                                        && isMCSSubgraph(ac1, mHolder.getCliqueMatrix().getValue(i, j))) {
                                    this.flagMatrix[i][j] = true;
                                } else if (eMolSize > 1 && pMolSize > 1
                                        && ac2.getAtomCount() <= ac1.getAtomCount()
                                        && isMCSSubgraph(ac2, mHolder.getCliqueMatrix().getValue(i, j))) {
                                    this.flagMatrix[i][j] = true;
                                } else {
                                    this.flagMatrix[i][j] = false;
                                    this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, MIN_VALUE);
                                    this.updatedHolder.getCliqueMatrix().setValue(i, j, MIN_VALUE);
                                    this.updatedHolder.getStereoMatrix().setValue(i, j, MIN_VALUE);
                                    this.updatedHolder.getCarbonOverlapMatrix().setValue(i, j, MIN_VALUE);
                                    this.updatedHolder.getFragmentMatrix().setValue(i, j, MAX_VALUE);
                                    this.updatedHolder.getEnergyMatrix().setValue(i, j, MAX_VALUE);
                                    this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, MIN_VALUE);
                                }
                            }
                        } else {
                            this.flagMatrix[i][j] = false;
                            this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getCliqueMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getStereoMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getCarbonOverlapMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getFragmentMatrix().setValue(i, j, MAX_VALUE);
                            this.updatedHolder.getEnergyMatrix().setValue(i, j, MAX_VALUE);
                            this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, MIN_VALUE);
                        }
                    }
                }
            }
        }

        private boolean PhaseTwoMatcher() throws IOException, CDKException {
            boolean flag = false;
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    double cliqueVal = this.updatedHolder.getCliqueMatrix().getValue(i, j);
                    double simVal = this.updatedHolder.getGraphSimilarityMatrix().getValue(i, j);
                    if (simVal != 1.0 && cliqueVal >= 1) {
                        if (flagMatrix[i][j]) {
                            flag = true;
                        } else {
                            this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getCliqueMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getStereoMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getCarbonOverlapMatrix().setValue(i, j, MIN_VALUE);
                            this.updatedHolder.getFragmentMatrix().setValue(i, j, MAX_VALUE);
                            this.updatedHolder.getEnergyMatrix().setValue(i, j, MAX_VALUE);
                            this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, MIN_VALUE);
                        }
                    }
                }
            }
            return flag;
        }

        private boolean isMCSSubgraph(IAtomContainer educt, double mcsSize) throws CDKException {
            return educt.getAtomCount() == mcsSize;
        }

        @Override
        public boolean isSubAndCompleteMatchFlag() { return SubAndCompleteFlag; }

        @Override
        public Holder getUpdatedHolder() { return updatedHolder; }
    }

    // ========== MinSelection ==========

    public static class MinSelection extends Selector implements IResult {

        private final static ILoggingTool LOGGER = createLoggingTool(MinSelection.class);
        private static final long serialVersionUID = 1908987778L;
        private final Holder mHolder;
        private final Holder updatedHolder;
        private boolean isSubstructure;

        public MinSelection(Holder mHolder, List<String> EdMapOrignal, List<String> PdMapOrignal)
                throws IOException, CloneNotSupportedException {
            this.mHolder = mHolder;
            this.updatedHolder = (Holder) mHolder.clone();
            rowSize = mHolder.getCliqueMatrix().getRowDimension();
            colSize = mHolder.getCliqueMatrix().getColumnDimension();
            try {
                PhaseOneMatcher();
                isSubstructure = PhaseTwoMatcher();
            } catch (IOException | CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        private void PhaseOneMatcher() throws IOException, CDKException {
            this.flagMatrix = new boolean[rowSize][colSize];
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    this.flagMatrix[i][j] = false;
                }
            }
            for (int i = 0; i < rowSize; i++) {
                IAtomContainer ac1 = this.mHolder.getReactionContainer().getEduct(i);
                for (int j = 0; j < colSize; j++) {
                    IAtomContainer ac2 = this.mHolder.getReactionContainer().getProduct(j);
                    double cliqueVal = mHolder.getCliqueMatrix().getValue(i, j);
                    if (cliqueVal > 0.0) {
                        double simVal = mHolder.getGraphSimilarityMatrix().getValue(i, j);
                        if (simVal > 0 && ac1 != null && ac2 != null) {
                            int eMolSize = ac1.getAtomCount();
                            int pMolSize = ac2.getAtomCount();
                            if (eMolSize == 1 && pMolSize == 1
                                    && (ac1.atoms().iterator().next().getSymbol().equals(ac2.atoms().iterator().next().getSymbol()))) {
                                this.flagMatrix[i][j] = true;
                            } else if (eMolSize > 1 && pMolSize > 1 && ac1.getAtomCount() <= ac2.getAtomCount()
                                    && isMCSSubgraph(ac1, mHolder.getCliqueMatrix().getValue(i, j))) {
                                this.flagMatrix[i][j] = true;
                            } else if (eMolSize > 1 && pMolSize > 1 && ac2.getAtomCount() <= ac1.getAtomCount()
                                    && isMCSSubgraph(ac2, mHolder.getCliqueMatrix().getValue(i, j))) {
                                this.flagMatrix[i][j] = true;
                            } else {
                                this.flagMatrix[i][j] = false;
                                this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, MIN_VALUE);
                                this.updatedHolder.getCliqueMatrix().setValue(i, j, MIN_VALUE);
                                this.updatedHolder.getStereoMatrix().setValue(i, j, MIN_VALUE);
                                this.updatedHolder.getCarbonOverlapMatrix().setValue(i, j, MIN_VALUE);
                                this.updatedHolder.getFragmentMatrix().setValue(i, j, MAX_VALUE);
                                this.updatedHolder.getEnergyMatrix().setValue(i, j, MAX_VALUE);
                                this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, MIN_VALUE);
                            }
                        }
                    }
                }
            }
        }

        private boolean PhaseTwoMatcher() throws IOException, CDKException {
            boolean flag = false;
            for (int i = 0; i < rowSize; i++) {
                IAtomContainer ac1 = this.updatedHolder.getReactionContainer().getEduct(i);
                for (int j = 0; j < colSize; j++) {
                    IAtomContainer ac2 = this.updatedHolder.getReactionContainer().getProduct(j);
                    double graphSim = this.updatedHolder.getGraphSimilarityMatrix().getValue(i, j);
                    if (flagMatrix[i][j]) {
                        double eMolSize = ac1.getAtomCount();
                        double pMolSize = ac2.getAtomCount();
                        if (eMolSize > 0 && pMolSize > 0) {
                            if (graphSim > 0.) {
                                boolean isMinorSubgraphColumn = isMinorSubgraphColumn(updatedHolder, i, j);
                                boolean isMinorSubgraphRow = isMinorSubgraphRow(updatedHolder, i, j);
                                if (isMinorSubgraphColumn && isMinorSubgraphRow) {
                                    double updatedGraphSimScore = 1.01 - (this.mHolder.getGraphSimilarityMatrix().getValue(i, j));
                                    double updatedFPSimScore = 1.01 - (this.mHolder.getFPSimilarityMatrix().getValue(i, j));
                                    this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, updatedGraphSimScore);
                                    this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, updatedFPSimScore);
                                    flag = true;
                                } else {
                                    this.updatedHolder.getGraphSimilarityMatrix().setValue(i, j, MIN_VALUE);
                                    this.updatedHolder.getFPSimilarityMatrix().setValue(i, j, MIN_VALUE);
                                }
                            }
                        }
                    }
                }
            }
            return flag;
        }

        @Override
        public boolean isSubAndCompleteMatchFlag() { return isSubstructure; }

        private boolean isMCSSubgraph(IAtomContainer educt, double mcsSize) throws CDKException {
            return educt.getAtomCount() == (int) mcsSize;
        }

        @Override
        public Holder getUpdatedHolder() { return updatedHolder; }

        static void printSimMatrix(Holder mh) {
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            LOGGER.debug("");
            LOGGER.debug("********* MATRIX **********");
            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;
                LOGGER.debug("Similarity Matrix");
                LOGGER.debug("\t\t");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
                }
                LOGGER.debug("");
                double val = 0;
                for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                    for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                        val = mh.getGraphSimilarityMatrix().getValue(i, j);
                        result = format.format(val);
                        LOGGER.debug("   " + result);
                    }
                    LOGGER.debug("");
                }
            } catch (IOException | CDKException e) {
                LOGGER.debug(" Parser Error: ");
            }
            LOGGER.debug("");
        }

        static void printCliqueMatrix(Holder mh) {
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            LOGGER.debug("");
            LOGGER.debug("********* MATRIX **********");
            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;
                LOGGER.debug("Clique Matrix");
                LOGGER.debug("\t\t");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
                }
                LOGGER.debug("");
                double val = 0;
                for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                    for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                        val = mh.getCliqueMatrix().getValue(i, j);
                        result = format.format(val);
                        LOGGER.debug("   " + result);
                    }
                    LOGGER.debug("");
                }
            } catch (IOException | CDKException e) {
                LOGGER.debug(" Parser Error: ");
            }
            LOGGER.debug("");
        }

        static void printFLAGMatrix(Holder mh, boolean[][] flagMatrix) {
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            LOGGER.debug("");
            LOGGER.debug("********* MATRIX **********");
            try {
                String result;
                LOGGER.debug("Flag Matrix");
                LOGGER.debug("\t\t");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
                }
                LOGGER.debug("");
                boolean val;
                for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                    LOGGER.debug(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                    for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                        LOGGER.debug("   " + flagMatrix[i][j]);
                    }
                    LOGGER.debug("");
                }
            } catch (IOException | CDKException e) {
                LOGGER.error("Parser Error", e.getMessage());
            }
            LOGGER.debug("");
        }
    }

    // ========== ReactionIsomorphismHandler ==========

    public static class ReactionIsomorphismHandler implements Serializable {

        private static final long serialVersionUID = 0x1bfce07abac99fL;
        private int rowSize = -1;
        private int colSize = -1;
        private boolean[][] flagSimilarityMatrix = null;
        private boolean[][] flagStereoMatrix = null;
        private boolean isomorphismFlag;
        private Holder matrixHolder;
        private final Holder matrixHolderWithSimilarityCheck;
        private final Holder matrixHolderWithStereoCheck;

        public ReactionIsomorphismHandler(Holder mHolder, List<String> EdMapOrignal, List<String> PdMapOrignal)
                throws Exception {
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
            for (int i = 0; i < rowSize; i++) {
                if (!flagStereoMatrix[i][i]) { RowT = false; break; }
            }
            for (int i = rowSize - 1; i >= 0; i--) {
                if (!flagStereoMatrix[i][i]) { ColT = false; break; }
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
            for (int i = 0; i < rowSize; i++) {
                if (!flagSimilarityMatrix[i][i]) { RowT = false; break; }
            }
            for (int i = rowSize - 1; i >= 0; i--) {
                if (!flagSimilarityMatrix[i][i]) { ColT = false; break; }
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
                    if (matrixHolder.getFPSimilarityMatrix().getValue(i, j) == 1.
                            && getTotalFormalCharge(ac1) == getTotalFormalCharge(ac2)) {
                        try {
                            AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
                            BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);
                            BaseMapping isomorphism = MAPPING_ENGINE.findSubstructure(
                                    ac1, ac2, atomMatcher, bondMatcher, false);
                            if (isomorphism.isSubgraph()) {
                                MAPPING_ENGINE.applyDefaultFilters(isomorphism);
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

        public boolean getIsomorphismFlag() { return isomorphismFlag; }
        public Holder getMatrixHolder() { return matrixHolder; }
        public void setMatrixHolder(Holder matrixHolder) { this.matrixHolder = matrixHolder; }
    }

    // ========== RuleBasedMappingHandler ==========

    public static final class RuleBasedMappingHandler implements Serializable {

        private static final long serialVersionUID = 88765671L;
        private final static ILoggingTool LOGGER = createLoggingTool(RuleBasedMappingHandler.class);
        private boolean ruleMatched = false;
        private Holder matrixHolder;
        private Holder matrixHolderClone;
        private final Map<Integer, Integer> matchedRowColoumn;
        private IAtomContainer smartsWater;
        private IAtomContainer smartsPhosphate;
        private IAtomContainer smartsSulphate;
        private IAtomContainer smartsL_Glutamate;
        private IAtomContainer smartsL_Glutamine;
        private IAtomContainer smartsL_Glutamine_clipped;
        private IAtomContainer smartsL_Glutamate_clipped;
        private IAtomContainer smartsTwoOxoglutarate;
        private IAtomContainer smartsD_Glutamate;
        private IAtomContainer smartsAcetate;
        private IAtomContainer smartsATP;
        private IAtomContainer smartsADP;
        private IAtomContainer smartsCoA;
        private IAtomContainer smartsAcetyl_CoA;
        private IAtomContainer smartsC00003;
        private IAtomContainer smartsC00006;
        private IAtomContainer smartsC00004;
        private IAtomContainer smartsC00005;
        private IAtomContainer smartsPyruvate;
        private IAtomContainer smartsAlanine;
        private IAtomContainer smartsNRule;
        private IAtomContainer smartsCRule;
        private IAtomContainer smartsDoublePhosphate;
        private IAtomContainer smartsC04666Rule;
        private IAtomContainer smartsC04916Rule;

        public RuleBasedMappingHandler(Holder matrixHolder, List<String> EdMapOrignal, List<String> PdMapOrignal)
                throws CDKException, IOException {
            setRulesSmiles();
            this.matrixHolder = matrixHolder;
            this.matchedRowColoumn = new HashMap<>();
            setRuleMatched(false);

            int smallestMatchedReactant = Integer.MAX_VALUE;
            int smallestMatchedProduct = Integer.MAX_VALUE;
            for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                IAtomContainer ac1 = new AtomContainer(this.matrixHolder.getReactionContainer().getEduct(i));
                ac1 = removeHydrogens(ac1);
                if (ac1.getAtomCount() >= getSmartsPhosphate().getAtomCount()
                        || ac1.getAtomCount() >= getSmartsSulphate().getAtomCount()) {
                    if (isMatch(getSmartsPhosphate(), ac1, false) || isMatch(getSmartsSulphate(), ac1, false)) {
                        if (smallestMatchedReactant > ac1.getAtomCount()) {
                            smallestMatchedReactant = ac1.getAtomCount();
                        }
                    }
                }
            }
            for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                IAtomContainer ac2 = new AtomContainer(this.matrixHolder.getReactionContainer().getProduct(j));
                ac2 = removeHydrogens(ac2);
                if (ac2.getAtomCount() >= getSmartsPhosphate().getAtomCount()
                        || ac2.getAtomCount() >= getSmartsSulphate().getAtomCount()) {
                    if (isMatch(getSmartsPhosphate(), ac2, false) || isMatch(getSmartsSulphate(), ac2, false)) {
                        if (smallestMatchedProduct > ac2.getAtomCount()) {
                            smallestMatchedProduct = ac2.getAtomCount();
                        }
                    }
                }
            }

            boolean phosphate_changed = phosphate_cleaved(this.matrixHolder.getReactionContainer().getEducts(),
                    this.matrixHolder.getReactionContainer().getProducts());

            try {
                for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                    IAtomContainer educt = this.matrixHolder.getReactionContainer().getEduct(i);
                    IAtomContainer ac1 = new AtomContainer(educt);
                    ac1 = removeHydrogens(ac1);
                    for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                        IAtomContainer product = this.matrixHolder.getReactionContainer().getProduct(j);
                        IAtomContainer ac2 = new AtomContainer(product);
                        ac2 = removeHydrogens(ac2);
                        if (this.matrixHolder.getCliqueMatrix().getValue(i, j) == 0) { continue; }
                        if (phosphate_changed && ac1.getAtomCount() == 1
                                && isMatch(getSmartsWater(), ac1, false) && isMatch(getSmartsPhosphate(), ac2, false)
                                && !isMatch(getSmartsDoublePhosphate(), ac2, false)
                                && ac2.getAtomCount() == smallestMatchedProduct) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        }
                        if (phosphate_changed && ac2.getAtomCount() == 1
                                && isMatch(getSmartsWater(), ac2, false) && isMatch(getSmartsPhosphate(), ac1, false)
                                && !isMatch(getSmartsDoublePhosphate(), ac1, false)
                                && ac1.getAtomCount() == smallestMatchedReactant) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        }
                        if (ac1.getAtomCount() == 1 && isMatch(getSmartsWater(), ac1, false)
                                && isMatch(getSmartsSulphate(), ac2, false) && ac2.getAtomCount() == smallestMatchedProduct) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if (ac2.getAtomCount() == 1 && isMatch(getSmartsWater(), ac2, false)
                                && isMatch(getSmartsSulphate(), ac1, false) && ac1.getAtomCount() == smallestMatchedReactant) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if (isMatch(getSmartsC04666Rule(), ac1, false) && isMatch(getSmartsC04916Rule(), ac2, false)
                                || (isMatch(getSmartsC04916Rule(), ac1, false) && isMatch(getSmartsC04666Rule(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                                && isMatch(getSmartsGlutamate(), ac1, false) && isMatch(getSmartsGlutamine(), ac2, false))
                                || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                                && isMatch(getSmartsGlutamine(), ac1, false) && isMatch(getSmartsGlutamate(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                                && isMatch(getSmartsGlutamateClipped(), ac1, false) && isMatch(getSmartsGlutamineClipped(), ac2, false))
                                || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                                && isMatch(getSmartsGlutamineClipped(), ac1, false) && isMatch(getSmartsGlutamateClipped(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac2.getAtomCount() == 10 && ac1.getAtomCount() == 10
                                && isMatch(getSmartsTwoOxoglutarate(), ac2, false) && isMatch(getSmartsD_Glutamate(), ac1, false))
                                || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                                && isMatch(getSmartsTwoOxoglutarate(), ac1, false) && isMatch(getSmartsD_Glutamate(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == 1 && isMatch(getSmartsWater(), ac1, false)
                                && ac2.getAtomCount() == getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac2, false))
                                || (ac2.getAtomCount() == 1 && isMatch(getSmartsWater(), ac2, false)
                                && ac1.getAtomCount() == getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac1, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == getSmartsATP().getAtomCount() && isMatch(getSmartsATP(), ac1, false)
                                && isMatch(getSmartsADP(), ac2, false))
                                || (ac1.getAtomCount() == getSmartsADP().getAtomCount() && isMatch(getSmartsADP(), ac1, false)
                                && isMatch(getSmartsATP(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == getSmartsCoA().getAtomCount() && isMatch(getSmartsCoA(), ac1, false)
                                && isMatch(getSmartsAcetyl_CoA(), ac2, false))
                                || (ac1.getAtomCount() == getSmartsAcetyl_CoA().getAtomCount() && isMatch(getSmartsAcetyl_CoA(), ac1, false)
                                && isMatch(getSmartsCoA(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == getSmartsC00003().getAtomCount() && isMatch(getSmartsC00003(), ac1, false)
                                && isMatch(getSmartsC00006(), ac2, false))
                                || (ac1.getAtomCount() == getSmartsC00006().getAtomCount() && isMatch(getSmartsC00006(), ac1, false)
                                && isMatch(getSmartsC00003(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == getSmartsC00004().getAtomCount() && isMatch(getSmartsC00004(), ac1, false)
                                && isMatch(getSmartsC00005(), ac2, false))
                                || (ac1.getAtomCount() == getSmartsC00005().getAtomCount() && isMatch(getSmartsC00005(), ac1, false)
                                && isMatch(getSmartsC00004(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if ((ac1.getAtomCount() == getSmartsPyruvate().getAtomCount() && isMatch(getSmartsPyruvate(), ac1, false)
                                && isMatch(getSmartsAlanine(), ac2, false))
                                || (ac1.getAtomCount() == getSmartsAlanine().getAtomCount() && isMatch(getSmartsAlanine(), ac1, false)
                                && isMatch(getSmartsPyruvate(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        } else if (isMatch(getSmartsNRule(), ac1, false) && isMatch(getSmartsCRule(), ac2, false)
                                || (isMatch(getSmartsCRule(), ac1, false) && isMatch(getSmartsNRule(), ac2, false))) {
                            setRuleMatched(true); matchedRowColoumn.put(i, j);
                        }
                    }
                }
            } catch (IOException | CDKException ex) {
                LOGGER.error(WARNING, "Error in Matching Rules", ex);
            }
            if (this.isMatchFound()) {
                try {
                    this.matrixHolderClone = (Holder) matrixHolder.clone();
                } catch (CloneNotSupportedException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                    for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                        if (this.matchedRowColoumn.containsKey(i) && this.matchedRowColoumn.get(i) == j) {
                            matrixHolderClone.getGraphSimilarityMatrix().setValue(i, j, 1.00);
                        } else {
                            matrixHolderClone.getGraphSimilarityMatrix().setValue(i, j, Double.MIN_VALUE);
                            matrixHolderClone.getCliqueMatrix().setValue(i, j, Double.MIN_VALUE);
                            matrixHolderClone.getStereoMatrix().setValue(i, j, Double.MIN_VALUE);
                            matrixHolderClone.getCarbonOverlapMatrix().setValue(i, j, Double.MIN_VALUE);
                            matrixHolderClone.getFragmentMatrix().setValue(i, j, Double.MAX_VALUE);
                            matrixHolderClone.getEnergyMatrix().setValue(i, j, Double.MAX_VALUE);
                            matrixHolderClone.getFPSimilarityMatrix().setValue(i, j, Double.MIN_VALUE);
                        }
                    }
                }
                this.matrixHolder = matrixHolderClone;
            }
        }

        public Holder getMatrixHolder() { return matrixHolder; }
        public boolean isMatchFound() { return ruleMatched; }
        private void setRuleMatched(boolean ruleMatched) { this.ruleMatched = ruleMatched; }

        private void setRulesSmiles() throws CDKException {
            SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
            smartsWater = smilesParser.parseSmiles("O");
            smartsPhosphate = smilesParser.parseSmiles("OP(O)(O)=O");
            smartsDoublePhosphate = smilesParser.parseSmiles("OP(O)(=O)OP(O)(O)=O");
            smartsSulphate = smilesParser.parseSmiles("O=S(=O)(O)O");
            smartsL_Glutamate = smilesParser.parseSmiles("N[C@@H](CCC(O)=O)C(O)=O");
            smartsL_Glutamine = smilesParser.parseSmiles("N[C@@H](CCC(N)=O)C(O)=O");
            smartsL_Glutamine_clipped = smilesParser.parseSmiles("O=[C]N.O=C(O)C(N)C[CH2]");
            smartsL_Glutamate_clipped = smilesParser.parseSmiles("O=[C]O.O=C(O)C(N)C[CH2]");
            smartsTwoOxoglutarate = smilesParser.parseSmiles("OC(=O)CCC(=O)C(O)=O");
            smartsD_Glutamate = smilesParser.parseSmiles("N[C@H](CCC(O)=O)C(O)=O");
            smartsAcetate = smilesParser.parseSmiles("CC(O)=O");
            smartsATP = smilesParser.parseSmiles("NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O");
            smartsADP = smilesParser.parseSmiles("NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O");
            smartsCoA = smilesParser.parseSmiles("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N)[C@@H](O)C(=O)NCCC(=O)NCCS");
            smartsAcetyl_CoA = smilesParser.parseSmiles("CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N");
            smartsC00003 = smilesParser.parseSmiles("NC(=O)C1=CC=C[N+](=C1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C(N)N=CN=C23)[C@@H](O)[C@H]1O");
            smartsC00006 = smilesParser.parseSmiles("NC(=O)C1=C[N+](=CC=C1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](OP(O)(O)=O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O");
            smartsC00004 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O");
            smartsC00005 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](OP(O)(O)=O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O");
            smartsAlanine = smilesParser.parseSmiles("[CH3][C](N)C(O)=O");
            smartsPyruvate = smilesParser.parseSmiles("[CH3][C](=O)C(O)=O");
            smartsNRule = smilesParser.parseSmiles("CC(C)[C@H](N)C(O)=O");
            smartsCRule = smilesParser.parseSmiles("CC(C)C(=O)C(O)=O");
            smartsC04666Rule = smilesParser.parseSmiles("O=P(O)(O)O[CH2].[CH]O.O[CH]C=1N=CNC1");
            smartsC04916Rule = smilesParser.parseSmiles("O=C(N)C=1N=CN(C1N=CNCC(=O)[CH]O)C(O[CH])C(O)[CH]O.O=P(O)(O)O[CH2].O=P(O)(O)O[CH2].[CH]O");
        }

        private IAtomContainer getSmartsWater() { return smartsWater; }
        private IAtomContainer getSmartsPhosphate() { return smartsPhosphate; }
        private IAtomContainer getSmartsDoublePhosphate() { return smartsDoublePhosphate; }
        private IAtomContainer getSmartsGlutamate() { return smartsL_Glutamate; }
        private IAtomContainer getSmartsGlutamine() { return smartsL_Glutamine; }
        private IAtomContainer getSmartsGlutamineClipped() { return smartsL_Glutamine_clipped; }
        private IAtomContainer getSmartsGlutamateClipped() { return smartsL_Glutamate_clipped; }
        private IAtomContainer getSmartsTwoOxoglutarate() { return smartsTwoOxoglutarate; }
        private IAtomContainer getSmartsD_Glutamate() { return smartsD_Glutamate; }
        private IAtomContainer getSmartsAcetate() { return smartsAcetate; }
        private IAtomContainer getSmartsSulphate() { return smartsSulphate; }
        public IAtomContainer getSmartsATP() { return smartsATP; }
        public IAtomContainer getSmartsADP() { return smartsADP; }
        public IAtomContainer getSmartsCoA() { return smartsCoA; }
        public IAtomContainer getSmartsAcetyl_CoA() { return smartsAcetyl_CoA; }
        public IAtomContainer getSmartsC00003() { return smartsC00003; }
        public IAtomContainer getSmartsC00006() { return smartsC00006; }
        public IAtomContainer getSmartsC00004() { return smartsC00004; }
        public IAtomContainer getSmartsC00005() { return smartsC00005; }
        public IAtomContainer getSmartsPyruvate() { return smartsPyruvate; }
        public IAtomContainer getSmartsAlanine() { return smartsAlanine; }
        public IAtomContainer getSmartsNRule() { return smartsNRule; }
        public IAtomContainer getSmartsCRule() { return smartsCRule; }
        public IAtomContainer getSmartsC04666Rule() { return smartsC04666Rule; }
        public IAtomContainer getSmartsC04916Rule() { return smartsC04916Rule; }

        private boolean phosphate_cleaved(Collection<IAtomContainer> molsE, Collection<IAtomContainer> molsP) {
            int countphosE = 0;
            int countphosP = 0;
            for (IAtomContainer ac : molsE) {
                try {
                    if (isMatch(getSmartsPhosphate(), ac, false)) {
                        countphosE += ac.getAtomCount();
                    }
                } catch (CDKException ex) {
                    Logger.getLogger(RuleBasedMappingHandler.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            for (IAtomContainer ac : molsP) {
                try {
                    if (isMatch(getSmartsPhosphate(), ac, false)) {
                        countphosP += ac.getAtomCount();
                    }
                } catch (CDKException ex) {
                    Logger.getLogger(RuleBasedMappingHandler.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            return countphosE != countphosP;
        }
    }
}
