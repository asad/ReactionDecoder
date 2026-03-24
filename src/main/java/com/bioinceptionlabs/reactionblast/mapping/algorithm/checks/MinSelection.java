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
package com.bioinceptionlabs.reactionblast.mapping.algorithm.checks;

import java.io.IOException;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;
import com.bioinceptionlabs.reactionblast.mapping.container.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.interfaces.IResult;

/**
 *
 * @author Syed Asad Rahman, BioInception
 * @contact asad.rahman@bioinceptionlabs.com
 */
public class MinSelection extends Selector implements IResult {

    private final static ILoggingTool LOGGER
            = createLoggingTool(MinSelection.class);
    private static final long serialVersionUID = 1908987778L;
    private final Holder mHolder;
    private final Holder updatedHolder;
    private boolean isSubstructure;

    /**
     *
     * @param mHolder
     * @param EdMapOrignal
     * @param PdMapOrignal
     * @throws IOException
     * @throws CloneNotSupportedException
     */
    public MinSelection(Holder mHolder,
            List<String> EdMapOrignal,
            List<String> PdMapOrignal) throws IOException, CloneNotSupportedException {

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
                //matrix.
                if (cliqueVal > 0.0) {
                    double simVal = mHolder.getGraphSimilarityMatrix().getValue(i, j);
                    if (simVal > 0 && ac1 != null && ac2 != null) {
                        int eMolSize = ac1.getAtomCount();
                        int pMolSize = ac2.getAtomCount();

                        /*
                         * This check will test for single atoms like H or O to be mapped first
                         */
                        if (eMolSize == 1 && pMolSize == 1
                                && (ac1.atoms().iterator().next().getSymbol().equals(ac2.atoms().iterator().next().getSymbol()))) {
                            this.flagMatrix[i][j] = true;
                        } /*
                         * This check will skip single atoms like H or O to be mapped with larger graph as that mightr
                         * be a flase postive match
                         */ else if (eMolSize > 1 && pMolSize > 1 && ac1.getAtomCount() <= ac2.getAtomCount()
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
//                    TO DO ASAD check if one has to choose substructure with lowest energy
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

    /**
     *
     * @return
     */
    @Override
    public boolean isSubAndCompleteMatchFlag() {
        return isSubstructure;
    }

    private boolean isMCSSubgraph(IAtomContainer educt, double mcsSize) throws CDKException {
        return educt.getAtomCount() == (int) mcsSize;
    }

    /**
     * @return the updatedHolder
     */
    @Override
    public Holder getUpdatedHolder() {
        return updatedHolder;
    }

    /**
     * Prints Similarity Matrix
     *
     * @param mh
     */
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

    /**
     * Prints Similarity Matrix
     *
     * @param mh
     */
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

    /**
     * Prints FLAG Matrix
     *
     * @param mh
     */
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
