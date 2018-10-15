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
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import static java.lang.System.out;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.interfaces.IResult;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class IsomorphismMin extends Selector implements IResult {

    private final static ILoggingTool LOGGER
            = createLoggingTool(IsomorphismMin.class);
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
    public IsomorphismMin(Holder mHolder,
            List<String> EdMapOrignal,
            List<String> PdMapOrignal) throws IOException, CloneNotSupportedException {

//        System.out.println("Check Subgraph IsomorphismMin ");
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

    private synchronized void PhaseOneMatcher() throws IOException, CDKException {
        this.flagMatrix = new boolean[rowSize][colSize];

        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                this.flagMatrix[i][j] = false;
            }
        }

//        System.out.println("PhaseOneMatcher\n");
//        printSimMatrix(mHolder);
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
//                            System.out.println("Subgraph sim " + simVal);
                            this.flagMatrix[i][j] = true;
                        } else if (eMolSize > 1 && pMolSize > 1 && ac2.getAtomCount() <= ac1.getAtomCount()
                                && isMCSSubgraph(ac2, mHolder.getCliqueMatrix().getValue(i, j))) {
//                            System.out.println("Subgraph sim " + simVal);
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

    private synchronized boolean PhaseTwoMatcher() throws IOException, CDKException {

//        System.out.println("**********Substructural Modification Part1: Matrix**************");
//        printFLAGMatrix(updatedHolder);
//        printSimMatrix(updatedHolder);
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
//                            System.out.println("check isMinorSubgraph " + ", " + ac1.getID() + " : " + ac2.getID());
                            boolean isMinorSubgraphColumn = isMinorSubgraphColumn(updatedHolder, i, j);
                            boolean isMinorSubgraphRow = isMinorSubgraphRow(updatedHolder, i, j);

//                            System.out.println("check isMinorSubgraphColumn " + ", "
//                                    + isMinorSubgraphColumn + " : isMinorSubgraphRow "
//                                    + isMinorSubgraphRow);
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
//        System.out.println("**********Substructural Modification Part2: Matrix**************");
//        printFLAGMatrix(updatedHolder);
//        printSimMatrix(updatedHolder);
        return flag;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized boolean isSubAndCompleteMatchFlag() {
        return isSubstructure;
    }

    private synchronized boolean isMCSSubgraph(IAtomContainer educt, double mcsSize) throws CDKException {
        return new Double(educt.getAtomCount()) == mcsSize;
    }

    /**
     * @return the updatedHolder
     */
    @Override
    public synchronized Holder getUpdatedHolder() {
        return updatedHolder;
    }

    /**
     * Prints Similarity Matrix
     *
     * @param mh
     */
    static void printSimMatrix(Holder mh) {
        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
        out.println();
        out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            out.println("Similarity Matrix");
            out.print("\t\t");
            for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                out.print(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            out.println();
            double val = 0;
            for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                out.print(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    val = mh.getGraphSimilarityMatrix().getValue(i, j);
                    result = format.format(val);
                    out.print("   " + result);
                }
                out.println();
            }
        } catch (IOException | CDKException e) {
            LOGGER.debug(" Parser Error: ");
        }
        out.println();

    }

    /**
     * Prints Similarity Matrix
     *
     * @param mh
     */
    static void printCliqueMatrix(Holder mh) {
        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
        out.println();
        out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            out.println("Clique Matrix");
            out.print("\t\t");
            for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                out.print(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            out.println();
            double val = 0;
            for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                out.print(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    val = mh.getCliqueMatrix().getValue(i, j);
                    result = format.format(val);
                    out.print("   " + result);
                }
                out.println();
            }
        } catch (IOException | CDKException e) {
            LOGGER.debug(" Parser Error: ");
        }
        out.println();

    }

    /**
     * Prints FLAG Matrix
     *
     * @param mh
     */
    static void printFLAGMatrix(Holder mh, boolean[][] flagMatrix) {
        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
        out.println();
        out.println("********* MATRIX **********");
        try {
            String result;
            out.println("Flag Matrix");
            out.print("\t\t");
            for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                out.print(" " + reactionStructureInformationContainer.getProduct(j).getID() + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            out.println();
            boolean val;
            for (int i = 0; i < reactionStructureInformationContainer.getEductCount(); i++) {
                out.print(" " + reactionStructureInformationContainer.getEduct(i).getID() + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < reactionStructureInformationContainer.getProductCount(); j++) {
                    out.print("   " + flagMatrix[i][j]);
                }
                out.println();
            }
        } catch (IOException | CDKException e) {
            LOGGER.error("Parser Error", e);
        }
        out.println();

    }
}
