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
import java.util.List;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.interfaces.IResult;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class IsomorphismMax extends Selector implements IResult {

    private final static ILoggingTool LOGGER
            = createLoggingTool(IsomorphismMax.class);
    private static final long serialVersionUID = 0x192aa60a59L;
    private final Holder mHolder;
    private final Holder updatedHolder;
    private boolean SubAndCompleteFlag;

    /**
     *
     * @param mHolder
     * @param EdMapOrignal
     * @param PdMapOrignal
     * @throws IOException
     * @throws Exception
     */
    public IsomorphismMax(Holder mHolder,
            List<String> EdMapOrignal,
            List<String> PdMapOrignal) throws IOException, Exception {

//        System.out.println("\n********* INSIDE IsomorphismMax**********");
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

    private synchronized void PhaseOneMatcher() throws IOException, CDKException {

//        System.out.println("\nPhase 1\n");
        boolean maxValueI;
        boolean maxValueJ;

        for (int i = 0; i < rowSize; i++) {
            IAtomContainer ac1 = this.mHolder.getReactionContainer().getEduct(i);
            for (int j = 0; j < colSize; j++) {
                IAtomContainer ac2 = this.mHolder.getReactionContainer().getProduct(j);
                double cliqueVal = mHolder.getCliqueMatrix().getValue(i, j);
                double simVal = mHolder.getGraphSimilarityMatrix().getValue(i, j);
                //matrix.
                if (cliqueVal > 0.0) {
                    maxValueI = isMajorSubgraphRow(mHolder, i, j);//mHolder.getCliqueMatrix().is_element_max_in_row(i, j); //, rowSize);
                    maxValueJ = isMajorSubgraphColumn(mHolder, i, j);//mHolder.getCliqueMatrix().is_element_max_in_column(i, j); //, colSize);
                    if (maxValueI || maxValueJ) {
                        if (simVal > 0 && ac1 != null && ac2 != null) {
                            //Here we check for bond count and in min we look for atom count
                            int eMolSize = ac1.getAtomCount();
                            int pMolSize = ac2.getAtomCount();

                            /*
                             * This check will test for single atoms like H or O to be mapped first
                             */
                            if (eMolSize == 1 && pMolSize == 1
                                    && (ac1.atoms().iterator().next().getSymbol()
                                            .equals(ac2.atoms().iterator().next().getSymbol()))) {
                                this.flagMatrix[i][j] = true;
                            } /*
                             * This check will skip single atoms like H or O to be mapped with larger graph as that
                             * mightr be a flase postivie match
                             */ else if (eMolSize > 1 && pMolSize > 1
                                    && ac1.getAtomCount() <= ac2.getAtomCount()
                                    && isMCSSubgraph(ac1, mHolder.getCliqueMatrix().getValue(i, j))) {
//                                    && isSubGraph(ac1, ac2)) {
//                            System.out.println("sim " + simVal);
                                this.flagMatrix[i][j] = true;
                            } else if (eMolSize > 1 && pMolSize > 1
                                    && ac2.getAtomCount() <= ac1.getAtomCount()
                                    && isMCSSubgraph(ac2, mHolder.getCliqueMatrix().getValue(i, j))) {
//                                    && isSubGraph(ac1, ac2)) {
//                            System.out.println("sim " + simVal);
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

    private synchronized boolean PhaseTwoMatcher() throws IOException, CDKException {

//        System.out.println("**********Substructural Modification Part2: Matrix**************");
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

    private synchronized boolean isMCSSubgraph(IAtomContainer educt, double mcsSize) throws CDKException {
        return educt.getAtomCount() == mcsSize;
    }

    private synchronized boolean isSubGraph(IAtomContainer educt, IAtomContainer product) {
        try {
            if (educt.getAtomCount() <= product.getAtomCount()) {
                Substructure s = new Substructure(educt, product, false, true, true, false);
                return s.isSubgraph();
            } else if (educt.getAtomCount() > product.getAtomCount()) {
                Substructure s = new Substructure(product, educt, false, true, true, false);
                return s.isSubgraph();
            }
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return false;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized boolean isSubAndCompleteMatchFlag() {
        return SubAndCompleteFlag;
    }

    /**
     * @return the updatedHolder
     */
    @Override
    public synchronized Holder getUpdatedHolder() {
        return updatedHolder;
    }
}
