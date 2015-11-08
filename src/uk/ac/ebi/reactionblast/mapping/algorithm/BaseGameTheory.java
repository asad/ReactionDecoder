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

package uk.ac.ebi.reactionblast.mapping.algorithm;

import java.io.IOException;
import java.io.Serializable;
import java.util.BitSet;
import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.graph.GraphMatcher;
import uk.ac.ebi.reactionblast.mapping.graph.MCSSolution;
import uk.ac.ebi.reactionblast.mapping.helper.Debugger;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.mapping.interfaces.IGameTheory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class BaseGameTheory extends Debugger implements IGameTheory, Serializable {

    private final static boolean DEBUG = false;
    private final static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(BaseGameTheory.class);
    private static final long serialVersionUID = 1698688633678282L;

    /**
     * Checks if a PseudoAtom is present
     *
     * @param atomContainer
     * @return
     */
    protected static synchronized boolean isPseudoAtoms(IAtomContainer atomContainer) {
        for (IAtom atoms : atomContainer.atoms()) {
            if (atoms instanceof IPseudoAtom || atoms instanceof PseudoAtom) {
                return true;
            }
        }
        return false;
    }

    @Override
    public synchronized String getSuffix() throws IOException {
        Calendar cal = new GregorianCalendar();
        int ms = cal.get(Calendar.YEAR);

        String suffix = String.valueOf(ms);
        ms = cal.get(Calendar.MONTH);
        suffix = suffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.DATE);
        suffix = suffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.HOUR);
        suffix = suffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.MINUTE);
        suffix = suffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.MILLISECOND);
        suffix = suffix.concat(String.valueOf(ms));
        //System.err.println("Suffix: " + suffix);
        return suffix;

    }

    /**
     *
     * @param mh matrix holder
     * @param removeHydrogen
     * @throws InterruptedException
     */
    @Override
    public synchronized void UpdateMatrix(Holder mh, boolean removeHydrogen) throws InterruptedException {
        try {
            if (DEBUG) {
                System.out.println("**********Updated Matrix And Calculate Similarity**************");
            }
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Collection<MCSSolution> mcsSolutions = null;
            try {
                mcsSolutions = GraphMatcher.matcher(mh);
            } catch (Exception e) {
                logger.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
            }
            for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
                for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
                    try {
                        IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                        IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                        if (DEBUG) {
                            System.out.println("mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) "
                                    + mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex));
                        }
                        if ((educt != null && product != null)
                                && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                                && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                                || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
                            if (reactionStructureInformation.isEductModified(substrateIndex)
                                    || reactionStructureInformation.isProductModified(productIndex)) {
                                refillMatrixWithNewData(mh, substrateIndex, productIndex, mcsSolutions);
                            } else {
                                refillMatrixWithOldData(mh, substrateIndex, productIndex);
                            }
                        } else {
                            mh.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getStereoMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getCliqueMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                        }
                    } catch (IOException | CDKException ex) {
                        logger.error(Level.SEVERE, null, ex);
                    }
                }
            }
        } catch (Exception e) {
            logger.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
        }

        try {
            /*
             * reset mapping flags to FALSE to allow remapping if needed
             *
             */
            resetFLAGS(mh);
        } catch (Exception ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @param mcsSolutions
     * @param mh matrix holder
     * @param removeHydrogen
     * @throws Exception
     */
    @Override
    public synchronized void UpdateMatrix(Collection<MCSSolution> mcsSolutions, Holder mh, boolean removeHydrogen) throws Exception {
        try {
//        System.out.println("**********Updated Matrix And Calculate Similarity**************");
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();

            for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
                for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {

                    IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                    IAtomContainer product = reactionStructureInformation.getProduct(productIndex);

                    if ((educt != null && product != null)
                            && (reactionStructureInformation.getEduct(substrateIndex).getAtomCount() > 0
                            && reactionStructureInformation.getProduct(productIndex).getAtomCount() > 0)
                            || mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) == -1) {
                        if (reactionStructureInformation.isEductModified(substrateIndex)
                                || reactionStructureInformation.isProductModified(productIndex)) {
                            refillMatrixWithNewData(mh, substrateIndex, productIndex, mcsSolutions);
                        } else {
                            refillMatrixWithOldData(mh, substrateIndex, productIndex);
                        }
                    } else {
                        mh.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getStereoMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getCliqueMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
                        mh.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                    }
                }
            }


            /*
             * reset mapping flags to FALSE to allow remapping if needed
             *
             */
            resetFLAGS(mh);
        } catch (Exception e) {
            logger.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
        }
    }

    private synchronized void refillMatrixWithNewData(
            Holder holder,
            int substrateIndex,
            int productIndex,
            Collection<MCSSolution> mcsSolutions) {
        if (DEBUG) {
            System.out.println("**********Generate MCS And Calculate Similarity**************");
        }
        try {
            ReactionContainer reactionContainer = holder.getReactionContainer();
            BestMatch initMcsAtom = holder.getBestMatchContainer();

            double stereoVal = 0.0;
            int fragmentVal = 0;
            double energyVal = 0.0;
            double graphSimilarity = 0.0;
            double mappingSize = 0.0;
            double fpSim = 0.0;

            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);

            if (DEBUG) {
                System.out.println("Get matches");
                System.out.print("Q " + educt.getID() + ": " + educt.getAtomCount());
                System.out.print(", P " + product.getID() + ": " + product.getAtomCount());
                System.out.println(", Matches " + mcsSolutions.size());
            }

            MCSSolution atomatomMapping = getMappings(substrateIndex, productIndex, educt, product, mcsSolutions);

            if (atomatomMapping == null) {
                throw new CDKException("atom-atom mapping is null");
            }

            if (DEBUG) {
                System.out.println("set matching atoms");
                System.out.print("Q " + educt.getID() + ": " + educt.getAtomCount());
                System.out.print(", P " + product.getID() + ": " + product.getAtomCount());
                System.out.println(", Matches " + mcsSolutions.size());
            }
            if (atomatomMapping.getStereoScore() != null) {
                stereoVal = atomatomMapping.getStereoScore();
            }
            if (atomatomMapping.getFragmentSize() != null) {
                fragmentVal = atomatomMapping.getFragmentSize();
            }
            if (atomatomMapping.getEnergy() != null) {
                energyVal = atomatomMapping.getEnergy();
            }

            AtomAtomMapping fam = atomatomMapping.getAtomAtomMapping();
            initMcsAtom.putBestMapping(substrateIndex, productIndex, fam);

            double ACount = educt.getAtomCount();
            double BCount = product.getAtomCount();

            mappingSize = atomatomMapping.getAtomAtomMapping().getCount();
            if (DEBUG) {
                System.out.println(substrateIndex + " KEY " + productIndex + ", MCS Mapping Size " + mappingSize);
            }
            graphSimilarity = mappingSize / (ACount + BCount - mappingSize);

            initMcsAtom.setTotalFragmentCount(substrateIndex, productIndex, fragmentVal);
            initMcsAtom.setBondEnergy(substrateIndex, productIndex, energyVal);
            initMcsAtom.setStereoScore(substrateIndex, productIndex, stereoVal);
            initMcsAtom.setGraphSimilarity(substrateIndex, productIndex, graphSimilarity);

            BitSet a = reactionContainer.getFingerPrintofEduct(substrateIndex);
            BitSet b = reactionContainer.getFingerPrintofProduct(productIndex);
            if (a != null && b != null) {
                try {
                    fpSim = Similarity.getTanimotoSimilarity(a, b);
                } catch (Exception ex) {
                    Logger.getLogger(BaseGameTheory.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);

        } catch (IOException | CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }

    private synchronized MCSSolution getMappings(
            int queryPosition,
            int targetPosition,
            IAtomContainer educt,
            IAtomContainer product,
            Collection<MCSSolution> mcsSolutions) throws CDKException {

        for (MCSSolution solution : mcsSolutions) {
            if (solution.getQueryPosition() == queryPosition
                    && solution.getTargetPosition() == targetPosition) {
                if (solution.getAtomAtomMapping().isEmpty()) {
                    Set<String> atomMaps = new HashSet<>();

                    for (IAtom a : educt.atoms()) {
                        atomMaps.add(a.getSymbol());
                    }
                    boolean mappingPossible = false;
                    for (IAtom a : product.atoms()) {
                        if (atomMaps.contains(a.getSymbol())) {
                            mappingPossible = true;
                        }
                    }
                    atomMaps.clear();
                    if (mappingPossible) {
                        if (DEBUG) {
                            System.out.println("Expected Mapping");
                            System.out.println(educt.getID() + " ED: " + SmilesGenerator.unique().aromatic().create(educt));
                            System.out.println(product.getID() + " PD: " + SmilesGenerator.unique().aromatic().create(product));
                        }
                        return quickMapping(educt, product, queryPosition, targetPosition);
                    }
                }
                return solution;
            }
        }

        return null;
    }

    private MCSSolution quickMapping(IAtomContainer educt, IAtomContainer product, int queryPosition, int targetPosition) {
        /*
         * This function is called as a backup emergency step to avoid null if matching is possible
         */
        Isomorphism mcsThread = new Isomorphism(educt, product, Algorithm.DEFAULT, false, false, false);
        mcsThread.setChemFilters(true, true, true);
        try {
            MCSSolution mcs = new MCSSolution(queryPosition, targetPosition, educt, product, mcsThread.getFirstAtomMapping());
            mcs.setEnergy(mcsThread.getEnergyScore(0));
            mcs.setFragmentSize(mcsThread.getFragmentSize(0));
            mcs.setStereoScore(mcsThread.getStereoScore(0));
            return mcs;
        } catch (Exception ex) {
            logger.error(Level.SEVERE, null, ex);
        }
        return null;
    }

    private void resetFLAGS(Holder mh) throws Exception {
        ReactionContainer reactionStructureInformation = mh.getReactionContainer();
        /*
         * Reset all the flags
         */
        for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
            for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
                reactionStructureInformation.setEductModified(substrateIndex, false);
                reactionStructureInformation.setProductModified(productIndex, false);
            }
        }
    }

    private void refillMatrixWithOldData(Holder holder, int substrateIndex, int productIndex) {
        //        System.out.println("**********Generate MCS And Calculate Similarity**************");
        try {
            ReactionContainer reactionContainer = holder.getReactionContainer();
            BestMatch initMcsAtom = holder.getBestMatchContainer();

            double stereoVal = 0.0;
            int fragmentVal = 0;
            double energyVal = 0.0;
            double graphSimilarity = 0.0;
            double mappingSize = 0.0;
            double fpSim = 0.0;

            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);

            if (initMcsAtom.containsKey(substrateIndex, productIndex)) {
                AtomAtomMapping bestAtomAtomMapping = initMcsAtom.getAtomMatch(substrateIndex, productIndex);

                if (bestAtomAtomMapping == null) {
                    throw new CDKException("atom-atom mapping is null");
                }
                stereoVal = initMcsAtom.getStereoScore(substrateIndex, productIndex);
                fragmentVal = initMcsAtom.getTotalFragmentCount(substrateIndex, productIndex);
                energyVal = initMcsAtom.getBondEnergy(substrateIndex, productIndex);

                double ACount = educt.getAtomCount();
                double BCount = product.getAtomCount();

                mappingSize = bestAtomAtomMapping.getCount();
//                System.out.println("KEY " + key + ", MCS Mapping Size " + mappingSize);
                graphSimilarity = mappingSize / (ACount + BCount - mappingSize);

                initMcsAtom.setTotalFragmentCount(substrateIndex, productIndex, fragmentVal);
                initMcsAtom.setBondEnergy(substrateIndex, productIndex, energyVal);
                initMcsAtom.setStereoScore(substrateIndex, productIndex, stereoVal);
                initMcsAtom.setGraphSimilarity(substrateIndex, productIndex, graphSimilarity);

                BitSet a = reactionContainer.getFingerPrintofEduct(substrateIndex);
                BitSet b = reactionContainer.getFingerPrintofProduct(productIndex);
                if (a != null && b != null) {
                    try {
                        fpSim = Similarity.getTanimotoSimilarity(a, b);
                    } catch (Exception ex) {
                        logger.error(Level.SEVERE, null, ex);
                    }
                }
            }

            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);
        } catch (CDKException ex) {
            logger.debug(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }
}
