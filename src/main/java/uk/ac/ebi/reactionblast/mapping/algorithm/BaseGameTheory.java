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
import static java.lang.String.valueOf;
import static java.lang.System.out;
import java.util.BitSet;
import java.util.Calendar;
import static java.util.Calendar.DATE;
import static java.util.Calendar.HOUR;
import static java.util.Calendar.MILLISECOND;
import static java.util.Calendar.MINUTE;
import static java.util.Calendar.MONTH;
import static java.util.Calendar.YEAR;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Set;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import static uk.ac.ebi.reactionblast.mapping.graph.GraphMatcher.matcher;
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
    private final static ILoggingTool LOGGER
            = createLoggingTool(BaseGameTheory.class);
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

    /**
     *
     * @return @throws IOException
     */
    @Override
    public synchronized String getSuffix() throws IOException {
        Calendar cal = new GregorianCalendar();
        int ms = cal.get(YEAR);

        String suffix = valueOf(ms);
        ms = cal.get(MONTH);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(DATE);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(HOUR);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(MINUTE);
        suffix = suffix.concat(valueOf(ms));
        ms = cal.get(MILLISECOND);
        suffix = suffix.concat(valueOf(ms));
        //System.LOGGER.debug("Suffix: " + suffix);
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
                out.println("**********Updated Matrix And Calculate Similarity**************");
            }
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            Collection<MCSSolution> mcsSolutions = null;
            try {
                mcsSolutions = matcher(mh);
            } catch (InterruptedException e) {
                LOGGER.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
            }
            for (int substrateIndex = 0; substrateIndex < reactionStructureInformation.getEductCount(); substrateIndex++) {
                for (int productIndex = 0; productIndex < reactionStructureInformation.getProductCount(); productIndex++) {
                    try {
                        IAtomContainer educt = reactionStructureInformation.getEduct(substrateIndex);
                        IAtomContainer product = reactionStructureInformation.getProduct(productIndex);
                        if (DEBUG) {
                            System.out.println("===================Mapped==================");
                            System.out.println(" educt " + educt.getID() + ", product " + product.getID());
                            System.out.println(" educt Index " + substrateIndex + ", product Index " + productIndex);
                            System.out.println(" educt count " + educt.getAtomCount() + ", product count " + product.getAtomCount());

                            System.out.println("mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex) "
                                    + mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex));
                            System.out.println("reactionStructureInformation.isEductModified(substrateIndex) "
                                    + reactionStructureInformation.isEductModified(substrateIndex));
                            System.out.println("reactionStructureInformation.isProductModified(productIndex) "
                                    + reactionStructureInformation.isProductModified(productIndex));
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
                            mh.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFragmentMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getEnergyMatrix().setValue(substrateIndex, productIndex, 0.0);
                            mh.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, 0.0);
                        }
                    } catch (IOException | CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                }
            }
        } catch (Exception e) {
            LOGGER.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
        }

        try {
            /*
             * reset mapping flags to FALSE to allow remapping if needed
             *
             */
            resetFLAGS(mh);
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
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
            if (DEBUG) {
                System.out.println("**********Updated Matrix And Calculate Similarity**************");
            }
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
                        mh.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, 0.0);
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
            LOGGER.error("Error in matching molecules, check Graph Matcher module! ", e.toString());
        }
    }

    private synchronized void refillMatrixWithNewData(
            Holder holder,
            int substrateIndex,
            int productIndex,
            Collection<MCSSolution> mcsSolutions) {
        if (DEBUG) {
            out.println("**********Generate NEW MCS And Calculate Similarity**************");
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
            double carbonCount = 0.0;

            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);

            if (DEBUG) {
                out.println("--Get matches--");
                out.print("Q " + educt.getID() + ": " + educt.getAtomCount());
                out.print(", P " + product.getID() + ": " + product.getAtomCount());
                out.println(", Matches " + mcsSolutions.size());
            }

            MCSSolution atomatomMapping = getMappings(substrateIndex, productIndex, educt, product, mcsSolutions);

            carbonCount = atomatomMapping.getAtomAtomMapping().getMappingsByAtoms().keySet().stream().filter((atom)
                    -> (atom.getSymbol().equalsIgnoreCase("C"))).map((IAtom _item) -> 1.0).reduce(carbonCount, (accumulator, _item)
                    -> accumulator + 1);

            if (atomatomMapping == null) {
                throw new CDKException("atom-atom mapping is null");
            }

            if (DEBUG) {
                out.println("set matching atoms");
                out.print("Q " + educt.getID() + ": " + educt.getAtomCount());
                out.print(", P " + product.getID() + ": " + product.getAtomCount());
                out.println(", Matches " + mcsSolutions.size());
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
                out.println(substrateIndex + " KEY " + productIndex + ", MCS Mapping Size " + mappingSize);
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
                    fpSim = getTanimotoSimilarity(a, b);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }

            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, carbonCount);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);

        } catch (IOException | CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    private synchronized MCSSolution getMappings(
            int queryPosition,
            int targetPosition,
            IAtomContainer educt,
            IAtomContainer product,
            Collection<MCSSolution> mcsSolutions) throws CDKException {

        if (DEBUG) {
            System.out.println("====++mcsSolutions++=== " + mcsSolutions.size());
        }
        if (mcsSolutions.isEmpty()) {
            return quickMapping(educt, product, queryPosition, targetPosition);
        }
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
                            SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
                            System.out.println(educt.getID() + " ED: " + smilesGenerator.create(educt));
                            System.out.println(product.getID() + " PD: " + smilesGenerator.create(product));
                        }
                        return quickMapping(educt, product, queryPosition, targetPosition);
                    }
                }

                if (DEBUG) {
                    System.out.println("--------======NEW MCS======------" + solution.getAtomAtomMapping().getCount());
                    System.out.println(" educt " + educt.getID() + ", product " + product.getID());
                    System.out.println(" educt Index " + queryPosition + ", product Index " + targetPosition);
                    System.out.println(" educt count " + educt.getAtomCount() + ", product count " + product.getAtomCount());
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
        if (DEBUG) {
            System.out.println("====Quick Mapping====");
        }
        try {
            //Isomorphism mcsThread = new Isomorphism(educt, product, Algorithm.VFLibMCS, false, false, false);
            Isomorphism isomorphism;
            isomorphism = new Isomorphism(educt, product, Algorithm.DEFAULT, false, false, false);
            isomorphism.setChemFilters(true, true, true);
            MCSSolution mcs = new MCSSolution(queryPosition, targetPosition, educt, product, isomorphism.getFirstAtomMapping());
            mcs.setEnergy(isomorphism.getEnergyScore(0));
            mcs.setFragmentSize(isomorphism.getFragmentSize(0));
            mcs.setStereoScore(isomorphism.getStereoScore(0));
            return mcs;
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
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
        if (DEBUG) {
            System.out.println("**********REFILL MCS And Calculate Similarity**************");
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
            double carbonCount = 0.0;

            IAtomContainer educt = reactionContainer.getEduct(substrateIndex);
            IAtomContainer product = reactionContainer.getProduct(productIndex);

            if (initMcsAtom.containsKey(substrateIndex, productIndex)) {
                AtomAtomMapping bestAtomAtomMapping = initMcsAtom.getAtomMatch(substrateIndex, productIndex);

                if (bestAtomAtomMapping == null) {
                    throw new CDKException("atom-atom mapping is null");
                }

                carbonCount = bestAtomAtomMapping.getMappingsByAtoms().keySet().stream().filter((atom) -> (atom.getSymbol().equalsIgnoreCase("C"))).map((_item) -> 1.0).reduce(carbonCount, (accumulator, _item) -> accumulator + 1);

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
                        fpSim = getTanimotoSimilarity(a, b);
                    } catch (Exception ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                }
            }

            holder.getCliqueMatrix().setValue(substrateIndex, productIndex, mappingSize);
            holder.getGraphSimilarityMatrix().setValue(substrateIndex, productIndex, graphSimilarity);
            holder.getStereoMatrix().setValue(substrateIndex, productIndex, stereoVal);
            holder.getCarbonOverlapMatrix().setValue(substrateIndex, productIndex, carbonCount);
            holder.getFragmentMatrix().setValue(substrateIndex, productIndex, fragmentVal);
            holder.getEnergyMatrix().setValue(substrateIndex, productIndex, energyVal);
            holder.getFPSimilarityMatrix().setValue(substrateIndex, productIndex, fpSim);
        } catch (CDKException ex) {
            LOGGER.debug(SEVERE, null, ex);
        } catch (IOException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }
}
