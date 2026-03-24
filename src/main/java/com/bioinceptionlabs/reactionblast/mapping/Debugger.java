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
package com.bioinceptionlabs.reactionblast.mapping;

import static java.io.File.separator;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.IAtomMapping;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.ChooseWinner;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.tools.BasicDebugger;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import com.bioinceptionlabs.reactionblast.tools.EBIMatrix;
import com.bioinceptionlabs.reactionblast.tools.ImageGenerator;
import static java.lang.System.getProperty;
import static java.text.NumberFormat.getInstance;

import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public abstract class Debugger extends BasicDebugger {

    private final static ILoggingTool LOGGER
            = createLoggingTool(Debugger.class);

    /**
     * Prints reactant and product atom container in the matrix
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printMatrixAtomContainer(Holder mh, List<String> EdMap, List<String> PdMap) {
        try {
            ReactionContainer _rSTMap = mh.getReactionContainer();
            StringBuilder sb = new StringBuilder();
            sb.append("<--------Atom Size in the Container-------->").append(NEW_LINE);
            for (int i = 0; i < EdMap.size(); i++) {
                sb.append("Educt ").append(EdMap.get(i)).append(" : ").append(_rSTMap.getEduct(i).getAtomCount()).append(NEW_LINE);
                if (!_rSTMap.getEduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getEduct(i), true, false);
                    sb.append("SMILES: ").append(sm.getCanonicalSMILES()).append(NEW_LINE);
                }
                printAtoms(_rSTMap.getEduct(i));
            }
            sb.append(NEW_LINE);
            for (int i = 0; i < PdMap.size(); i++) {
                sb.append("Product ").append(PdMap.get(i)).append(" : ").append(_rSTMap.getProduct(i).getAtomCount()).append(NEW_LINE);
                if (!_rSTMap.getProduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getProduct(i), true, false);
                    sb.append("SMILES: ").append(sm.getCanonicalSMILES()).append(NEW_LINE);
                }
                printAtoms(_rSTMap.getProduct(i));
            }
            LOGGER.debug(sb.toString());
        } catch (IOException | CDKException | CloneNotSupportedException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     * Prints Clique Matrix
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printCliqueMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {

        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            sb.append("Clique Matrix").append(NEW_LINE);
            sb.append("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                sb.append(" ").append(PdMap.get(j)).append(":(").append(reactionStructureInformationContainer.getProduct(j).getAtomCount()).append(")");
            }
            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i)).append(":(").append(reactionStructureInformationContainer.getEduct(i).getAtomCount()).append(")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getCliqueMatrix().getValue(i, j);
                    result = format.format(val);
                    sb.append("   ").append(result);
                }
                sb.append(NEW_LINE);
            }
        } catch (IOException | CDKException e) {
            LOGGER.debug("Parser Error" + e);
        }
        LOGGER.debug(sb.toString());
    }

    /**
     * Prints Similarity Matrix
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printSimMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            sb.append("Similarity Matrix").append(NEW_LINE);
            sb.append("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                sb.append(" ").append(PdMap.get(j)).append(":(").append(reactionStructureInformationContainer.getProduct(j).getAtomCount()).append(")");
            }
            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i)).append(":(").append(reactionStructureInformationContainer.getEduct(i).getAtomCount()).append(")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getGraphSimilarityMatrix().getValue(i, j);
                    result = format.format(val);
                    sb.append("   ").append(result);
                }
                sb.append(NEW_LINE);
            }
        } catch (IOException | CDKException e) {
            LOGGER.debug("Parser Error" + e);
        }
        LOGGER.debug(sb.toString());

    }

    /**
     *
     * @param winner
     * @param EdMap
     * @param PdMap
     */
    protected void printFlagMatrix(ChooseWinner winner, List<String> EdMap, List<String> PdMap) {

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);
        boolean[][] FlagMatrix = winner.getFlagMatrix();
        sb.append("Flag Matrix").append(NEW_LINE);
        sb.append("\t\t");
        PdMap.forEach((PdMap1) -> {
            sb.append("  ").append(PdMap1).append(" ");
        });

        sb.append(NEW_LINE);
        for (int i = 0; i < EdMap.size(); i++) {
            sb.append(" ").append(EdMap.get(i));
            for (int j = 0; j < PdMap.size(); j++) {
                if (FlagMatrix[i][j]) {
                    sb.append("      ").append(1).append("  ");
                } else {
                    sb.append("      ").append(0).append("  ");
                }

            }
            sb.append(NEW_LINE);
        }
        LOGGER.debug(sb.toString());

    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printStereoMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix StereoMatrix = mh.getStereoMatrix();

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            sb.append("Stereo Matrix").append(NEW_LINE);
            sb.append("\t\t");

            PdMap.forEach((PdMap1) -> {
                sb.append(" ").append(PdMap1);
            });

            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = StereoMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    sb.append("   ").append(result);
                }

                sb.append(NEW_LINE);
            }

        } catch (Exception e) {
            LOGGER.debug("Parser Error" + e);
        }

        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printFragmentMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix fragmentMatrix = mh.getFragmentMatrix();

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            sb.append("Fragment Matrix").append(NEW_LINE);
            sb.append("\t\t");

            PdMap.forEach((PdMap1) -> {
                sb.append(" ").append(PdMap1);
            });

            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = fragmentMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    sb.append("   ").append(result);
                }

                sb.append(NEW_LINE);
            }

        } catch (Exception e) {
            LOGGER.debug("Parser Error" + e);
        }

        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printCarbonMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix carbonMatrix = mh.getCarbonOverlapMatrix();

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            sb.append("Fragment Matrix").append(NEW_LINE);
            sb.append("\t\t");

            PdMap.forEach((PdMap1) -> {
                sb.append(" ").append(PdMap1);
            });

            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = carbonMatrix.getValue(i, j);
                    result = format.format(val);
                    sb.append("   ").append(result);
                }

                sb.append(NEW_LINE);
            }

        } catch (Exception e) {
            LOGGER.debug("Parser Error" + e);
        }

        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printEnergyMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix energyMatrixProfile = mh.getEnergyMatrix();

        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("********* MATRIX **********").append(NEW_LINE);

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            sb.append("Energy Matrix").append(NEW_LINE);
            sb.append("\t\t");

            PdMap.forEach((PdMap1) -> {
                sb.append("\t").append(PdMap1);
            });

            sb.append(NEW_LINE);
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                sb.append("\t").append(EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = energyMatrixProfile.getValue(i, j);
                    result = format.format(val);
                    sb.append("\t").append(result);
                }

                sb.append(NEW_LINE);
            }

        } catch (Exception e) {
            LOGGER.debug("Parser Error" + e);
        }

        LOGGER.debug(sb.toString());
    }

    /**
     * Print Graph matching solutions
     *
     * @param comparison
     * @param mol1
     * @param mol2
     */
    protected void printGraphMatching(IAtomMapping comparison, IAtomContainer mol1, IAtomContainer mol2) {
        int count_final_sol = 0;
        StringBuilder sb = new StringBuilder();
        sb.append("Output of the final Mappings: ").append(NEW_LINE);
        sb.append("Mol1: ").append(mol1.getID()).append(NEW_LINE);
        sb.append("Mol2: ").append(mol2.getID()).append(NEW_LINE);
        try {
            if (comparison.getMappingCount() > 0) {

                for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                    int final_solution_size = final_solution.getCount();
                    sb.append("Final mapping Nr. ").append(++count_final_sol)
                            .append(" Size:").append(final_solution_size).append(NEW_LINE);

                    final int solIndex = count_final_sol;
                    final_solution.getMappingsByAtoms().entrySet().forEach((mapping) -> {
                        IAtom eAtom = mapping.getKey();
                        IAtom pAtom = mapping.getValue();

                        sb.append(mol1.indexOf(eAtom) + 1).append(" ").append(mol2.indexOf(pAtom) + 1).append(NEW_LINE);

                        sb.append(eAtom.getSymbol()).append(" ")
                                .append(pAtom.getSymbol()).append(NEW_LINE);
                    });
                    sb.append("").append(NEW_LINE);

                    sb.append("Stereo Match: ").append(comparison.getStereoScore(count_final_sol - 1)).append(NEW_LINE);
                    sb.append("Stereo different: ").append(comparison.isStereoMisMatch()).append(NEW_LINE);
                    sb.append("Fragment Size: ").append(comparison.getFragmentSize(count_final_sol - 1)).append(NEW_LINE);
                }

                sb.append("").append(NEW_LINE);
            }
        } catch (Exception ex) {
            LOGGER.debug("Parser Error" + ex);
        }
        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param outPutFileName
     * @param query
     * @param target
     * @param smsd
     */
    protected void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, Isomorphism smsd) {

        ImageGenerator imageGenerator = new ImageGenerator();

        ////set the format right for the Tanimoto score (only two digits printed)
        NumberFormat nf = getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        LOGGER.debug("Output of the final Mappings: ");
        int counter = 1;
        for (AtomAtomMapping mapping : smsd.getAllAtomMapping()) {

            String tanimoto = nf.format(smsd.getTanimotoSimilarity());
            String stereo = "NA";
            if (smsd.getStereoScore(counter - 1) != null) {
                stereo = nf.format(smsd.getStereoScore(counter - 1));
            }
            String label = "Scores [" + "Tanimoto: " + tanimoto + ", Stereo: " + stereo + "]";
            try {
                imageGenerator.addImages(query, target, label, mapping);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            counter++;
        }
        String filePNG = getProperty("user.dir") + separator + outPutFileName;
        imageGenerator.createImage(filePNG, "Query", "Target");
    }

    /**
     *
     * @param mh
     * @param substrateIndex
     * @param productIndex
     * @throws CloneNotSupportedException
     * @throws IOException
     * @throws CDKException
     */
    protected void printSMILES(Holder mh, int substrateIndex, int productIndex)
            throws CloneNotSupportedException, IOException, CDKException {
        ReactionContainer reactionStructureInformation = mh.getReactionContainer();
        String cdkSmilesE = new CDKSMILES(reactionStructureInformation.getEduct(substrateIndex), false, false).getCanonicalSMILES();
        String cdkSmilesP = new CDKSMILES(reactionStructureInformation.getProduct(productIndex), false, false).getCanonicalSMILES();

        StringBuilder sb = new StringBuilder();
        sb.append("A: ").append(reactionStructureInformation.getEduct(substrateIndex).getID()).append(" ").append(cdkSmilesE)
                .append(" B: ").append(reactionStructureInformation.getProduct(productIndex).getID()).append(" ").append(cdkSmilesP).append(NEW_LINE);

        sb.append("A: ").append(reactionStructureInformation.getEduct(substrateIndex).getAtomCount())
                .append(" B: ").append(reactionStructureInformation.getProduct(productIndex).getAtomCount()).append(NEW_LINE);

        sb.append(" GetValue: ").append(mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex))
                .append(", ").append(mh.getStereoMatrix().getValue(substrateIndex, productIndex));

        LOGGER.debug(sb.toString());
    }
}
