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
package uk.ac.ebi.reactionblast.mapping.helper;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.IAtomMapping;
import uk.ac.ebi.reactionblast.mapping.algorithm.BaseGameTheory;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.algorithm.checks.ChooseWinner;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import uk.ac.ebi.reactionblast.tools.CDKSMILES;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;
import uk.ac.ebi.reactionblast.tools.ImageGenerator;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class Debugger extends BasicDebugger {

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
            System.out.println("<--------Atom Size in the Container-------->");
            for (int i = 0; i < EdMap.size(); i++) {
                System.out.println("Educt " + EdMap.get(i) + " : " + _rSTMap.getEduct(i).getAtomCount());
                if (!_rSTMap.getEduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getEduct(i), true, false);
                    System.out.println("SMILES: " + sm.getCanonicalSMILES());
                }
                printAtoms(_rSTMap.getEduct(i));
            }
            System.out.println();
            for (int i = 0; i < PdMap.size(); i++) {
                System.out.println("Product " + PdMap.get(i) + " : " + _rSTMap.getProduct(i).getAtomCount());
                if (!_rSTMap.getProduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getProduct(i), true, false);
                    System.out.println("SMILES: " + sm.getCanonicalSMILES());
                }
                printAtoms(_rSTMap.getProduct(i));
            }
            System.out.println();
        } catch (IOException | CDKException | CloneNotSupportedException ex) {
            Logger.getLogger(BaseGameTheory.class.getName()).log(Level.SEVERE, null, ex);
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
        System.out.println();
        System.out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            System.out.println("Clique Matrix");
            System.out.print("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                System.out.print(" " + PdMap.get(j) + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            System.out.println();
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                System.out.print(" " + EdMap.get(i) + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getCliqueMatrix().getValue(i, j);
                    result = format.format(val);
                    System.out.print("   " + result);
                }
                System.out.println();
            }
        } catch (IOException | CDKException e) {
            System.err.println("Parser Error");
            e.printStackTrace();
        }
        System.out.println();
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
        System.out.println();
        System.out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            System.out.println("Similarity Matrix");
            System.out.print("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                System.out.print(" " + PdMap.get(j) + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            System.out.println();
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                System.out.print(" " + EdMap.get(i) + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getGraphSimilarityMatrix().getValue(i, j);
                    result = format.format(val);
                    System.out.print("   " + result);
                }
                System.out.println();
            }
        } catch (IOException | CDKException e) {
            System.err.println("Parser Error");
            e.printStackTrace();
        }
        System.out.println();

    }

    /**
     *
     * @param Max
     * @param EdMap
     * @param PdMap
     */
    protected void printFlagMatrix(ChooseWinner Max, List<String> EdMap, List<String> PdMap) {

        System.out.println();
        System.out.println("********* MATRIX **********");
        boolean[][] FlagMatrix = Max.getFlagMatrix();
        System.out.println("Flag Matrix");
        System.out.print("\t\t");
        for (String PdMap1 : PdMap) {
            System.out.print("  " + PdMap1 + " ");
        }

        System.out.println();
        for (int i = 0; i < EdMap.size(); i++) {
            System.out.print(" " + EdMap.get(i));
            for (int j = 0; j < PdMap.size(); j++) {
                if (FlagMatrix[i][j]) {
                    System.out.print("      " + 1 + "  ");
                } else {
                    System.out.print("      " + 0 + "  ");
                }

            }
            System.out.println();
        }

    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printStereoMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix StereoMatrix = mh.getStereoMatrix();

        System.out.println();
        System.out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            System.out.println("Stereo Matrix");
            System.out.print("\t\t");

            for (String PdMap1 : PdMap) {
                System.out.print(" " + PdMap1);
            }

            System.out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                System.out.print(" " + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = StereoMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    System.out.print("   " + result);
                }

                System.out.println();
            }

        } catch (Exception e) {
            System.err.println("Parser Error");
            e.printStackTrace();
        }

        System.out.println();
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printFragmentMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix FragmentMatrix = mh.getFragmentMatrix();

        System.out.println();
        System.out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            System.out.println("Fragment Matrix");
            System.out.print("\t\t");

            for (String PdMap1 : PdMap) {
                System.out.print(" " + PdMap1);
            }

            System.out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                System.out.print(" " + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = FragmentMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    System.out.print("   " + result);
                }

                System.out.println();
            }

        } catch (Exception e) {
            System.err.println("Parser Error");
            e.printStackTrace();
        }

        System.out.println();
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printEnergyMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix energyMatrixProfile = mh.getEnergyMatrix();

        System.out.println();
        System.out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            System.out.println("Energy Matrix");
            System.out.print("\t\t");

            for (String PdMap1 : PdMap) {
                System.out.print("\t" + PdMap1);
            }

            System.out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                System.out.print("\t" + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = energyMatrixProfile.getValue(i, j);
                    result
                            = format.format(val);
                    System.out.print("\t" + result);
                }

                System.out.println();
            }

        } catch (Exception e) {
            System.err.println("Parser Error");
            e.printStackTrace();
        }

        System.out.println();
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
        System.out.println("Output of the final Mappings: ");
        System.out.println("Mol1: " + mol1.getID());
        System.out.println("Mol2: " + mol2.getID());
        try {
            if (comparison.getMappingCount() > 0) {

                for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                    int final_solution_size = final_solution.getCount();
                    System.out.println("Final mapping Nr. " + ++count_final_sol
                            + " Size:" + final_solution_size);

                    for (Map.Entry<IAtom, IAtom> mapping : final_solution.getMappingsByAtoms().entrySet()) {
                        IAtom eAtom = mapping.getKey();
                        IAtom pAtom = mapping.getValue();

                        System.out.println((mol1.getAtomNumber(eAtom) + 1) + " " + (mol2.getAtomNumber(pAtom) + 1));

                        System.out.println(eAtom.getSymbol() + " "
                                + pAtom.getSymbol());
                    }
                    System.out.println("");

                    System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol - 1));
                    System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol - 1));
                }

                System.out.println("");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
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
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        System.out.println("Output of the final Mappings: ");
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
                Logger.getLogger(Debugger.class.getName()).log(Level.SEVERE, null, ex);
            }
            counter++;
        }
        String filePNG = System.getProperty("user.dir") + File.separator + outPutFileName;
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

        System.out.println("A: " + reactionStructureInformation.getEduct(substrateIndex).getID() + " " + cdkSmilesE
                + " B: " + reactionStructureInformation.getProduct(productIndex).getID() + " " + cdkSmilesP);

        System.out.println("A: " + reactionStructureInformation.getEduct(substrateIndex).getAtomCount()
                + " B: " + reactionStructureInformation.getProduct(productIndex).getAtomCount());

        System.out.println(
                " GetValue: " + mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex)
                + ", " + mh.getStereoMatrix().getValue(substrateIndex, productIndex));
    }
}
