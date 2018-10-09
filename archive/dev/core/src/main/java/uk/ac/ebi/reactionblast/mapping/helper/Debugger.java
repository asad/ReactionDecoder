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

import static java.io.File.separator;
import java.io.IOException;
import static java.lang.System.err;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import static java.text.NumberFormat.getInstance;
import java.util.List;
import java.util.Map;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Logger.getLogger;
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
            out.println("<--------Atom Size in the Container-------->");
            for (int i = 0; i < EdMap.size(); i++) {
                out.println("Educt " + EdMap.get(i) + " : " + _rSTMap.getEduct(i).getAtomCount());
                if (!_rSTMap.getEduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getEduct(i), true, false);
                    out.println("SMILES: " + sm.getCanonicalSMILES());
                }
                printAtoms(_rSTMap.getEduct(i));
            }
            out.println();
            for (int i = 0; i < PdMap.size(); i++) {
                out.println("Product " + PdMap.get(i) + " : " + _rSTMap.getProduct(i).getAtomCount());
                if (!_rSTMap.getProduct(i).isEmpty()) {
                    CDKSMILES sm = new CDKSMILES(_rSTMap.getProduct(i), true, false);
                    out.println("SMILES: " + sm.getCanonicalSMILES());
                }
                printAtoms(_rSTMap.getProduct(i));
            }
            out.println();
        } catch (IOException | CDKException | CloneNotSupportedException ex) {
            getLogger(BaseGameTheory.class.getName()).log(SEVERE, null, ex);
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
        out.println();
        out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            out.println("Clique Matrix");
            out.print("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                out.print(" " + PdMap.get(j) + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            out.println();
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                out.print(" " + EdMap.get(i) + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getCliqueMatrix().getValue(i, j);
                    result = format.format(val);
                    out.print("   " + result);
                }
                out.println();
            }
        } catch (IOException | CDKException e) {
            err.println("Parser Error");
            e.printStackTrace();
        }
        out.println();
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
        out.println();
        out.println("********* MATRIX **********");
        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;
            out.println("Similarity Matrix");
            out.print("\t\t");
            for (int j = 0; j < PdMap.size(); j++) {
                out.print(" " + PdMap.get(j) + ":(" + reactionStructureInformationContainer.getProduct(j).getAtomCount() + ")");
            }
            out.println();
            double val;
            for (int i = 0; i < EdMap.size(); i++) {
                out.print(" " + EdMap.get(i) + ":(" + reactionStructureInformationContainer.getEduct(i).getAtomCount() + ")");
                for (int j = 0; j < PdMap.size(); j++) {
                    val = mh.getGraphSimilarityMatrix().getValue(i, j);
                    result = format.format(val);
                    out.print("   " + result);
                }
                out.println();
            }
        } catch (IOException | CDKException e) {
            err.println("Parser Error");
            e.printStackTrace();
        }
        out.println();

    }

    /**
     *
     * @param Max
     * @param EdMap
     * @param PdMap
     */
    protected void printFlagMatrix(ChooseWinner Max, List<String> EdMap, List<String> PdMap) {

        out.println();
        out.println("********* MATRIX **********");
        boolean[][] FlagMatrix = Max.getFlagMatrix();
        out.println("Flag Matrix");
        out.print("\t\t");
        for (String PdMap1 : PdMap) {
            out.print("  " + PdMap1 + " ");
        }

        out.println();
        for (int i = 0; i < EdMap.size(); i++) {
            out.print(" " + EdMap.get(i));
            for (int j = 0; j < PdMap.size(); j++) {
                if (FlagMatrix[i][j]) {
                    out.print("      " + 1 + "  ");
                } else {
                    out.print("      " + 0 + "  ");
                }

            }
            out.println();
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

        out.println();
        out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            out.println("Stereo Matrix");
            out.print("\t\t");

            for (String PdMap1 : PdMap) {
                out.print(" " + PdMap1);
            }

            out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                out.print(" " + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = StereoMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    out.print("   " + result);
                }

                out.println();
            }

        } catch (Exception e) {
            err.println("Parser Error");
            e.printStackTrace();
        }

        out.println();
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printFragmentMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix FragmentMatrix = mh.getFragmentMatrix();

        out.println();
        out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            out.println("Fragment Matrix");
            out.print("\t\t");

            for (String PdMap1 : PdMap) {
                out.print(" " + PdMap1);
            }

            out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                out.print(" " + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = FragmentMatrix.getValue(i, j);
                    result
                            = format.format(val);
                    out.print("   " + result);
                }

                out.println();
            }

        } catch (Exception e) {
            err.println("Parser Error");
            e.printStackTrace();
        }

        out.println();
    }

    /**
     *
     * @param mh
     * @param EdMap
     * @param PdMap
     */
    protected void printEnergyMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
        EBIMatrix energyMatrixProfile = mh.getEnergyMatrix();

        out.println();
        out.println("********* MATRIX **********");

        try {
            NumberFormat format = new DecimalFormat("0.00");
            String result;

            out.println("Energy Matrix");
            out.print("\t\t");

            for (String PdMap1 : PdMap) {
                out.print("\t" + PdMap1);
            }

            out.println();
            double val;
            for (int i = 0; i
                    < EdMap.size(); i++) {
                out.print("\t" + EdMap.get(i));
                for (int j = 0; j
                        < PdMap.size(); j++) {
                    val = energyMatrixProfile.getValue(i, j);
                    result
                            = format.format(val);
                    out.print("\t" + result);
                }

                out.println();
            }

        } catch (Exception e) {
            err.println("Parser Error");
            e.printStackTrace();
        }

        out.println();
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
        out.println("Output of the final Mappings: ");
        out.println("Mol1: " + mol1.getID());
        out.println("Mol2: " + mol2.getID());
        try {
            if (comparison.getMappingCount() > 0) {

                for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                    int final_solution_size = final_solution.getCount();
                    out.println("Final mapping Nr. " + ++count_final_sol
                            + " Size:" + final_solution_size);

                    for (Map.Entry<IAtom, IAtom> mapping : final_solution.getMappingsByAtoms().entrySet()) {
                        IAtom eAtom = mapping.getKey();
                        IAtom pAtom = mapping.getValue();

                        out.println((mol1.indexOf(eAtom) + 1) + " " + (mol2.indexOf(pAtom) + 1));

                        out.println(eAtom.getSymbol() + " "
                                + pAtom.getSymbol());
                    }
                    out.println("");

                    out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol - 1));
                    out.println("Stereo different: " + comparison.isStereoMisMatch());
                    out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol - 1));
                }

                out.println("");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
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

        out.println("A: " + reactionStructureInformation.getEduct(substrateIndex).getID() + " " + cdkSmilesE
                + " B: " + reactionStructureInformation.getProduct(productIndex).getID() + " " + cdkSmilesP);

        out.println("A: " + reactionStructureInformation.getEduct(substrateIndex).getAtomCount()
                + " B: " + reactionStructureInformation.getProduct(productIndex).getAtomCount());

        out.println(
                " GetValue: " + mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex)
                + ", " + mh.getStereoMatrix().getValue(substrateIndex, productIndex));
    }
}
