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
import java.io.Serializable;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;
import org.openscience.cdk.AtomContainer;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogens;

import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public final class RuleBasedMappingHandler implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static boolean DEBUG1 = false;
    private final static boolean DEBUG2 = false;
    private static final long serialVersionUID = 88765671L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(RuleBasedMappingHandler.class);

    /*
     * Flags
     */
    private boolean ruleMatched = false;
    private Holder matrixHolder;
    private Holder matrixHolderClone;
    private final Map<Integer, Integer> matchedRowColoumn;

    /*
     * SMARTS parser
     */
    private IAtomContainer smartsWater;
    private IAtomContainer smartsPhosphate;
    private IAtomContainer smartsSulphate;
    private IAtomContainer smartsL_Glutamate;
    private IAtomContainer smartsL_Glutamine;
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

    /**
     *
     * @param matrixHolder
     * @param EdMapOrignal
     * @param PdMapOrignal
     * @throws CDKException
     * @throws IOException
     */
    public RuleBasedMappingHandler(Holder matrixHolder, List<String> EdMapOrignal, List<String> PdMapOrignal) throws CDKException, IOException {
        setRulesSmiles();
        if (DEBUG1) {
            out.println("Mapping Rules Checked");
        }
        this.matrixHolder = matrixHolder;
        this.matchedRowColoumn = new HashMap<>();
        setRuleMatched(false);

        int smallestMatchedReactant = Integer.MAX_VALUE;
        int smallestMatchedProduct = Integer.MAX_VALUE;
        for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
            IAtomContainer ac1 = new AtomContainer(this.matrixHolder.getReactionContainer().getEduct(i));
            ac1 = removeHydrogens(ac1);
            if (DEBUG1) {
                out.println("Educt " + unique().create(ac1));
            }

            if (ac1.getAtomCount() >= getSmartsPhosphate().getAtomCount()
                    || ac1.getAtomCount() >= getSmartsSulphate().getAtomCount()) {
                if (isMatch(getSmartsPhosphate(), ac1) || isMatch(getSmartsSulphate(), ac1)) {
                    if (smallestMatchedReactant > ac1.getAtomCount()) {
                        smallestMatchedReactant = ac1.getAtomCount();
                    }
                }
            }
        }
        if (DEBUG1) {
            out.println("smallestMatchedReactant " + smallestMatchedReactant);
        }
        for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
            IAtomContainer ac2 = new AtomContainer(this.matrixHolder.getReactionContainer().getProduct(j));
            ac2 = removeHydrogens(ac2);
            if (DEBUG1) {
                out.println("Product " + unique().create(ac2));
            }
            if (ac2.getAtomCount() >= getSmartsPhosphate().getAtomCount()
                    || ac2.getAtomCount() >= getSmartsSulphate().getAtomCount()) {
                if (isMatch(getSmartsPhosphate(), ac2) || isMatch(getSmartsSulphate(), ac2)) {
                    if (smallestMatchedProduct > ac2.getAtomCount()) {
                        smallestMatchedProduct = ac2.getAtomCount();
                    }
                }
            }
        }
        if (DEBUG1) {
            out.println("smallestMatchedProduct " + smallestMatchedProduct);
            out.println(NEW_LINE + NEW_LINE + "----------------------" + NEW_LINE + NEW_LINE + NEW_LINE);
        }
        try {
            for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                IAtomContainer educt = this.matrixHolder.getReactionContainer().getEduct(i);
                IAtomContainer ac1 = new AtomContainer(educt);
                ac1 = removeHydrogens(ac1);
                if (DEBUG2) {
                    out.println(NEW_LINE + NEW_LINE + NEW_LINE + "Educt " + unique().create(ac1));
                    out.println("Educt found " + ac1.getAtomCount());
                }

                for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                    IAtomContainer product = this.matrixHolder.getReactionContainer().getProduct(j);
                    IAtomContainer ac2 = new AtomContainer(product);
                    ac2 = removeHydrogens(ac2);

                    if (DEBUG2) {
                        out.println("Product " + unique().create(ac2));
                        out.println("Product found " + ac2.getAtomCount());
                    }
                    if (DEBUG2) {
                        out.println("Match 1 " + isMatch(getSmartsWater(), ac1));
                        out.println("Match 2 " + isMatch(getSmartsPhosphate(), ac2));
                        out.println("Query " + ac1.getAtomCount());
                        out.println("Target " + ac2.getAtomCount());
                        out.println("smallest R  " + smallestMatchedReactant);
                        out.println("smallest P  " + smallestMatchedProduct);
                    }

                    if (this.matrixHolder.getCliqueMatrix().getValue(i, j) == 0) {
                        continue;
                    }

                    /*
                    Rule 1_A water and Phosphate
                     */
                    if (ac1.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac1)
                            && isMatch(getSmartsPhosphate(), ac2)
                            && ac2.getAtomCount() == smallestMatchedProduct) {
                        if (DEBUG2) {
                            out.println("Match ");
                            out.println("smallest R phosphate " + smallestMatchedReactant);
                            out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println(" Rule 1 water and Phosphate");
                        }

                    } else /*
                        Rule 1_B phophate and water
                     */ if (ac2.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac2)
                            && isMatch(getSmartsPhosphate(), ac1)
                            && ac1.getAtomCount() == smallestMatchedReactant) {
                        if (DEBUG2) {
                            out.println("Match ");
                            out.println("smallest R phosphate " + smallestMatchedReactant);
                            out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println(" Rule 1 phosphate and water");
                        }

                    }
                    /*
                    Rule 1_C water and Sulphate
                     */
                    if (ac1.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac1)
                            && isMatch(getSmartsSulphate(), ac2)
                            && ac2.getAtomCount() == smallestMatchedProduct) {
                        if (DEBUG2) {
                            out.println("Match ");
                            out.println("smallest R phosphate " + smallestMatchedReactant);
                            out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println(" Rule 1 water and Sulphate");
                        }

                    } else /*
                        Rule 1_D Sulphate and water
                     */ if (ac2.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac2)
                            && isMatch(getSmartsSulphate(), ac1)
                            && ac1.getAtomCount() == smallestMatchedReactant) {
                        if (DEBUG2) {
                            out.println("Match ");
                            out.println("smallest R phosphate " + smallestMatchedReactant);
                            out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println(" Rule 1 Sulphate and water");
                        }

                    }/*
                        Rule 2 L_Glutamate and L_Glutamine
                     */ else if ((ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsGlutamate(), ac1) && isMatch(getSmartsGlutamine(), ac2))
                            || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsGlutamine(), ac1) && isMatch(getSmartsGlutamate(), ac2))) {
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println("Rule 2 L-Glutamate with L-Glutamine found");
                        }
                    } /*
                        Rule 3 D_Glutamate and TwoOxoglutarate
                     */ else if ((ac2.getAtomCount() == 10 && ac1.getAtomCount() == 10
                            && isMatch(getSmartsTwoOxoglutarate(), ac2) && isMatch(getSmartsD_Glutamate(), ac1))
                            || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsTwoOxoglutarate(), ac1) && isMatch(getSmartsD_Glutamate(), ac2))) {

                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println("Rule 3 D-Glutamate with 2-Oxoglutarate found");
                        }

                    }/*
                        Rule 4 water and Acetate (exact match)
                     */ else if ((ac1.getAtomCount() == 1 && isMatch(getSmartsWater(), ac1)
                            && ac2.getAtomCount() == getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac2))
                            || (ac2.getAtomCount() == 1 && isMatch(getSmartsWater(), ac2)
                            && ac1.getAtomCount() == getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac1))) {
                        if (DEBUG1) {
                            out.println("Rule 4 Water and Acetate found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }/*
                        Rule 5 ADP_ATP
                     */ else if ((ac1.getAtomCount() == getSmartsATP().getAtomCount() && isMatch(getSmartsATP(), ac1)
                            && isMatch(getSmartsADP(), ac2))
                            || (ac1.getAtomCount() == getSmartsADP().getAtomCount() && isMatch(getSmartsADP(), ac1)
                            && isMatch(getSmartsATP(), ac2))) {
                        if (DEBUG1) {
                            out.println("Rule 5 ADP_ATP found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }/*
                        Rule 6 CoA_Acetyl_CoA
                     */ else if ((ac1.getAtomCount() == getSmartsCoA().getAtomCount() && isMatch(getSmartsCoA(), ac1)
                            && isMatch(getSmartsAcetyl_CoA(), ac2))
                            || (ac1.getAtomCount() == getSmartsAcetyl_CoA().getAtomCount() && isMatch(getSmartsAcetyl_CoA(), ac1)
                            && isMatch(getSmartsCoA(), ac2))) {
                        if (DEBUG1) {
                            out.println("Rule 6 CoA_Acetyl_CoA found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }/*
                        Rule 7 C00003_C00006
                     */ else if ((ac1.getAtomCount() == getSmartsC00003().getAtomCount() && isMatch(getSmartsC00003(), ac1)
                            && isMatch(getSmartsC00006(), ac2))
                            || (ac1.getAtomCount() == getSmartsC00006().getAtomCount() && isMatch(getSmartsC00006(), ac1)
                            && isMatch(getSmartsC00003(), ac2))) {
                        if (DEBUG1) {
                            out.println("Rule 7 C00003_C00006 found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }/*
                        Rule 8 C00004_C00005
                     */ else if ((ac1.getAtomCount() == getSmartsC00004().getAtomCount() && isMatch(getSmartsC00004(), ac1)
                            && isMatch(getSmartsC00005(), ac2))
                            || (ac1.getAtomCount() == getSmartsC00005().getAtomCount() && isMatch(getSmartsC00005(), ac1)
                            && isMatch(getSmartsC00004(), ac2))) {
                        if (DEBUG1) {
                            out.println("Rule 8 C00004_C00005 found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    } /*
                        Rule 9 C00022_C00041
                     */ else if ((ac1.getAtomCount() == getSmartsPyruvate().getAtomCount() && isMatch(getSmartsPyruvate(), ac1)
                            && isMatch(getSmartsAlanine(), ac2))
                            || (ac1.getAtomCount() == getSmartsAlanine().getAtomCount() && isMatch(getSmartsAlanine(), ac1)
                            && isMatch(getSmartsPyruvate(), ac2))) {
                        if (DEBUG1) {
                            out.println("Rule 9 C00022_C00041 found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }/*
                        Rule 10 N_C
                     */ else if (isMatch(getSmartsNRule(), ac1) && isMatch(getSmartsCRule(), ac2)
                            || (isMatch(getSmartsCRule(), ac1) && isMatch(getSmartsNRule(), ac2))) {
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG1) {
                            out.println("Rule 10 N with C found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
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
                LOGGER.debug("ERROR: Matrix Holder clone error");
                LOGGER.error(SEVERE, null, ex);
            }
            for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                    /*
                    * Keep the mapping for the smart matches
                     */
                    if (this.matchedRowColoumn.containsKey(i) && this.matchedRowColoumn.get(i) == j) {
                        /*
                        * reset the mapping for the unmatched
                         */
                        matrixHolderClone.getGraphSimilarityMatrix().setValue(i, j, 1.00);
                    } else /*
                        * reset the mapping for the unmatched
                     */ {
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
            /*
            * Assign the new matrix
             */
            this.matrixHolder = matrixHolderClone;
        }
    }

    /**
     * @return the matrixHolder
     */
    public synchronized Holder getMatrixHolder() {
        return matrixHolder;
    }

    /**
     * @return the ruleMatched
     */
    public synchronized boolean isMatchFound() {
        return ruleMatched;
    }

    /**
     * @param ruleMatched the ruleMatched to set
     */
    private synchronized void setRuleMatched(boolean ruleMatched) {
        this.ruleMatched = ruleMatched;
    }

    private boolean isMatch(IAtomContainer ac1, IAtomContainer ac2) {
        if (ac1.getAtomCount() <= ac2.getAtomCount()) {
            try {
                Substructure s = new Substructure(ac1, ac2, true, true, false, false);
                if (DEBUG2) {
                    out.println("ac1 " + ac1.getAtomCount());
                    out.println("ac2 " + ac2.getAtomCount());
                    out.println("Sub " + s.isSubgraph());
                    out.println("score " + s.getTanimotoSimilarity());
                }
                return s.isSubgraph();
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        } else {
            try {
                Substructure s = new Substructure(ac2, ac1, true, true, false, false);
                if (DEBUG2) {
                    out.println("ac2 " + ac2.getAtomCount());
                    out.println("ac1 " + ac1.getAtomCount());
                    out.println("Sub " + s.isSubgraph());
                    out.println("score " + s.getTanimotoSimilarity());
                }
                return s.isSubgraph();
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        return false;
    }

    private void setRulesSmiles() throws CDKException {
        /*
         * water with phosphate
         */

        final String waterSMILES = "O";
        final String phosphateSMILES = "OP(O)(O)=O";

        /*
         * Acetate C00033
         */
        final String acetate = "CC(O)=O";
        final String sulphateSMILES = "O=S(=O)(O)O";
        /*
         * L-Glutamate with L-Glutamine
         */
        final String lGlutamate = "N[C@@H](CCC(O)=O)C(O)=O";
        final String lGlutamine = "N[C@@H](CCC(N)=O)C(O)=O";
        /*
         * 2-Oxoglutarate to D-Glutamate
         */
        final String twoOxoglutarate = "OC(=O)CCC(=O)C(O)=O";
        final String dGlutamate = "N[C@H](CCC(O)=O)C(O)=O";


        /*
         * ADP_ATP
         */
        final String ATP = "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O";
        final String ADP = "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O";

        /*
         * CoA_Acetyl_CoA
         */
        final String CoA = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N)[C@@H](O)C(=O)NCCC(=O)NCCS";
        final String Acetyl_CoA = "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N";


        /*
         * C00003_C00006
         */
        final String C00003 = "NC(=O)C1=CC=C[N+](=C1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C(N)N=CN=C23)[C@@H](O)[C@H]1O";
        final String C00006 = "NC(=O)C1=C[N+](=CC=C1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](OP(O)(O)=O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O";

        /*
         * C00004_C00005
         */
        final String C00004 = "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O";
        final String C00005 = "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](OP(O)(O)=O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O";

        /*
         * C00022_C00041
         */
        final String C00022 = "[CH3][C](=O)C(O)=O";
        final String C00041 = "[CH3][C](N)C(O)=O";

        /*
         * N_C CC(C)[C@H](N)C(O)=O>>CC(C)C(=O)C(O)=O
         */
        final String N = "CC(C)[C@H](N)C(O)=O";
        final String C = "CC(C)C(=O)C(O)=O";

        SmilesParser smilesParser = new SmilesParser(getInstance());
        /*
         * Rule 1
         */
        smartsWater = smilesParser.parseSmiles(waterSMILES);
        smartsPhosphate = smilesParser.parseSmiles(phosphateSMILES);
        smartsSulphate = smilesParser.parseSmiles(sulphateSMILES);

        /*
         * Rule 2
         */
        smartsL_Glutamate = smilesParser.parseSmiles(lGlutamate);
        smartsL_Glutamine = smilesParser.parseSmiles(lGlutamine);

        /*
         * Rule 3
         */
        smartsTwoOxoglutarate = smilesParser.parseSmiles(twoOxoglutarate);
        smartsD_Glutamate = smilesParser.parseSmiles(dGlutamate);

        /*
         * Rule 4 water tends to attack acetate when Phophate is not present
         */
        smartsAcetate = smilesParser.parseSmiles(acetate);

        /*
         * Rule 5 AT_ADP
         */
        smartsATP = smilesParser.parseSmiles(ATP);
        smartsADP = smilesParser.parseSmiles(ADP);

        /*
         * Rule 6 CoA_Acetyl_CoA
         */
        smartsCoA = smilesParser.parseSmiles(CoA);
        smartsAcetyl_CoA = smilesParser.parseSmiles(Acetyl_CoA);

        /*
         * Rule 7 C00003_C00006
         */
        smartsC00003 = smilesParser.parseSmiles(C00003);
        smartsC00006 = smilesParser.parseSmiles(C00006);

        /*
         * Rule 8 C00004_C00005
         */
        smartsC00004 = smilesParser.parseSmiles(C00004);
        smartsC00005 = smilesParser.parseSmiles(C00005);

        /*
         * Rule 9 Alanine_Pyruvate
         */
        smartsAlanine = smilesParser.parseSmiles(C00041);
        smartsPyruvate = smilesParser.parseSmiles(C00022);

        /*
         * Rule 10 Valine_Isoleucine
         */
        smartsNRule = smilesParser.parseSmiles(N);
        smartsCRule = smilesParser.parseSmiles(C);

    }

    /**
     * @return the smartsWater
     */
    private IAtomContainer getSmartsWater() {
        return smartsWater;
    }

    /**
     * @return the smartsPhosphate
     */
    private IAtomContainer getSmartsPhosphate() {
        return smartsPhosphate;
    }

    /**
     * @return the smartsL_Glutamate
     */
    private IAtomContainer getSmartsGlutamate() {
        return smartsL_Glutamate;
    }

    /**
     * @return the smartsL_Glutamine
     */
    private IAtomContainer getSmartsGlutamine() {
        return smartsL_Glutamine;
    }

    /**
     * @return the smartsTwoOxoglutarate
     */
    private IAtomContainer getSmartsTwoOxoglutarate() {
        return smartsTwoOxoglutarate;
    }

    /**
     * @return the smartsD_Glutamate
     */
    private IAtomContainer getSmartsD_Glutamate() {
        return smartsD_Glutamate;
    }

    /**
     * @return the smartsAcetate
     */
    private IAtomContainer getSmartsAcetate() {
        return smartsAcetate;
    }

    /**
     * @return the smartsSulphate
     */
    private IAtomContainer getSmartsSulphate() {
        return smartsSulphate;
    }

    /**
     * @return the smartsATP
     */
    public IAtomContainer getSmartsATP() {
        return smartsATP;
    }

    /**
     * @return the smartsADP
     */
    public IAtomContainer getSmartsADP() {
        return smartsADP;
    }

    /**
     * @return the smartsCoA
     */
    public IAtomContainer getSmartsCoA() {
        return smartsCoA;
    }

    /**
     * @return the smartsAcetyl_CoA
     */
    public IAtomContainer getSmartsAcetyl_CoA() {
        return smartsAcetyl_CoA;
    }

    /**
     * @return the smartsC00003
     */
    public IAtomContainer getSmartsC00003() {
        return smartsC00003;
    }

    /**
     * @return the smartsC00006
     */
    public IAtomContainer getSmartsC00006() {
        return smartsC00006;
    }

    /**
     * @return the smartsC00004
     */
    public IAtomContainer getSmartsC00004() {
        return smartsC00004;
    }

    /**
     * @return the smartsC00005
     */
    public IAtomContainer getSmartsC00005() {
        return smartsC00005;
    }

    /**
     * @return the smartsPyruvate
     */
    public IAtomContainer getSmartsPyruvate() {
        return smartsPyruvate;
    }

    /**
     * @return the smartsAlanine
     */
    public IAtomContainer getSmartsAlanine() {
        return smartsAlanine;
    }

    /**
     * @return the smartsNRule
     */
    public IAtomContainer getSmartsNRule() {
        return smartsNRule;
    }

    /**
     * @return the smartsCRule
     */
    public IAtomContainer getSmartsCRule() {
        return smartsCRule;
    }

}
