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
package uk.ac.ebi.reactionblast.mapping.algorithm.checks;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public final class RuleBasedMappingHandler implements Serializable {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 88765671L;

    /*
     * Flags
     */
    private boolean ruleMatched = false;
    private Holder matrixHolder;
    private Holder matrixHolderClone;
    private final Map<Integer, Integer> matchedRowColoumn;

    /**
     *
     * @param matrixHolder
     * @param EdMapOrignal
     * @param PdMapOrignal
     * @throws CDKException
     * @throws IOException
     */
    public RuleBasedMappingHandler(
            Holder matrixHolder,
            List<String> EdMapOrignal,
            List<String> PdMapOrignal) throws CDKException, IOException {
        setRulesSmiles();
        if (DEBUG) {
            System.out.println("Mapping Rules Checked");
        }
        this.matrixHolder = matrixHolder;
        this.matchedRowColoumn = new HashMap<>();
        setRuleMatched(false);

        int smallestMatchedReactant = Integer.MAX_VALUE;
        int smallestMatchedProduct = Integer.MAX_VALUE;
        for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
            IAtomContainer ac1 = new AtomContainer(this.matrixHolder.getReactionContainer().getEduct(i));
            ac1 = ExtAtomContainerManipulator.removeHydrogens(ac1);
            if (DEBUG) {
                System.out.println("Educt " + SmilesGenerator.unique().create(ac1));
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
        if (DEBUG) {
            System.out.println("smallestMatchedReactant " + smallestMatchedReactant);
        }
        for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
            IAtomContainer ac2 = new AtomContainer(this.matrixHolder.getReactionContainer().getProduct(j));
            ac2 = ExtAtomContainerManipulator.removeHydrogens(ac2);
            if (DEBUG) {
                System.out.println("Product " + SmilesGenerator.unique().create(ac2));
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
        if (DEBUG) {
            System.out.println("smallestMatchedProduct " + smallestMatchedProduct);
            System.out.println("\n\n----------------------\n\n\n");
        }
        try {
            for (int i = 0; i < this.matrixHolder.getReactionContainer().getEductCount(); i++) {
                IAtomContainer educt = this.matrixHolder.getReactionContainer().getEduct(i);
                IAtomContainer ac1 = new AtomContainer(educt);
                ac1 = ExtAtomContainerManipulator.removeHydrogens(ac1);
                if (DEBUG) {
                    System.out.println("\n\n\nEduct " + SmilesGenerator.unique().create(ac1));
                    System.out.println("Educt found " + ac1.getAtomCount());
                }

                for (int j = 0; j < this.matrixHolder.getReactionContainer().getProductCount(); j++) {
                    IAtomContainer product = this.matrixHolder.getReactionContainer().getProduct(j);
                    IAtomContainer ac2 = new AtomContainer(product);
                    ac2 = ExtAtomContainerManipulator.removeHydrogens(ac2);

                    if (DEBUG) {
                        System.out.println("Product " + SmilesGenerator.unique().create(ac2));
                        System.out.println("Product found " + ac2.getAtomCount());
                    }
                    if (DEBUG) {
                        System.out.println("Match 1 Water " + isMatch(getSmartsWater(), ac1));
                        System.out.println("Match 2 Phos " + isMatch(getSmartsPhosphate(), ac2));
                        System.out.println("Water " + ac1.getAtomCount());
                        System.out.println("Phos " + ac2.getAtomCount());
                        System.out.println("smallest R phosphate " + smallestMatchedReactant);
                        System.out.println("smallest P phosphate " + smallestMatchedProduct);
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
                        if (DEBUG) {
                            System.out.println("Match ");
                            System.out.println("smallest R phosphate " + smallestMatchedReactant);
                            System.out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println(" Rule 1 water and Phosphate");
                        }

                    } else /*
                     Rule 1_B phophate and water
                     */ if (ac2.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac2)
                            && isMatch(getSmartsPhosphate(), ac1)
                            && ac1.getAtomCount() == smallestMatchedReactant) {
                        if (DEBUG) {
                            System.out.println("Match ");
                            System.out.println("smallest R phosphate " + smallestMatchedReactant);
                            System.out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println(" Rule 1 phosphate and water");
                        }

                    }
                    /*
                     Rule 1_C water and Sulphate
                     */
                    if (ac1.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac1)
                            && isMatch(getSmartsSulphate(), ac2)
                            && ac2.getAtomCount() == smallestMatchedProduct) {
                        if (DEBUG) {
                            System.out.println("Match ");
                            System.out.println("smallest R phosphate " + smallestMatchedReactant);
                            System.out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println(" Rule 1 water and Sulphate");
                        }

                    } else /*
                     Rule 1_D Sulphate and water
                     */ if (ac2.getAtomCount() == 1
                            && isMatch(getSmartsWater(), ac2)
                            && isMatch(getSmartsSulphate(), ac1)
                            && ac1.getAtomCount() == smallestMatchedReactant) {
                        if (DEBUG) {
                            System.out.println("Match ");
                            System.out.println("smallest R phosphate " + smallestMatchedReactant);
                            System.out.println("smallest P phosphate " + smallestMatchedProduct);
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println(" Rule 1 Sulphate and water");
                        }

                    }/*
                     Rule 2 L_Glutamate and L_Glutamine
                     */ else if ((ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsL_Glutamate(), ac1) && isMatch(getSmartsL_Glutamine(), ac2))
                            || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsL_Glutamine(), ac1) && isMatch(getSmartsL_Glutamate(), ac2))) {
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println("Rule 2 L-Glutamate with L-Glutamine found");
                        }
                    } /*
                     Rule 3 D_Glutamate and TwoOxoglutarate
                     */ else if ((ac2.getAtomCount() == 10 && ac1.getAtomCount() == 10
                            && isMatch(getSmartsTwoOxoglutarate(), ac2) && isMatch(getSmartsD_Glutamate(), ac1))
                            || (ac1.getAtomCount() == 10 && ac2.getAtomCount() == 10
                            && isMatch(getSmartsTwoOxoglutarate(), ac1) && isMatch(getSmartsD_Glutamate(), ac2))) {

                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                        if (DEBUG) {
                            System.out.println("Rule 3 D-Glutamate with 2-Oxoglutarate found");
                        }

                    }/*
                     Rule 4 water and Acetate
                     */ else if ((ac1.getAtomCount() == 1 && isMatch(getSmartsWater(), ac1)
                            && ac2.getAtomCount() >= getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac2))
                            || (ac2.getAtomCount() == 1 && isMatch(getSmartsWater(), ac2)
                            && ac1.getAtomCount() >= getSmartsAcetate().getAtomCount() && isMatch(getSmartsAcetate(), ac1))) {
                        if (DEBUG) {
                            System.out.println("Rule 4 Water and Acetate found");
                        }
                        setRuleMatched(true);
                        matchedRowColoumn.put(i, j);
                    }
                }
            }
        } catch (IOException | CDKException ex) {
            Logger.getLogger(RuleBasedMappingHandler.class.getName()).
                    log(Level.WARNING, "Error in Matching Rules", ex);
        }
        if (this.isMatchFound()) {
            try {
                this.matrixHolderClone = (Holder) matrixHolder.clone();
            } catch (CloneNotSupportedException ex) {
                System.err.println("ERROR: Matrix Holder clone error");
                Logger.getLogger(RuleBasedMappingHandler.class.getName()).log(Level.SEVERE, null, ex);
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
                        matrixHolderClone.getStereoMatrix().setValue(i, j, -Double.MAX_VALUE);
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
                if (DEBUG) {
                    System.out.println("ac1 " + ac1.getAtomCount());
                    System.out.println("ac2 " + ac2.getAtomCount());
                    System.out.println("Sub " + s.isSubgraph());
                    System.out.println("score " + s.getTanimotoSimilarity());
                }
                return s.isSubgraph();
            } catch (CDKException ex) {
                Logger.getLogger(RuleBasedMappingHandler.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            try {
                Substructure s = new Substructure(ac2, ac1, true, true, false, false);
                if (DEBUG) {
                    System.out.println("ac2 " + ac2.getAtomCount());
                    System.out.println("ac1 " + ac1.getAtomCount());
                    System.out.println("Sub " + s.isSubgraph());
                    System.out.println("score " + s.getTanimotoSimilarity());
                }
                return s.isSubgraph();
            } catch (CDKException ex) {
                Logger.getLogger(RuleBasedMappingHandler.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return false;
    }


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

    public void setRulesSmiles() throws CDKException {
        /*
         * water with phosphate
         */

        final String waterSMILES = "O";
        final String phosphateSMILES = "OP(O)(O)=O";

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
         Acetate C00033
         */
        final String acetate = "CC(O)=O";

        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
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
    private IAtomContainer getSmartsL_Glutamate() {
        return smartsL_Glutamate;
    }

    /**
     * @return the smartsL_Glutamine
     */
    private IAtomContainer getSmartsL_Glutamine() {
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
    private static final Logger LOG = Logger.getLogger(RuleBasedMappingHandler.class.getName());

}
