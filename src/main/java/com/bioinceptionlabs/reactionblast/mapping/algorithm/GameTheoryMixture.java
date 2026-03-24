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
/**
 * @RCSfile: GameTheoryMixture.java
 *
 * @Author: Syed Asad Rahman @Date: 2004/06/3 @Revision: 1.10
 *
 * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
 *
 * @Contact: asad.rahman@bioinceptionlabs.com
 *
 * @This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version. All we ask is that proper credit is given for our
 * work, which includes - but is not limited to - adding the above copyright
 * notice to the beginning of your source code files, and to any copyright
 * notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *
 */
package com.bioinceptionlabs.reactionblast.mapping.algorithm;

//~--- non-JDK imports --------------------------------------------------------
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.ChooseWinner;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.ReactionIsomorphismHandler;
import com.bioinceptionlabs.reactionblast.mapping.container.MoleculeMoleculeMapping;
import com.bioinceptionlabs.reactionblast.mapping.container.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.container.helper.MolMapping;
import com.bioinceptionlabs.reactionblast.mapping.graph.GraphMatching;
import com.bioinceptionlabs.reactionblast.mapping.interfaces.AbstractGraphMatching;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import com.bioinceptionlabs.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.tools.labelling.SmilesMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.RuleBasedMappingHandler;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.Selector;

final class GameTheoryMixture extends BaseGameTheory {

    private static final int MAX_MAPPING_ITERATIONS = 100;
    private static final long serialVersionUID = 1808979786969868698L;
    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(GameTheoryMixture.class);
    private final List<String> eductList;
    private final List<String> productList;
    private Holder mh;
    private final ChooseWinner winner;
    private final IReaction reaction;
    private final String RID;
    private final String _dirSuffix;
    private final boolean removeHydrogen;
    private MoleculeMoleculeMapping reactionMolMapping = null;
    private Map<Integer, IAtomContainer> educts = null;
    private Map<Integer, IAtomContainer> products = null;
    private int delta = 0;
    private Integer stepIndex = 0;
    private final ICanonicalMoleculeLabeller canonLabeler;

    //~--- constructors -------------------------------------------------------
    GameTheoryMixture(
            IReaction reaction,
            boolean removeHydrogen,
            Map<Integer, IAtomContainer> _educts,
            Map<Integer, IAtomContainer> _products,
            GameTheoryMatrix rpsh)
            throws Exception {
        LOGGER.debug("I am MIXTURE");
        this.canonLabeler = new SmilesMoleculeLabeller();
        this.removeHydrogen = removeHydrogen;
        this.reaction = reaction;
        this.educts = _educts;
        this.products = _products;
        this.RID = reaction.getID();
        this.eductList = new ArrayList<>(rpsh.getEductCounter());
        this.productList = new ArrayList<>(rpsh.getProductCounter());
        this.mh = rpsh.getMatrixHolder();

        setReactionMolMapping(rpsh.getReactionMolMapping());

        winner = new ChooseWinner(eductList, productList);
        this._dirSuffix = super.getSuffix();

        ReactionIsomorphismHandler RIH = new ReactionIsomorphismHandler(mh, eductList, productList);
        if (RIH.getIsomorphismFlag()) {
//            LOGGER.debug("ISOMORPHISM");
//            printSimMatrix(mh, eductList, productList);
//            printCliqueMatrix(mh, eductList, productList);
//            printStereoMatrix();
//            printFragmentMatrix();
//            printEnergyMatrix(mh, eductList, productList);
            mh = RIH.getMatrixHolder();
            GenerateIsoMorphismMapping();
        } else {
            GenerateMapping(false);
        }
    }
//~--- methods ------------------------------------------------------------

    private void GenerateIsoMorphismMapping() throws Exception {

        winner.searchWinners(educts, products, mh);

        if (winner.getFlag()) {

//            LOGGER.debug("**********Updated Mapping**************");
            UpdateMapping();
//            LOGGER.debug("**********Updated Matrix**************");
            UpdateMatrix(mh, removeHydrogen);
//            LOGGER.debug("**********Generate Mapping**************");
            GenerateMapping(false);
        }
    }

    private void GenerateMapping(boolean flag) throws Exception {
        boolean ruleMatchingFlag = flag;
        int iteration = 0;
        boolean continueMapping = true;
        while (continueMapping && iteration < MAX_MAPPING_ITERATIONS) {
            if (!ruleMatchingFlag) {//First map the biggest fragment the call rules
                RuleBasedMappingHandler ruleBasedMappingHandler
                        = new RuleBasedMappingHandler(mh, eductList, productList);
                if (ruleBasedMappingHandler.isMatchFound()) {
                    LOGGER.debug("Rule Based Mapping Handler Match Found");
                    mh = Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                }
                ruleMatchingFlag = true;
            }

            winner.searchWinners(educts, products, mh);
            if (winner.getFlag()) {

                LOGGER.debug("**********Updated Mapping**************");
                UpdateMapping();

                // Early termination: check if all atoms are already mapped
                ReactionContainer rc = mh.getReactionContainer();
                boolean allMapped = true;
                for (int i = 0; i < rc.getEductCount() && allMapped; i++) {
                    if (rc.getEduct(i).getAtomCount() > 0) allMapped = false;
                }
                for (int j = 0; j < rc.getProductCount() && allMapped; j++) {
                    if (rc.getProduct(j).getAtomCount() > 0) allMapped = false;
                }
                if (allMapped) {
                    LOGGER.debug("Early termination: perfect match — all atoms mapped");
                    break;
                }

                // Early termination if no remaining mappable pairs
                boolean hasRemainingPairs = false;
                for (int i = 0; i < mh.getGraphSimilarityMatrix().getRowDimension(); i++) {
                    for (int j = 0; j < mh.getGraphSimilarityMatrix().getColumnDimension(); j++) {
                        if (mh.getGraphSimilarityMatrix().getValue(i, j) > 0) {
                            hasRemainingPairs = true;
                            break;
                        }
                    }
                    if (hasRemainingPairs) break;
                }
                if (!hasRemainingPairs) {
                    LOGGER.debug("Early termination: no remaining mappable pairs");
                    break;
                }

                LOGGER.debug("**********Updated Matrix**************");
                UpdateMatrix(mh, removeHydrogen);
                LOGGER.debug("**********Generate Mapping**************");
                iteration++;
            } else {
                continueMapping = false;
            }
        }
    }

    private void UpdateMapping() throws Exception {
        boolean[][] FlagMatrix = winner.getFlagMatrix();

        ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();

        for (int iIndex = 0; iIndex < reactionStructureInformationContainer.getEductCount(); iIndex++) {
            for (int jIndex = 0; jIndex < reactionStructureInformationContainer.getProductCount(); jIndex++) {
                int substrateIndex = iIndex;
                int productIndex = jIndex;
                IAtomContainer ac1 = reactionStructureInformationContainer.getEduct(substrateIndex);
                IAtomContainer ac2 = reactionStructureInformationContainer.getProduct(productIndex);

                if (FlagMatrix[substrateIndex][productIndex]) {

                    // updateFlag=true;
                    BitSet A = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);    // A=EDUCT
                    BitSet B = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);    // B=PRODUCT

                    /*
                     * Choose this function if you want JMCS to run
                     */
                    ac1.setID(this.eductList.get(substrateIndex));
                    ac2.setID(this.productList.get(productIndex));

                    AbstractGraphMatching GM = new GraphMatching(RID, ac1, ac2, _dirSuffix, removeHydrogen);
                    boolean mcsMatch = GM.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, A, B);
//                    LOGGER.debug("Mol Size E: " + ac1.getAtomCount() + " , Mol Size P: " + ac2.getAtomCount());
                    if (mcsMatch) {
//                        LOGGER.debug(eductList.get(substrateIndex) + " <=> " + productList.get(productIndex));
                        delta += GM.removeMatchedAtomsAndUpdateAAM(reaction);
                        List<MolMapping> rMap = getReactionMolMapping().
                                getMapping(RID, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                        for (MolMapping map : rMap) {
                            map.setReactionMapping(true);
                            IAtomContainer mol = GM.getMatchedPart();
                            mol = canonLabeler.getCanonicalMolecule(mol);
                            CDKSMILES cdkSmiles = new CDKSMILES(mol, true, false);
                            map.setMatchedSMILES(cdkSmiles.getCanonicalSMILES(), ++stepIndex);
                        }
                    }
                    IAtomContainer RemainingEduct = GM.getRemainingEduct();
                    IAtomContainer RemainingProduct = GM.getRemainingProduct();

                    reactionStructureInformationContainer.putEduct(substrateIndex, RemainingEduct);
                    reactionStructureInformationContainer.putProduct(productIndex, RemainingProduct);
                    reactionStructureInformationContainer.setEductModified(substrateIndex, true);
                    reactionStructureInformationContainer.setProductModified(productIndex, true);
                }
            }
        }
    }

    /**
     * @return the reactionMolMapping
     */
    @Override
    public MoleculeMoleculeMapping getReactionMolMapping() {
        return reactionMolMapping;
    }

    /**
     * @param reactionMolMapping the reactionMolMapping to set
     */
    @Override
    public void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) {
        this.reactionMolMapping = reactionMolMapping;
    }

    /**
     * @return the delta
     */
    @Override
    public int getDelta() {
        return delta;
    }
}
