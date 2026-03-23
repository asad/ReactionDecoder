/*
 * Copyright (C) 2003-2020 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
 * @Author: Syed Asad Rahman
 * @Date: 2004/06/3
 * @Revision: 1.10
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
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.ChooseWinner;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.ReactionIsomorphismHandler;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.RuleBasedMappingHandler;
import com.bioinceptionlabs.reactionblast.mapping.container.MoleculeMoleculeMapping;
import com.bioinceptionlabs.reactionblast.mapping.container.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.container.helper.MolMapping;
import com.bioinceptionlabs.reactionblast.mapping.graph.GraphMatching;
import com.bioinceptionlabs.reactionblast.mapping.interfaces.AbstractGraphMatching;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import com.bioinceptionlabs.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.tools.labelling.SmilesMoleculeLabeller;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.Selector;
import static java.util.Collections.synchronizedList;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.checks.MinSelection;

final class GameTheoryMin extends BaseGameTheory {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 1808979786969868698L;
    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(GameTheoryMin.class);
    private final List<String> eductList;
    private final List<String> productList;
    private Holder mh;
    private final ChooseWinner winner;
    private final IReaction reaction;
    private final String reactionName;
    private final String _dirSuffix;
    private final boolean removeHydrogen;
    private MoleculeMoleculeMapping reactionMolMapping = null;
    private Map<Integer, IAtomContainer> educts = null;
    private Map<Integer, IAtomContainer> products = null;
    private int delta = 0;
    private Integer stepIndex = 0;
    private final ICanonicalMoleculeLabeller canonLabeler;

    private int counter = 0;

    //~--- constructors -------------------------------------------------------
    GameTheoryMin(
            IReaction reaction,
            boolean removeHydrogen,
            Map<Integer, IAtomContainer> _educts,
            Map<Integer, IAtomContainer> _products,
            GameTheoryMatrix rpsh)
            throws Exception {
        LOGGER.debug("I am MIN MIX");
        this.counter = 1;
        this.canonLabeler = new SmilesMoleculeLabeller();
        this.removeHydrogen = removeHydrogen;
        this.reaction = reaction;
        this.educts = _educts;
        this.products = _products;
        this.reactionName = reaction.getID();
        this.eductList = synchronizedList(rpsh.getEductCounter());
        this.productList = synchronizedList(rpsh.getProductCounter());
        this.mh = rpsh.getMatrixHolder();

        setReactionMolMapping(rpsh.getReactionMolMapping());

        winner = new ChooseWinner(eductList, productList);
        this._dirSuffix = super.getSuffix();

        ReactionIsomorphismHandler RIH = new ReactionIsomorphismHandler(mh, eductList, productList);

        if (RIH.getIsomorphismFlag()) {
            LOGGER.debug("ISOMORPHISM");
            mh = RIH.getMatrixHolder();
            GenerateIsoMorphismMapping();
        } else {
            GenerateMapping(false);
        }
    }

    private synchronized void GenerateIsoMorphismMapping() throws Exception {

        winner.searchWinners(educts, products, mh);

        if (winner.getFlag()) {

            LOGGER.debug("**********Updated Mapping**************");
            UpdateMapping();
            LOGGER.debug("**********Updated Matrix**************");
            UpdateMatrix(mh, removeHydrogen);
            LOGGER.debug("**********Generate Mapping**************");
            GenerateMapping(false);
        }
    }

    private synchronized void GenerateMapping(boolean flag) throws Exception {
        boolean ruleMatchingFlag = flag;
        int maxIterations = 100;
        int iteration = 0;
        boolean continueMapping = true;
        while (continueMapping && iteration < maxIterations) {
            this.counter++;
            if (DEBUG) {
                LOGGER.debug("**********Orignal Matrix**************");
                printMatrixAtomContainer(mh, eductList, productList);
                printSimMatrix(mh, eductList, productList);
                printCliqueMatrix(mh, eductList, productList);
//            printStereoMatrix(mh, eductList, productList);
//            printFragmentMatrix(mh, eductList, productList);
//            printEnergyMatrix(mh, eductList, productList);
            }
            boolean conditionmet = false;
            if (!ruleMatchingFlag) {
                LOGGER.debug("CHECK Rule Based Mapping Handler Match");
                RuleBasedMappingHandler ruleBasedMappingHandler = new RuleBasedMappingHandler(mh, eductList, productList);
                if (ruleBasedMappingHandler.isMatchFound()) {
                    LOGGER.debug("Rule Based Mapping Handler Match Found");
                    mh = Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
                    conditionmet = true;
                }
                ruleMatchingFlag = true;
                LOGGER.debug("DONE CHECK Rule Based Mapping Handler");
            }

            if (!conditionmet && counter <= 5) {
                LOGGER.debug("call counter " + counter);
                LOGGER.debug("Subgraph/Exact Match Test");
                MinSelection select
                        = new MinSelection(mh, eductList, productList);
                if (select.isSubAndCompleteMatchFlag()) {
                    LOGGER.debug("Subgraph/Exact Match");
                    mh = select.getUpdatedHolder();
                }
            }

            if (DEBUG) {
                LOGGER.debug("**********Modified Matrix**************");
//            printMatrixAtomContainer(mh, eductList, productList);
                printSimMatrix(mh, eductList, productList);
                printCliqueMatrix(mh, eductList, productList);
//            printStereoMatrix(mh, eductList, productList);
//            printFragmentMatrix(mh, eductList, productList);
//            printEnergyMatrix(mh, eductList, productList);
            }
            winner.searchWinners(educts, products, mh);
            if (DEBUG) {
                printFlagMatrix(winner, eductList, productList);
            }
            if (winner.getFlag()) {
                LOGGER.debug("**********Updated Mapping**************");
                UpdateMapping();
                LOGGER.debug("**********Updated Matrix**************");
                UpdateMatrix(mh, removeHydrogen);
                LOGGER.debug("**********Generate Mapping**************");
                iteration++;
            } else {
                continueMapping = false;
            }
        }
    }

    private synchronized void UpdateMapping() throws Exception {
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
                    BitSet a_BitSet = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);    // a_BitSet=EDUCT
                    BitSet b_BitSet = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);    // b_BitSet=PRODUCT

                    /*
                     * Choose this function if you want JMCS to run
                     */
                    ac1.setID(this.eductList.get(substrateIndex));
                    ac2.setID(this.productList.get(productIndex));

                    AbstractGraphMatching graphMatching = new GraphMatching(reactionName, ac1, ac2, _dirSuffix, removeHydrogen);
                    boolean mcsMatch = graphMatching.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, a_BitSet, b_BitSet);
                    LOGGER.debug("Mol Size E: " + ac1.getAtomCount() + " , Mol Size P: " + ac2.getAtomCount());
                    if (mcsMatch) {
                        LOGGER.debug(eductList.get(substrateIndex) + " <=> " + productList.get(productIndex));
                        delta += graphMatching.removeMatchedAtomsAndUpdateAAM(reaction);
                        List<MolMapping> rMap = getReactionMolMapping().
                                getMapping(reactionName, this.eductList.get(substrateIndex), this.productList.get(productIndex));
                        for (MolMapping map : rMap) {
                            map.setReactionMapping(true);
                            IAtomContainer mol = graphMatching.getMatchedPart();
                            mol = canonLabeler.getCanonicalMolecule(mol);
                            CDKSMILES cdkSmiles = new CDKSMILES(mol, true, false);
                            map.setMatchedSMILES(cdkSmiles.getCanonicalSMILES(), ++stepIndex);
                        }
                    }
                    IAtomContainer remainingEduct = graphMatching.getRemainingEduct();
                    IAtomContainer remainingProduct = graphMatching.getRemainingProduct();

                    reactionStructureInformationContainer.putEduct(substrateIndex, remainingEduct);
                    reactionStructureInformationContainer.putProduct(productIndex, remainingProduct);
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
    public synchronized MoleculeMoleculeMapping getReactionMolMapping() {
        return reactionMolMapping;
    }

    /**
     * @param reactionMolMapping the reactionMolMapping to set
     */
    @Override
    public synchronized void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping) {
        this.reactionMolMapping = reactionMolMapping;
    }

    /**
     * @return the delta
     */
    @Override
    public synchronized int getDelta() {
        return delta;
    }
}
