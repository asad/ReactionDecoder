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
/**
 * @RCSfile: GameTheoryMixture.java
 *
 * @Author: Syed Asad Rahman @Date: 2004/06/3 @Revision: 1.10
 *
 * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
 *
 * @Contact: asad@ebi.ac.uk
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
package uk.ac.ebi.reactionblast.mapping.algorithm;

//~--- non-JDK imports --------------------------------------------------------
import static java.lang.System.out;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.mapping.algorithm.checks.ChooseWinner;
import uk.ac.ebi.reactionblast.mapping.algorithm.checks.ReactionIsomorphismHandler;
import uk.ac.ebi.reactionblast.mapping.algorithm.checks.RuleBasedMappingHandler;
import uk.ac.ebi.reactionblast.mapping.algorithm.checks.Selector;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;
import uk.ac.ebi.reactionblast.mapping.container.ReactionContainer;
import uk.ac.ebi.reactionblast.mapping.container.helper.MolMapping;
import uk.ac.ebi.reactionblast.mapping.graph.GraphMatching;
import uk.ac.ebi.reactionblast.mapping.interfaces.AbstractGraphMatching;
import uk.ac.ebi.reactionblast.tools.CDKSMILES;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import uk.ac.ebi.reactionblast.tools.labelling.SmilesMoleculeLabeller;
import static java.util.Collections.synchronizedList;

final class GameTheoryMixture extends BaseGameTheory {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 1808979786969868698L;
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
        if (DEBUG) {
            out.println("I am MIXTURE");
        }
        this.canonLabeler = new SmilesMoleculeLabeller();
        this.removeHydrogen = removeHydrogen;
        this.reaction = reaction;
        this.educts = _educts;
        this.products = _products;
        this.RID = reaction.getID();
        this.eductList = synchronizedList(rpsh.getEductCounter());
        this.productList = synchronizedList(rpsh.getProductCounter());
        this.mh = rpsh.getMatrixHolder();

        setReactionMolMapping(rpsh.getReactionMolMapping());

        winner = new ChooseWinner(eductList, productList);
        this._dirSuffix = super.getSuffix();

        ReactionIsomorphismHandler RIH = new ReactionIsomorphismHandler(mh, eductList, productList);
//        System.out.println("\n--------------\nprintAtomContainerSize");
//        printMatrixAtomContainer(mh, eductList, productList);

        if (RIH.getIsomorphismFlag()) {
//            System.out.println("ISOMORPHISM");
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

    private synchronized void GenerateIsoMorphismMapping() throws Exception {

        winner.searchWinners(educts, products, mh);

        if (winner.getFlag()) {

//            System.out.println("**********Updated Mapping**************");
            UpdateMapping();
//            System.out.println("**********Updated Matrix**************");
            UpdateMatrix(mh, removeHydrogen);
//            System.out.println("**********Generate Mapping**************");
            GenerateMapping(false);
        }
    }

    private synchronized void GenerateMapping(boolean flag) throws Exception {
        boolean ruleMatchingFlag = flag;
        if (DEBUG) {
            printMatrixAtomContainer(mh, eductList, productList);
            out.println("**********Orignal Matrix**************");
            printSimMatrix(mh, eductList, productList);
            printCliqueMatrix(mh, eductList, productList);
//        printStereoMatrix(mh, eductList, productList);
//        printFragmentMatrix(mh, eductList, productList);
//        printEnergyMatrix(mh, eductList, productList);
        }

        if (!ruleMatchingFlag) {
            RuleBasedMappingHandler ruleBasedMappingHandler = new RuleBasedMappingHandler(mh, eductList, productList);
            if (ruleBasedMappingHandler.isMatchFound()) {
                if (DEBUG) {
                    out.println("Rule Based Mapping Handler Match Found");
                }
                mh = Selector.modifyMatrix(ruleBasedMappingHandler.getMatrixHolder());
            }
            ruleMatchingFlag = true;
        }

        winner.searchWinners(educts, products, mh);
        if (DEBUG) {
            printFlagMatrix(winner, eductList, productList);
        }
        if (winner.getFlag()) {

            if (DEBUG) {
                System.out.println("**********Updated Mapping**************");
            }
            UpdateMapping();
            if (DEBUG) {
                System.out.println("**********Updated Matrix**************");
            }
            UpdateMatrix(mh, removeHydrogen);
            if (DEBUG) {
                System.out.println("**********Generate Mapping**************");
            }
            GenerateMapping(ruleMatchingFlag);
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
                    BitSet A = reactionStructureInformationContainer.getFingerPrintofEduct(substrateIndex);    // A=EDUCT
                    BitSet B = reactionStructureInformationContainer.getFingerPrintofProduct(productIndex);    // B=PRODUCT

                    /*
                     * Choose this function if you want JMCS to run
                     */
                    ac1.setID(this.eductList.get(substrateIndex));
                    ac2.setID(this.productList.get(productIndex));

                    AbstractGraphMatching GM = new GraphMatching(RID, ac1, ac2, _dirSuffix, removeHydrogen);
                    boolean mcsMatch = GM.mcsMatch(mh, removeHydrogen, substrateIndex, productIndex, A, B);
//                    System.out.println("Mol Size E: " + ac1.getAtomCount() + " , Mol Size P: " + ac2.getAtomCount());
                    if (mcsMatch) {
//                        System.out.println(eductList.get(substrateIndex) + " <=> " + productList.get(productIndex));
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
