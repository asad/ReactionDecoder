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

import java.util.*;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockPair;
import uk.ac.ebi.reactionblast.mapping.blocks.MappingGraph;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;

public class RBlastReaction {

    private IReaction reaction;
    private List<BlockPair> blockPairs;
    private List<IBond> bondsCleavedInReactant;
    private List<IBond> bondsFormedInProduct;
    private List<IBond> bondsOrderChangedInProduct;
    private List<IBond> bondsOrderChangedInReactant;
    private List<IBond> bondsStereoChangedInProduct;
    private List<IBond> bondsStereoChangedInReactant;
    private Map<IAtom, IStereoAndConformation> atomStereoProductMap;
    private Map<IAtom, IStereoAndConformation> atomStereoReactantMap;

    public RBlastReaction(IReaction reaction, boolean useMappingGraph) {
        this.reaction = reaction;
        if (useMappingGraph) {
            MappingGraph mappingGraph = new MappingGraph(reaction);
            blockPairs = mappingGraph.createBlockPairs(reaction);
        } else {
            blockPairs = new ArrayList<>();
        }
        setupBondChangeLists(reaction);
    }

    public RBlastReaction(IReaction reaction, BlockPair pair) {
        this.reaction = reaction;
        blockPairs = new ArrayList<>();
        blockPairs.add(pair);
        setupBondChangeLists(reaction);
    }

    public RBlastReaction(IReaction reaction, List<BlockPair> pairs) {
        this.reaction = reaction;
        blockPairs = new ArrayList<>();
        blockPairs.addAll(pairs);
        setupBondChangeLists(reaction);
    }

    private void setupBondChangeLists(IReaction reaction) {
        bondsCleavedInReactant = new ArrayList<>();
        bondsFormedInProduct = new ArrayList<>();
        bondsOrderChangedInReactant = new ArrayList<>();
        bondsOrderChangedInProduct = new ArrayList<>();
        bondsStereoChangedInReactant = new ArrayList<>();
        bondsStereoChangedInProduct = new ArrayList<>();
        atomStereoProductMap = new HashMap<>();
        atomStereoReactantMap = new HashMap<>();
        setupBondChangeLists(reaction.getReactants(), true);
        setupBondChangeLists(reaction.getProducts(), false);
        setupAtomChangeLists(reaction.getReactants(), atomStereoProductMap);
        setupAtomChangeLists(reaction.getProducts(), atomStereoReactantMap);
    }

    private void setupAtomChangeLists(IAtomContainerSet molSet, Map<IAtom, IStereoAndConformation> mapToAddTo) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            for (IAtom atom : atomContainer.atoms()) {
                if (atom.getProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION) != null) {
                    IStereoAndConformation flag = atom.getProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION);
                    if (flag == IStereoAndConformation.R) {
                        mapToAddTo.put(atom, IStereoAndConformation.R);
                    } else if (flag == IStereoAndConformation.S) {
                        mapToAddTo.put(atom, IStereoAndConformation.S);
                    }
                    if (flag == IStereoAndConformation.Z) {
                        mapToAddTo.put(atom, IStereoAndConformation.Z);
                    } else if (flag == IStereoAndConformation.E) {
                        mapToAddTo.put(atom, IStereoAndConformation.E);
                    }
                }
            }
            for (IBond bond : atomContainer.bonds()) {
                for (IAtom atom : bond.atoms()) {
                    if (atom.getProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION) != null) {
                        IStereoAndConformation flag = atom.getProperty(ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION);
                        if (flag == IStereoAndConformation.R) {
                            mapToAddTo.put(atom, IStereoAndConformation.R);
                        } else if (flag == IStereoAndConformation.S) {
                            mapToAddTo.put(atom, IStereoAndConformation.S);
                        }
                        if (flag == IStereoAndConformation.Z) {
                            mapToAddTo.put(atom, IStereoAndConformation.Z);
                        } else if (flag == IStereoAndConformation.E) {
                            mapToAddTo.put(atom, IStereoAndConformation.E);
                        }
                    }
                }
            }
        }
    }

    private void setupBondChangeLists(IAtomContainerSet molSet, boolean isReactant) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            setupBondChangeLists(atomContainer, isReactant);
        }
    }

    private void setupBondChangeLists(IAtomContainer atomContainer, boolean isReactant) {
        for (IBond bond : atomContainer.bonds()) {
            if (bond.getProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) != null) {
                ECBLAST_BOND_CHANGE_FLAGS bondChangeType = bond.getProperty(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION);
                if (bondChangeType == ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED) {
                    bondsCleavedInReactant.add(bond);
                } else if (bondChangeType == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED) {
                    bondsFormedInProduct.add(bond);
                } else if (bondChangeType == ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER) {
                    if (isReactant) {
                        bondsOrderChangedInReactant.add(bond);
                    } else {
                        bondsOrderChangedInProduct.add(bond);
                    }
                } else if (bondChangeType == ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO) {
                    if (isReactant) {
                        bondsStereoChangedInReactant.add(bond);
                    } else {
                        bondsStereoChangedInProduct.add(bond);
                    }
                } else if (bondChangeType == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED_OR_CLEAVED) {
                    if (isReactant) {
                        bondsCleavedInReactant.add(bond);
                    } else {
                        bondsFormedInProduct.add(bond);
                    }
                }
            }
        }
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsStereoChangedInProduct() {
        return Collections.unmodifiableList(bondsStereoChangedInProduct);
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsStereoChangedInReactant() {
        return Collections.unmodifiableList(bondsStereoChangedInReactant);
    }

    /**
     *
     * @return
     */
    public IReaction getReaction() {
        return reaction;
    }

    /**
     *
     * @param reaction
     */
    public void setReaction(IReaction reaction) {
        this.reaction = reaction;
    }

    /**
     *
     * @return
     */
    public List<BlockPair> getMappedSubgraphs() {
        return Collections.unmodifiableList(blockPairs);
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsCleavedInReactant() {
        return Collections.unmodifiableList(bondsCleavedInReactant);
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsFormedInProduct() {
        return Collections.unmodifiableList(bondsFormedInProduct);
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsOrderChangedInProduct() {
        return Collections.unmodifiableList(bondsOrderChangedInProduct);
    }

    /**
     *
     * @return
     */
    public List<IBond> getBondsOrderChangedInReactant() {
        return Collections.unmodifiableList(bondsOrderChangedInReactant);
    }

    /**
     *
     * @return
     */
    public Map<IAtom, IStereoAndConformation> getAtomStereoProductMap() {
        return Collections.unmodifiableMap(atomStereoProductMap);
    }

    /**
     *
     * @return
     */
    public Map<IAtom, IStereoAndConformation> getAtomStereoReactantMap() {
        return Collections.unmodifiableMap(atomStereoReactantMap);
    }
}
