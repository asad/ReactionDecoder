/*
 * ReactionContainer - consolidated mapping container classes.
 * Merged: Key, MolMapping, BestMatchContainer, CDKReactionBuilder, HydrogenFreeFingerPrintContainer, MoleculeMoleculeMapping into ReactionContainer
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.FingerprintGenerator;
import com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.IFingerprintGenerator;
import com.bioinceptionlabs.reactionblast.mapping.BestMatch;
import com.bioinceptionlabs.reactionblast.tools.MoleculeTools.AtomContainerSetComparator;
import com.bioinceptionlabs.reactionblast.tools.MoleculeTools.BasicDebugger;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import com.bioinception.smsd.core.SearchEngine;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.Substructure;
import static com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.FingerprintGenerator.getFingerprinterSize;
import static com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.Similarity.getTanimotoSimilarity;
import static java.util.Collections.sort;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.smsd.ExtAtomContainerManipulator.aromatizeMolecule;
import static org.openscience.smsd.ExtAtomContainerManipulator.cloneWithIDs;
import static org.openscience.smsd.ExtAtomContainerManipulator.fixDativeBonds;
import static org.openscience.smsd.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogens;


/**
 *
 *
 * @Author: Syed Asad Rahman
 * @Contact: asad.rahman@bioinceptionlabs.com
 * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
 * @Date: 2004/06/3
 * @RCSfile: ReactionContainer.java,v
 * @Revision: 1.10
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
 */
/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
interface IKey extends Comparable<ReactionContainer.Key>, Comparator<ReactionContainer.Key> {

    @Override
    int compareTo(ReactionContainer.Key t);

    @Override
    boolean equals(Object o);

    @Override
    int hashCode();

    @Override
    String toString();

    @Override
    int compare(ReactionContainer.Key t, ReactionContainer.Key t1);
}

public class ReactionContainer implements Cloneable, Serializable {

    static final long serialVersionUID = 17278639972837695L;
    /*
     * Singleton Pattern Implementation
     */
    private final Map<Integer, IAtomContainer> eAtomContainerMap;
    private final Map<Integer, IAtomContainer> pAtomContainerMap;
    private final Map<Integer, BitSet> eFingerPrintMap;
    private final Map<Integer, BitSet> pFingerPrintMap;
    private final Map<Integer, Boolean> eductContainerModificationMap;
    private final Map<Integer, Boolean> productContainerModificationMap;
    private final IFingerprintGenerator fpr;

    //~--- constructors -------------------------------------------------------
    /**
     *
     * @throws Exception
     */
    public ReactionContainer() throws Exception {
        eAtomContainerMap = new TreeMap<>();
        pAtomContainerMap = new TreeMap<>();
        eFingerPrintMap = new TreeMap<>();
        pFingerPrintMap = new TreeMap<>();
        eductContainerModificationMap = new TreeMap<>();
        productContainerModificationMap = new TreeMap<>();
        fpr = new FingerprintGenerator();
    }

    @Override
    public String toString() {
        return "ReactionContainer{" + "eAtomContainerMap=" + eAtomContainerMap
                + ", pAtomContainerMap=" + pAtomContainerMap
                + ", eFingerPrintMap=" + eFingerPrintMap
                + ", pFingerPrintMap=" + pFingerPrintMap
                + ", eductContainerModificationMap=" + eductContainerModificationMap
                + ", productContainerModificationMap=" + productContainerModificationMap
                + ", fpr=" + fpr + '}';
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    public void Clear() throws IOException {
        eAtomContainerMap.clear();
        pAtomContainerMap.clear();
        eFingerPrintMap.clear();
        pFingerPrintMap.clear();
        eductContainerModificationMap.clear();
        productContainerModificationMap.clear();
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    public void eraseEduct(int Key) throws IOException {
        eAtomContainerMap.remove(Key);
        eFingerPrintMap.remove(Key);
        eductContainerModificationMap.remove(Key);
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    public void eraseProduct(int Key) throws IOException {
        pAtomContainerMap.remove(Key);
        pFingerPrintMap.remove(Key);
        productContainerModificationMap.remove(Key);
    }

    /**
     *
     * @return
     */
    public Integer getEductCount() {
        return eAtomContainerMap.size();
    }

    /**
     *
     * @return
     */
    public Integer getProductCount() {
        return pAtomContainerMap.size();
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public IAtomContainer getEduct(int value)
            throws IOException, CDKException {
        return eAtomContainerMap.containsKey(value) ? eAtomContainerMap.get(value) : null;

    }

    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public IAtomContainer getProduct(int value)
            throws IOException, CDKException {
        return pAtomContainerMap.containsKey(value) ? pAtomContainerMap.get(value) : null;
    }

    /**
     *
     * @return @throws java.io.IOException
     * @throws CDKException
     */
    public Collection<IAtomContainer> getEducts()
            throws IOException, CDKException {
        return eAtomContainerMap.values();

    }

    /**
     *
     * @return @throws java.io.IOException
     * @throws CDKException
     */
    public Collection<IAtomContainer> getProducts()
            throws IOException, CDKException {
        return pAtomContainerMap.values();
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public boolean isEductModified(int value)
            throws IOException, CDKException {
        return eductContainerModificationMap.containsKey(value)
                ? eductContainerModificationMap.get(value) : false;

    }

    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public boolean isProductModified(int value)
            throws IOException, CDKException {
        return productContainerModificationMap.containsKey(value)
                ? productContainerModificationMap.get(value) : false;
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param index
     * @param educt
     * @throws java.io.IOException
     * @throws Exception
     */
    public void putEduct(int index, IAtomContainer educt)
            throws IOException, Exception {
        eAtomContainerMap.put(index, educt);
        if (educt.getAtomCount() == 0) {
            setFingerPrintofEduct(index, new BitSet(getFingerprinterSize()));
        } else {
            BitSet fp = fpr.getFingerprint(educt);
            setFingerPrintofEduct(index, (BitSet) fp.clone());
        }
    }

    /**
     *
     * @param index
     * @param product
     * @throws java.io.IOException
     * @throws Exception
     */
    public void putProduct(int index, IAtomContainer product)
            throws IOException, Exception {
        pAtomContainerMap.put(index, product);
        if (product.getAtomCount() == 0) {
            setFingerPrintofProduct(index, new BitSet(getFingerprinterSize()));
        } else {
            BitSet fp = fpr.getFingerprint(product);
            setFingerPrintofProduct(index, (BitSet) fp.clone());
        }
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param index
     * @param flag
     * @throws java.io.IOException
     * @throws Exception
     */
    public void setEductModified(int index, boolean flag)
            throws IOException, Exception {
        eductContainerModificationMap.put(index, flag);
    }

    /**
     *
     * @param index
     * @param flag
     * @throws java.io.IOException
     * @throws Exception
     */
    public void setProductModified(int index, boolean flag)
            throws IOException, Exception {
        productContainerModificationMap.put(index, flag);
    }

    /**
     *
     * @param eudMap
     * @throws java.io.IOException
     * @throws Exception
     */
    public void putAllEduct(TreeMap<Integer, IAtomContainer> eudMap)
            throws IOException, Exception {
        eAtomContainerMap.putAll(eudMap);
        for (Map.Entry<Integer, IAtomContainer> map : eudMap.entrySet()) {
            BitSet fp = fpr.getFingerprint(map.getValue());
            setFingerPrintofEduct(map.getKey(), fp);
        }
    }

    /**
     *
     * @param prodMap
     * @throws java.io.IOException
     * @throws Exception
     */
    public void putAllProduct(TreeMap<Integer, IAtomContainer> prodMap)
            throws IOException, Exception {
        pAtomContainerMap.putAll(prodMap);
        for (Map.Entry<Integer, IAtomContainer> map : prodMap.entrySet()) {
            BitSet fp = fpr.getFingerprint(map.getValue());
            setFingerPrintofProduct(map.getKey(), fp);
        }
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public BitSet getFingerPrintofEduct(int value)
            throws IOException, CDKException {
        return eFingerPrintMap.containsKey(value) ? eFingerPrintMap.get(value) : null;
    }

    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public BitSet getFingerPrintofProduct(int value)
            throws IOException, CDKException {
        return pFingerPrintMap.containsKey(value) ? pFingerPrintMap.get(value) : null;
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param value
     * @param Edu
     * @throws java.io.IOException
     */
    private void setFingerPrintofEduct(int value, BitSet edu)
            throws IOException {
        eFingerPrintMap.put(value, edu);

    }

    /**
     *
     * @param value
     * @param Prod
     * @throws java.io.IOException
     */
    private void setFingerPrintofProduct(int value, BitSet prod) throws IOException {
        pFingerPrintMap.put(value, prod);

    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class Key implements IKey, Serializable {

        private static final long serialVersionUID = 92392372979041041L;
        private final int sourceIndex;
        private final int targetIndex;

        /**
         *
         * @param sourceIndex
         * @param targetIndex
         */
        public Key(int sourceIndex, int targetIndex) {
            this.sourceIndex = sourceIndex;
            this.targetIndex = targetIndex;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Key:");
            sb.append(this.sourceIndex);
            sb.append(":");
            sb.append(this.targetIndex);
            return sb.toString();
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final Key other = (Key) obj;
            if (this.sourceIndex != other.sourceIndex) {
                return false;
            }
            return this.targetIndex == other.targetIndex;
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 97 * hash + this.sourceIndex;
            hash = 97 * hash + this.targetIndex;
            return hash;
        }

        @Override
        public int compareTo(Key t) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;
            String key1 = this.sourceIndex + "_" + this.targetIndex;
            String key2 = t.sourceIndex + "_" + t.targetIndex;

            if (key1.equals(key2)) {
                return EQUAL;
            } else {
                return key1.compareTo(key2);
            }
        }

        /**
         *
         * @param t
         * @param t1
         * @return
         */
        @Override
        public int compare(Key t, Key t1) {
            return t.compareTo(t1);
        }
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class MolMapping extends Object implements Serializable {

        private static final long serialVersionUID = 1738327023703717L;
        private final String mol1;
        private final String mol2;
        private final Integer indexI;
        private final Integer indexJ;
        private Integer indexStep = 0;
        private boolean keggMapping;
        private boolean rBLASTMapping;
        private IAtomContainer matchedMol = null;
        private String matchedSMILES = null;

        /**
         *
         * @param mol1
         * @param mol2
         * @param indexI
         * @param indexJ
         */
        public MolMapping(String mol1, String mol2, Integer indexI, Integer indexJ) {
            this.mol1 = mol1;
            this.mol2 = mol2;
            this.indexI = indexI;
            this.indexJ = indexJ;
            this.keggMapping = false;
            this.rBLASTMapping = false;
        }

        /**
         *
         * @return
         */
        public String getQuery() {
            return mol1;
        }

        /**
         *
         * @return
         */
        public String getTarget() {
            return mol2;
        }

        /**
         * @return the keggMapping
         */
        public boolean isKeggMapping() {
            return keggMapping;
        }

        /**
         * @param keggMapping the keggMapping to set
         */
        public void setKeggMapping(boolean keggMapping) {
            this.keggMapping = keggMapping;
        }

        /**
         * @return the rBLASTMapping
         */
        public boolean isrBLASTMapping() {
            return rBLASTMapping;
        }

        /**
         * @param rBLASTMapping the rBLASTMapping to set
         */
        public void setReactionMapping(boolean rBLASTMapping) {
            this.rBLASTMapping = rBLASTMapping;
        }

        /**
         * @return the indexI
         */
        public Integer getIndexI() {
            return indexI;
        }

        /**
         * @return the indexJ
         */
        public Integer getIndexJ() {
            return indexJ;
        }

        /**
         * @return the matchedMol
         */
        public IAtomContainer getMatchedMol() {
            return matchedMol;
        }

        /**
         * @param matchedMol the matchedMol to set
         */
        public void setMatchedMol(IAtomContainer matchedMol) {
            this.matchedMol = matchedMol;
        }

        /**
         * @return the matchedSMILES
         */
        public String getMatchedSMILES() {
            return matchedSMILES;
        }

        /**
         * @param matchedSMILES the matchedSMILES to set
         * @param step index for the mapping process
         */
        public void setMatchedSMILES(String matchedSMILES, Integer step) {
            this.matchedSMILES = matchedSMILES;
            setIndexStep(step);
        }

        /**
         * @return the indexStep
         */
        public Integer getIndexStep() {
            return indexStep;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final MolMapping other = (MolMapping) obj;
            if ((this.mol1 == null) ? (other.mol1 != null) : !this.mol1.equals(other.mol1)) {
                return false;
            }
            if ((this.mol2 == null) ? (other.mol2 != null) : !this.mol2.equals(other.mol2)) {
                return false;
            }
            if (!Objects.equals(this.indexI, other.indexI) && (this.indexI == null || !this.indexI.equals(other.indexI))) {
                return false;
            }
            if (!Objects.equals(this.indexJ, other.indexJ) && (this.indexJ == null || !this.indexJ.equals(other.indexJ))) {
                return false;
            }
            if (!Objects.equals(this.indexStep, other.indexStep) && (this.indexStep == null || !this.indexStep.equals(other.indexStep))) {
                return false;
            }
            if (this.keggMapping != other.keggMapping) {
                return false;
            }
            if (this.rBLASTMapping != other.rBLASTMapping) {
                return false;
            }
            if (this.matchedMol != other.matchedMol && (this.matchedMol == null || !this.matchedMol.equals(other.matchedMol))) {
                return false;
            }
            return !((this.matchedSMILES == null) ? (other.matchedSMILES != null) : !this.matchedSMILES.equals(other.matchedSMILES));
        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 83 * hash + (this.mol1 != null ? this.mol1.hashCode() : 0);
            hash = 83 * hash + (this.mol2 != null ? this.mol2.hashCode() : 0);
            hash = 83 * hash + (this.indexI != null ? this.indexI.hashCode() : 0);
            hash = 83 * hash + (this.indexJ != null ? this.indexJ.hashCode() : 0);
            hash = 83 * hash + (this.indexStep != null ? this.indexStep.hashCode() : 0);
            hash = 83 * hash + (this.keggMapping ? 1 : 0);
            hash = 83 * hash + (this.rBLASTMapping ? 1 : 0);
            hash = 83 * hash + (this.matchedMol != null ? this.matchedMol.hashCode() : 0);
            hash = 83 * hash + (this.matchedSMILES != null ? this.matchedSMILES.hashCode() : 0);
            return hash;
        }

        /**
         * @param indexStep the indexStep to set
         */
        private void setIndexStep(Integer indexStep) {
            this.indexStep = indexStep;
        }

        @Override
        public String toString() {
            return "MolMapping{" + "mol1=" + mol1 + ", mol2=" + mol2 + ", indexI=" + indexI + ", indexJ=" + indexJ + ", indexStep=" + indexStep + ", keggMapping=" + keggMapping + ", rBLASTMapping=" + rBLASTMapping + ", matchedMol=" + matchedMol + ", matchedSMILES=" + matchedSMILES + '}';
        }
    }



    //~--- classes --------------------------------------
    /**
     *
     *
     * @Author: Syed Asad Rahman: asad.rahman@bioinceptionlabs.com
     * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
     * @Date: 2004/06/3
     * @Revision: 1.10
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
     */
    public static class BestMatchContainer extends BestMatch implements Serializable {

        private final static ILoggingTool LOGGER
                = createLoggingTool(BestMatchContainer.class);
        private static final long serialVersionUID = 10947239472059259L;
        private final Map<IKey, AtomAtomMapping> mcsAtomMap;
        private final Map<IKey, Integer> fragmentCount;
        private final Map<IKey, Double> bondBreakingEnergy;
        private final Map<IKey, Double> stereoScore;
        private final Map<IKey, Double> similarity;

        //~--- constructors -------------------------------------------------------
        /**
         *
         */
        public BestMatchContainer() {
            mcsAtomMap = new HashMap<>();
            fragmentCount = new HashMap<>();
            bondBreakingEnergy = new HashMap<>();
            stereoScore = new HashMap<>();
            similarity = new HashMap<>();
        }

        @Override
        public String toString() {
            return "BestMatchContainer{" + "mcsAtomMap=" + mcsAtomMap + ", fragmentCount=" + fragmentCount + ", bondBreakingEnergy=" + bondBreakingEnergy + ", stereoScore=" + stereoScore + ", similarity=" + similarity + '}';
        }

        //~--- methods ------------------------------------------------------------
        /**
         *
         * @throws java.io.IOException
         */
        @Override
        public void Clear() throws IOException {
            mcsAtomMap.clear();
            fragmentCount.clear();
            bondBreakingEnergy.clear();
            stereoScore.clear();
            similarity.clear();
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @throws java.io.IOException
         */
        @Override
        public void Erase(int indexI, int indexJ) throws IOException {
            IKey key = new Key(indexI, indexJ);
            mcsAtomMap.remove(key);
            fragmentCount.remove(key);
            bondBreakingEnergy.remove(key);
            stereoScore.remove(key);
            similarity.remove(key);
        }

        //~--- get methods --------------------------------------------------------
        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         * @throws IOException
         */
        @Override
        public AtomAtomMapping getAtomMatch(int indexI, int indexJ)
                throws IOException {
            IKey key = new Key(indexI, indexJ);
            if (mcsAtomMap.containsKey(key)) {
                return mcsAtomMap.get(key);
            } else {
                try {
                    throw new CDKException("Key not found:" + key + " in " + mcsAtomMap.keySet());
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
            return null;
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         */
        @Override
        public double getBondEnergy(int indexI, int indexJ) {
            IKey key = new Key(indexI, indexJ);
            return bondBreakingEnergy.containsKey(key) ? bondBreakingEnergy.get(key) : 0.;
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         */
        @Override
        public int getTotalFragmentCount(int indexI, int indexJ) {
            IKey key = new Key(indexI, indexJ);
            return fragmentCount.containsKey(key) ? fragmentCount.get(key) : 0;
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @param value
         */
        @Override
        public void setBondEnergy(int indexI, int indexJ, double value) {
            IKey key = new Key(indexI, indexJ);
            bondBreakingEnergy.put(key, value);
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @param value
         */
        @Override
        public void setTotalFragmentCount(int indexI, int indexJ, Integer value) {
            IKey key = new Key(indexI, indexJ);
            fragmentCount.put(key, value);
        }

        //~--- set methods --------------------------------------------------------
        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @param matchingAtoms
         * @throws IOException
         */
        @Override
        public void putBestMapping(int indexI, int indexJ, AtomAtomMapping matchingAtoms)
                throws IOException {
            IKey key = new Key(indexI, indexJ);
            mcsAtomMap.put(key, matchingAtoms);
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         * @throws IOException
         */
        @Override
        public boolean containsKey(int indexI, int indexJ)
                throws IOException {
            IKey key = new Key(indexI, indexJ);
            return mcsAtomMap.containsKey(key);
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @param stereoVal
         */
        @Override
        public void setStereoScore(int indexI, int indexJ, double stereoVal) {
            IKey key = new Key(indexI, indexJ);
            stereoScore.put(key, stereoVal);
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         */
        @Override
        public double getStereoScore(int indexI, int indexJ) {
            IKey key = new Key(indexI, indexJ);
            return stereoScore.containsKey(key) == true ? stereoScore.get(key) : 0.0;
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @param _simVal
         */
        @Override
        public void setGraphSimilarity(int indexI, int indexJ, double _simVal) {
            IKey key = new Key(indexI, indexJ);
            similarity.put(key, _simVal);
        }

        /**
         * String IKey = substrateIndex + "_" + productIndex;
         *
         * @param indexI
         * @param indexJ
         * @return
         */
        @Override
        public double getGraphSimilarity(int indexI, int indexJ) {
            IKey key = new Key(indexI, indexJ);
            return similarity.containsKey(key) == true ? similarity.get(key) : 0.0;
        }
    }



    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class CDKReactionBuilder extends BasicDebugger implements Serializable {

        private static final long serialVersionUID = 19869866609698L;
        private final static ILoggingTool LOGGER
                = createLoggingTool(CDKReactionBuilder.class);
        private final IReactionSet reactionSet;
        private int moleculeCounter = 0; //Counter to create Unique Molecules
        private final Map<String, Double> stoichiometryMap;
        private final Map<String, BitSet> fingerprintMap;
        private final Map<String, IAtomContainer> moleculeMap;

        /**
         *
         * @throws java.lang.Exception
         */
        public CDKReactionBuilder() throws Exception {
            reactionSet = SilentChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
            stoichiometryMap = new HashMap<>();
            fingerprintMap = new HashMap<>();
            moleculeMap = new HashMap<>();
        }

        @Override
        public String toString() {
            return "CDKReactionBuilder{" + "reactionSet=" + reactionSet + ", moleculeCounter="
                    + moleculeCounter + ", stoichiometryMap=" + stoichiometryMap + ", fingerprintMap="
                    + fingerprintMap + ", moleculeMap=" + moleculeMap + '}';
        }

        /**
         *
         * @param reactionSet
         * @throws java.lang.Exception
         */
        public void standardize(IReactionSet reactionSet) throws Exception {
            for (IReaction reaction : reactionSet.reactions()) {
                IReaction standardizedReaction = standardize(reaction);
                reactionSet.addReaction(standardizedReaction);
            }
        }

        /**
         *
         * @param reaction
         * @return
         * @throws Exception
         */
        public IReaction standardize(IReaction reaction) throws Exception {
            int old_atom_rank_index_reactant = 1;
            int old_atom_rank_index_product = 1;
            List<IAtomContainer> _metabolites = new ArrayList<>();
            IReaction standardizedReaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);

            String reactionID = reaction.getID();
            int reactionCounter = 1;
            if (reactionID == null) {
                reactionID = "R" + Long.toString(reactionCounter++);
                reaction.setID(reactionID);
            }

            _metabolites.clear();

            standardizedReaction.setID(reactionID);

            stoichiometryMap.clear();

            Double tempStoic;

            LOGGER.debug("standardize reaction module phase 1");
            for (IAtomContainer mol : reaction.getReactants().atomContainers()) {
                String id = mol.getID() == null || mol.getID().trim().isEmpty() ? null : mol.getID();
                tempStoic = 1.0;
                if (reaction.getReactantCoefficient(mol) > 0) {
                    tempStoic = reaction.getReactantCoefficient(mol);
                }

                IAtomContainer gMol = cloneWithIDs(mol);

                /*
                 * Set old Atom Index
                 */
                for (IAtom a : gMol.atoms()) {
                    if (a.getProperties() == null) {
                        a.addProperties(new HashMap<>());
                    }
                    a.setProperty("OLD_RANK", old_atom_rank_index_reactant++);
                }
                LOGGER.debug("standardize reaction module phase 1.1.1");
                fixDativeBonds(gMol);
                LOGGER.debug("standardize reaction module phase 1.1.2");
                percieveAtomTypesAndConfigureAtoms(gMol);
                IAtomContainer molWithH = gMol;
                //= ExtAtomContainerManipulator.addExplicitH(gMol);
                aromatizeMolecule(molWithH);

                LOGGER.debug(id + " standardize reaction module phase 1.2");

                if (id == null || id.isEmpty()) {
                    molWithH = setProperty(molWithH);
                } else {
                    molWithH.setID(id);
                }

                if (stoichiometryMap.containsKey(molWithH.getID())) {
                    tempStoic += stoichiometryMap.get(molWithH.getID());
                    stoichiometryMap.put(molWithH.getID(), tempStoic);
                } else {
                    stoichiometryMap.put(molWithH.getID(), tempStoic);
                    _metabolites.add(molWithH);
                }
            }

            try {
                Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
                sort(_metabolites, comparator);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }

            setReactantMolecule(standardizedReaction, _metabolites);
            _metabolites.clear();

            LOGGER.debug("standardize reaction module phase 2");
            LOGGER.debug("");
            LOGGER.debug("****************************");
            LOGGER.debug("");
            for (IAtomContainer mol : reaction.getProducts().atomContainers()) {
                String id = mol.getID() == null || mol.getID().trim().isEmpty() ? null : mol.getID();
                tempStoic = 1.0;
                if (reaction.getProductCoefficient(mol) > 0) {
                    tempStoic = reaction.getProductCoefficient(mol);
                }
                IAtomContainer gMol = cloneWithIDs(mol);

                /*
                 * Set old Atom Index
                 */
                for (IAtom a : gMol.atoms()) {
                    if (a.getProperties() == null) {
                        a.addProperties(new HashMap<>());
                    }
                    a.setProperty("OLD_RANK", old_atom_rank_index_product++);
                }
                fixDativeBonds(gMol);
                percieveAtomTypesAndConfigureAtoms(gMol);
                IAtomContainer molWithH = gMol;
                //= ExtAtomContainerManipulator.addExplicitH(gMol);
                aromatizeMolecule(molWithH);

                if (id == null) {
                    molWithH = setProperty(molWithH);
                } else {
                    molWithH.setID(id);
                }
                if (stoichiometryMap.containsKey(molWithH.getID())) {
                    tempStoic += stoichiometryMap.get(molWithH.getID());
                    stoichiometryMap.put(molWithH.getID(), tempStoic);
                } else {
                    stoichiometryMap.put(molWithH.getID(), tempStoic);
                    _metabolites.add(molWithH);
                }
            }

            try {
                Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
                sort(_metabolites, comparator);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }

            setProductMolecule(standardizedReaction, _metabolites);
            _metabolites.clear();
            //As per IntEnz 0 for undefined direction, 1 for LR, 2 for RL and 3 for bidirectional
            //As per CDK BIDIRECTION 1, Forward 2, Backward 0

            reactionSet.addReaction(standardizedReaction);

            //BIDIRECTION 1, Forward 2, Backward 0
            if (reaction.getDirection() != null) {
                standardizedReaction.setDirection(reaction.getDirection());
            } else {
                standardizedReaction.setDirection(BIDIRECTIONAL);
            }
            fingerprintMap.clear();
            moleculeMap.clear();
            stoichiometryMap.clear();

            LOGGER.debug("standardize reaction module end");
            return standardizedReaction;
        }

        private IAtomContainer setProperty(IAtomContainer molecule) throws Exception {
            /*
             * If ID is NULL or empty please assign it to null
             */
            String molID = molecule.getID() == null
                    || molecule.getID().isEmpty() ? null : molecule.getID();
            try {
                try {
                    if (molecule.getAtomCount() > 0) {
                        IFingerprintGenerator fpr = new FingerprintGenerator();
                        BitSet fingerprint_Present_Mol = fpr.getFingerprint(molecule);
                        /*
                        Single Atom fingerprints
                         */
                        if (fingerprint_Present_Mol.isEmpty()) {
                            long[] fp = SearchEngine.pathFingerprint(molecule, 7, 1024);
                            fingerprint_Present_Mol = com.bioinceptionlabs.reactionblast.fingerprints.ReactionFingerprinter.longArrayToBitSet(fp);
                        }
                        //Loop for Unique Mol ID Creation
                        if (!fingerprint_Present_Mol.isEmpty()) {
                            if (!isValuePresent(fingerprint_Present_Mol)) {
                                if (molID == null) {
                                    moleculeCounter += 1;
                                    int val = moleculeCounter + 100000;
                                    String Temp = Integer.toString(val);
                                    molID = Temp.replaceFirst("1", "M");
                                    molecule.setID(molID);
                                }
                                fingerprintMap.put(molID, fingerprint_Present_Mol);
                                moleculeMap.put(molID, molecule);
                            } else if (isValuePresent(fingerprint_Present_Mol)
                                    && isAtomContainerPresent(getMoleculeID(fingerprint_Present_Mol), molecule)) {
                                if (molID == null) {
                                    molID = getMoleculeID(fingerprint_Present_Mol);
                                    molecule.setID(molID);
                                }
                            } else {
                                if (molID == null) {
                                    moleculeCounter += 1;
                                    int val = moleculeCounter + 100000;
                                    String Temp = Integer.toString(val);
                                    molID = Temp.replaceFirst("1", "M");
                                    molecule.setID(molID);
                                }
                                fingerprintMap.put(molID, fingerprint_Present_Mol);
                                moleculeMap.put(molID, molecule);
                            }
                        } else {
                            LOGGER.debug("error: Fingerprint can't be generated for this molecule " + SmilesGenerator.generic().create(molecule));
                        }
                    } else {
                        LOGGER.debug("error: Mol file should contain atleast one atom! " + SmilesGenerator.generic().create(molecule));
                    }
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, " Error in setting mol id: ", ex.getMessage());
                }
                if (molecule.getID() == null) {
                    try {
                        throw new CDKException("Mol ID is NULL");
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, "Mol is can't be set ", ex.getMessage());
                    }
                }

            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            return molecule;
        }

        private void setReactantMolecule(IReaction IR, Collection<IAtomContainer> metabolites) {

            Iterator<IAtomContainer> it = metabolites.iterator();

            while (it.hasNext()) {
                IAtomContainer mol = it.next();
                mol.setProperty("STOICHIOMETRY", stoichiometryMap.get(mol.getID()));
                IR.addReactant(mol, stoichiometryMap.get(mol.getID()));
            }

            metabolites.clear();
            stoichiometryMap.clear();
        }

        private void setProductMolecule(IReaction IR, Collection<IAtomContainer> metabolites) {

            Iterator<IAtomContainer> it = metabolites.iterator();
            while (it.hasNext()) {
                IAtomContainer mol = it.next();
                mol.setProperty("STOICHIOMETRY", stoichiometryMap.get(mol.getID()));
                IR.addProduct(mol, stoichiometryMap.get(mol.getID()));
            }

            metabolites.clear();
            stoichiometryMap.clear();
        }

        /**
         *
         * @param value
         * @throws java.io.IOException
         * @return
         */
        private boolean isValuePresent(BitSet value) throws IOException, Exception {
            for (BitSet bitset : fingerprintMap.values()) {
                if (getTanimotoSimilarity(value, bitset) == 1.0) {
                    return true;
                }
            }
            return false;
        }

        /**
         *
         * @param bitset
         * @return
         * @throws java.io.IOException
         */
        private String getMoleculeID(BitSet bitset) throws IOException {
            String Key = null;
            for (Map.Entry<String, BitSet> map : fingerprintMap.entrySet()) {
                String key = map.getKey();
                try {
                    if (getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                        Key = key;
                        break;
                    }
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
            return Key;
        }

        /**
         *
         * @param key
         * @param molecule
         * @return
         * @throws Exception
         */
        private boolean isAtomContainerPresent(String key, IAtomContainer molecule) throws Exception {
            try {
                boolean flag = moleculeMap.containsKey(key);
                if (flag && molecule.getAtomCount() > 0) {
                    IAtomContainer molFromContainer = moleculeMap.get(key);
                    return isIdentical(molecule, molFromContainer, true);
                }
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            return false;
        }

        /**
         *
         * @param queryMol_org
         * @param targetMol_org
         * @param removeHydrogen
         * @return
         * @throws Exception
         */
        private boolean isIdentical(IAtomContainer queryMol_org, IAtomContainer targetMol_org, boolean removeHydrogen) throws Exception {

            IAtomContainer queryMol = queryMol_org.clone();
            IAtomContainer targetMol = targetMol_org.clone();

            if (removeHydrogen) {
                queryMol = removeHydrogens(queryMol);
                percieveAtomTypesAndConfigureAtoms(queryMol);
                aromatizeMolecule(queryMol);
                targetMol = removeHydrogens(targetMol);
                percieveAtomTypesAndConfigureAtoms(targetMol);
                aromatizeMolecule(targetMol);
            }

            if (queryMol.getAtomCount() == 1 && targetMol.getAtomCount() == 1) {
                IAtom a = queryMol.atoms().iterator().next();
                IAtom b = targetMol.atoms().iterator().next();
                return a.getSymbol().equalsIgnoreCase(b.getSymbol())
                        && Objects.equals(a.getFormalCharge(), b.getFormalCharge())
                        && queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount();
            }
            Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
            Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();

            int leftHandAtomCount = 0;

            for (IAtom a : queryMol.atoms()) {
                if (a.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                    atomUniqueCounter1.put(a.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                    atomUniqueCounter1.put(a.getSymbol(), counter);
                }
                leftHandAtomCount++;
            }
            int rightHandAtomCount = 0;

            for (IAtom b : targetMol.atoms()) {
                if (b.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                    atomUniqueCounter2.put(b.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                    atomUniqueCounter2.put(b.getSymbol(), counter);
                }
                rightHandAtomCount++;
            }
            LOGGER.debug("atomUniqueCounter1 " + leftHandAtomCount);
            LOGGER.debug("atomUniqueCounter2 " + rightHandAtomCount);

            if (leftHandAtomCount != rightHandAtomCount) {
                LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
                return false;
            } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
                LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
                return false;
            }

            return atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())
                    ? queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount()
                    ? isSubgraphIdentical(queryMol, targetMol, removeHydrogen) : false : false;
        }

        private boolean isSubgraphIdentical(IAtomContainer qMol, IAtomContainer tMol, boolean removeHydrogen) throws CDKException, IOException, CloneNotSupportedException {

            IAtomContainer mol1 = qMol.clone();
            IAtomContainer mol2 = tMol.clone();

            if (removeHydrogen) {
                mol1 = removeHydrogens(mol1);
                percieveAtomTypesAndConfigureAtoms(mol1);
                aromatizeMolecule(mol1);
                mol2 = removeHydrogens(mol2);
                percieveAtomTypesAndConfigureAtoms(mol2);
                aromatizeMolecule(mol2);
            }
            if (mol1.getAtomCount() != mol2.getAtomCount()) {
                return false;
            }
            AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(true, true);
            BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(true, true);

            Substructure mcs = new Substructure(mol1, mol2, atomMatcher, bondMatcher, false);
            mcs.setChemFilters(true, true, true);
            return mcs.isSubgraph() && !mcs.isStereoMisMatch();
        }
    }



    //~--- classes ----------------------------------------------------------------
    /**
     * @RCSfile: atomMapperTool.java,v
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
    public static class HydrogenFreeFingerPrintContainer implements Serializable {

        //define the FINGER_SIZE of the fingerprint
        //NOTE: this should be a multiple of 64 and preferably not 1024 or 2048
        //as for these values we often get the random numbers for one-atom or
        //two-atom paths the same!
        private static final int FINGER_SIZE = 64 * 30;
        //depth search is set to 6, if not given explicitly as a parameter to Fingerprinter
        private static final int DEPTH_SEARCH = 8;
        private static final long serialVersionUID = 987987606869669691L;
        private final static ILoggingTool LOGGER
                = createLoggingTool(HydrogenFreeFingerPrintContainer.class);

        /**
         *
         * @return
         */
        public static int getFingerPrintSize() {
            return FINGER_SIZE;
        }

        /**
         *
         * @return
         */
        public static int getFingerPrintDepth() {
            return DEPTH_SEARCH;
        }
        private final Map<String, BitSet> fingerPrintMap;

        //~--- constructors -------------------------------------------------------
        /**
         * HydrogenFreeFingerPrintContainer container
         */
        public HydrogenFreeFingerPrintContainer() {
            fingerPrintMap = new TreeMap<>();
        }

        public String toString() {
            return "HydrogenFreeFingerPrintContainer{" + "fingerPrintMap=" + fingerPrintMap + '}';
        }

        //~--- methods ------------------------------------------------------------
        /**
         *
         * @throws java.io.IOException
         */
        public void Clear() throws IOException {
            fingerPrintMap.clear();
        }

        /**
         *
         * @param Key
         * @throws java.io.IOException
         */
        public void Erase(String Key) throws IOException {
            fingerPrintMap.remove(Key);
        }

        /**
         *
         * @param Key
         * @param Value
         * @throws java.io.IOException
         */
        public void put(String Key, BitSet Value) throws IOException {
            try {
                fingerPrintMap.put(Key, Value);
            } catch (Exception e) {
                LOGGER.error(SEVERE, null, e);
            }
        }

        //~--- get methods --------------------------------------------------------
        /**
         *
         * @param Key
         * @throws java.io.IOException
         * @return
         */
        public BitSet getFingerPrint(String Key) throws IOException {
            return fingerPrintMap.get(Key);
        }

        /**
         *
         * @param bitset
         * @return
         * @throws java.io.IOException
         */
        public String getMoleculeID(BitSet bitset) throws IOException {
            String Key = null;
            for (Map.Entry<String, BitSet> map : fingerPrintMap.entrySet()) {
                String key = map.getKey();
                try {
                    if (getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                        Key = key;
                        break;
                    }
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }

            }
            //LOGGER.debug("Error: Unable to Find AtomContainer ID!!!");
            return Key;
        }

        /**
         *
         * @throws java.io.IOException
         * @return
         */
        public Map<String, BitSet> getFingerPrintMap()
                throws IOException {
            return new TreeMap<>(fingerPrintMap);
        }

        /**
         *
         * @param Key
         * @throws java.io.IOException
         * @return
         */
        public boolean isKeyPresent(String Key) throws IOException {
            return fingerPrintMap.containsKey(Key);
        }

        //~--- set methods --------------------------------------------------------
        /**
         *
         * @param Key
         * @param Value
         * @throws java.io.IOException
         */
        public void setValue(String Key, BitSet Value) throws IOException {
            fingerPrintMap.put(Key, Value);
        }

        /**
         *
         * @param value
         * @throws java.io.IOException
         * @return
         */
        public boolean isValuePresent(BitSet value) throws IOException {
            for (BitSet bitset : fingerPrintMap.values()) {
                try {
                    if (getTanimotoSimilarity(value, bitset) == 1.0) {
                        return true;
                    }
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }

            }
            return false;
        }

        /**
         *
         * @return
         */
        public int getSize() {
            return fingerPrintMap.size();
        }

        public boolean isEmpty() throws IOException {
            return fingerPrintMap.isEmpty();
        }

        public void write() throws IOException {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

    }



    /**
     * @RCSfile: atomMapperTool.java,v
     *
     * @Author: Syed Asad Rahman
     * @Date: 2009/06/3
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
    public static class MoleculeMoleculeMapping implements Serializable {

        private static final long serialVersionUID = 1094750239472059259L;

        //~--- fields -------------------------------------------------------------
        private final Map<String, List<MolMapping>> reactant_product_mapping_map;

        //~--- constructors -------------------------------------------------------
        /**
         *
         */
        public MoleculeMoleculeMapping() {
            reactant_product_mapping_map = new HashMap<>();
        }

        @Override
        public String toString() {
            return "MoleculeMoleculeMapping{" + "reactant_product_mapping_map=" + reactant_product_mapping_map + '}';
        }

        /**
         *
         * @throws java.io.IOException
         */
        public void Clear() throws IOException {
            reactant_product_mapping_map.clear();
        }

        /**
         *
         * @param Key
         * @throws java.io.IOException
         */
        public void Erase(String Key) throws IOException {
            reactant_product_mapping_map.remove(Key);
        }

        /**
         *
         * @param Key
         * @return
         * @throws java.io.IOException
         */
        public boolean isPresent(String Key)
                throws IOException {
            return reactant_product_mapping_map.containsKey(Key);
        }

        //~--- set methods --------------------------------------------------------
        /**
         *
         * @param Key
         * @param Value
         * @throws java.io.IOException
         */
        public void setMolMappings(String Key, List<MolMapping> Value) throws
                IOException {
            reactant_product_mapping_map.put(Key, Value);
            // Stores Reaction ID and RPAIR ID as Value in ArrayList
        }

        /**
         *
         * @param RID
         * @return
         * @throws CDKException
         */
        public List<MolMapping> getMolMappings(String RID) throws CDKException {
            return reactant_product_mapping_map.containsKey(RID) == true ? reactant_product_mapping_map.get(RID) : null;
        }

        /**
         *
         * @return Reaction count with RPAIR
         */
        public long getCount() {
            return reactant_product_mapping_map.size();
        }

        /**
         *
         * @return
         */
        public Set<String> getKeySet() {
            return reactant_product_mapping_map.keySet();
        }

        /**
         *
         * @return
         */
        public Set<Map.Entry<String, List<MolMapping>>> getEntrySet() {
            return reactant_product_mapping_map.entrySet();
        }

        /**
         *
         * @param reactionID
         * @param rName
         * @param pName
         * @return
         */
        public List<MolMapping> getMapping(String reactionID, String rName, String pName) {
            List<MolMapping> mMap = reactant_product_mapping_map.get(reactionID);
            List<MolMapping> mappedMap = new ArrayList<>();
            for (MolMapping map : mMap) {
                if ((map.getTarget().equalsIgnoreCase(rName) && map.getQuery().equalsIgnoreCase(pName))
                        || (map.getTarget().equalsIgnoreCase(pName) && map.getQuery().equalsIgnoreCase(rName))) {
                    mappedMap.add(map);
                }
            }
            return mappedMap;
        }
    }


}
