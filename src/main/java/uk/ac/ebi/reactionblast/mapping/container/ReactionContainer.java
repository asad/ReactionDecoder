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
package uk.ac.ebi.reactionblast.mapping.container;

import java.io.IOException;
import java.io.Serializable;
import java.util.BitSet;
import static java.util.Collections.synchronizedMap;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import static uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator.getFingerprinterSize;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;

/**
 *
 *
 * @Author: Syed Asad Rahman
 * @Contact: asad@ebi.ac.uk
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
        eAtomContainerMap = synchronizedMap(new TreeMap<>());
        pAtomContainerMap = synchronizedMap(new TreeMap<>());
        eFingerPrintMap = synchronizedMap(new TreeMap<>());
        pFingerPrintMap = synchronizedMap(new TreeMap<>());
        eductContainerModificationMap = synchronizedMap(new TreeMap<>());
        productContainerModificationMap = synchronizedMap(new TreeMap<>());
        fpr = new FingerprintGenerator();
    }

    @Override
    public String toString() {
        return "ReactionContainer{" + "eAtomContainerMap=" + eAtomContainerMap + ", pAtomContainerMap=" + pAtomContainerMap + ", eFingerPrintMap=" + eFingerPrintMap + ", pFingerPrintMap=" + pFingerPrintMap + ", eductContainerModificationMap=" + eductContainerModificationMap + ", productContainerModificationMap=" + productContainerModificationMap + ", fpr=" + fpr + '}';
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    public synchronized void Clear() throws IOException {
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
    public synchronized void eraseEduct(int Key) throws IOException {
        eAtomContainerMap.remove(Key);
        eFingerPrintMap.remove(Key);
        eductContainerModificationMap.remove(Key);
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    public synchronized void eraseProduct(int Key) throws IOException {
        pAtomContainerMap.remove(Key);
        pFingerPrintMap.remove(Key);
        productContainerModificationMap.remove(Key);
    }

    /**
     *
     * @return
     */
    public synchronized Integer getEductCount() {
        return eAtomContainerMap.size();
    }

    /**
     *
     * @return
     */
    public synchronized Integer getProductCount() {
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
    public synchronized IAtomContainer getEduct(int value)
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
    public synchronized IAtomContainer getProduct(int value)
            throws IOException, CDKException {
        return pAtomContainerMap.containsKey(value) ? pAtomContainerMap.get(value) : null;
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param value
     * @return
     * @throws java.io.IOException
     * @throws CDKException
     */
    public synchronized boolean isEductModified(int value)
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
    public synchronized boolean isProductModified(int value)
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
    public synchronized void putEduct(int index, IAtomContainer educt)
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
    public synchronized void putProduct(int index, IAtomContainer product)
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
    public synchronized void setEductModified(int index, boolean flag)
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
    public synchronized void setProductModified(int index, boolean flag)
            throws IOException, Exception {
        productContainerModificationMap.put(index, flag);
    }

    /**
     *
     * @param eudMap
     * @throws java.io.IOException
     * @throws Exception
     */
    public synchronized void putAllEduct(TreeMap<Integer, IAtomContainer> eudMap)
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
    public synchronized void putAllProduct(TreeMap<Integer, IAtomContainer> prodMap)
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
    public synchronized BitSet getFingerPrintofEduct(int value)
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
    public synchronized BitSet getFingerPrintofProduct(int value)
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
    private synchronized void setFingerPrintofEduct(int value, BitSet edu)
            throws IOException {
        eFingerPrintMap.put(value, edu);

    }

    /**
     *
     * @param value
     * @param Prod
     * @throws java.io.IOException
     */
    private synchronized void setFingerPrintofProduct(int value, BitSet prod) throws IOException {
        pFingerPrintMap.put(value, prod);

    }

    @Override
    public synchronized Object clone() throws CloneNotSupportedException {
        return super.clone();
    }
}
