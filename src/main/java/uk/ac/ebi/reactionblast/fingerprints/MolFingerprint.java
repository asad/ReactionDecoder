/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.fingerprints;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.System.arraycopy;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import static uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator.getFingerprinterSize;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;

/**
 *
 * @author lorenzo 2007-2008
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MolFingerprint implements Comparable<MolFingerprint>,
        Comparator<MolFingerprint> {

    static final String NEW_LINE = getProperty("line.separator");
    private static final long serialVersionUID = 7057060562283378622L;

    private static synchronized MolFingerprint or(boolean[] boolArray1, boolean[] boolArray2) throws CDKException {
        if (boolArray1.length != boolArray2.length) {
            throw new CDKException("EBIFingerprint.or(boolean[], boolean[]): array with different dimensions.");
        }
        MolFingerprint res = new MolFingerprint(boolArray1);
        for (int i = 0; i < boolArray1.length; i++) {
            if (boolArray2[i] == true) {
                res.setBit(i, true);
            }
        }
        return res;
    }

    private static synchronized MolFingerprint and(boolean[] boolArray1, boolean[] boolArray2) throws CDKException {
        if (boolArray1.length != boolArray2.length) {
            throw new CDKException("EBIFingerprint.and(boolean[], boolean[]): array with different dimensions.");
        }
        MolFingerprint res = new MolFingerprint(boolArray1.length);
        for (int i = 0; i < boolArray1.length; i++) {
            if ((boolArray1[i] == true) && (boolArray2[i] == true)) {
                res.setBit(i, true);
            } else {
                res.setBit(i, false);
            }
        }
        return res;
    }
    private boolean[] arrayFingerprint = null;
    private BitSet bitsetFingerprint = null;
    private final FingerprintGenerator hashedFP = new FingerprintGenerator();

    /**
     *
     */
    public MolFingerprint() {
        this.arrayFingerprint = new boolean[0];
    }

    /**
     * Prepare the target molecule for analysis.
     * <p/>
     * We perform ring perception and aromaticity detection and set up the
     * appropriate properties. Right now, this function is called each time we
     * need to do a query and this is inefficient.
     *
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    /**
     *
     * @param mol
     * @throws CDKException
     */
    public MolFingerprint(IAtomContainer mol) throws CDKException {
        this();
        try {
            this.bitsetFingerprint = hashedFP.getFingerprint(mol);
            this.set(this.bitsetFingerprint);
            arrayFingerprint = new boolean[getFingerprinterSize()];
            for (int i = 0; i < getFingerprinterSize(); i++) {
                arrayFingerprint[i] = (this.bitsetFingerprint.get(i));
            }
        } catch (CDKException e) {
            throw new CDKException("Failed to create CDKMolecularDescriptor "
                    + "while constructing EBIFingerprint " + mol.getAtomCount() + "," + NEW_LINE + e.getMessage());
        }
    }

    /**
     *
     * @param fgrprt
     */
    public MolFingerprint(BitSet fgrprt) {
        this();
        arrayFingerprint = new boolean[fgrprt.size()];
        for (int i = 0; i < fgrprt.length(); i++) {
            arrayFingerprint[i] = (fgrprt.get(i));
        }
    }

    /**
     *
     * @param length
     */
    public MolFingerprint(int length) {
        this();
        arrayFingerprint = new boolean[length];
        set(false);
    }

    /**
     *
     * @param fgprt
     * @throws CDKException
     */
    public MolFingerprint(boolean[] fgprt) throws CDKException {
        this();
        arrayFingerprint = new boolean[fgprt.length];
        arraycopy(fgprt, 0, arrayFingerprint, 0, fgprt.length);
    }

    /**
     *
     * @param molFingerprint
     * @throws CDKException
     */
    public MolFingerprint(MolFingerprint molFingerprint) throws CDKException {
        this();
        arrayFingerprint = new boolean[molFingerprint.getBooleanArray().length];
        arraycopy(molFingerprint.getBooleanArray(), 0, this.arrayFingerprint, 0, arrayFingerprint.length);
    }

    private synchronized void set(boolean value) {
        for (int i = 0; i < arrayFingerprint.length; i++) {
            arrayFingerprint[i] = value;
        }
    }

    private synchronized void set(BitSet bitset) {
        arrayFingerprint = new boolean[bitset.size()];
        for (int i = 0; i < bitset.length(); i++) {
            arrayFingerprint[i] = (bitset.get(i));
        }
    }

    /**
     * Returns binary arrayFingerprint as bitset
     *
     * @return
     */
    public synchronized BitSet getBitSet() {
        BitSet bts = new BitSet(arrayFingerprint.length);
        for (int i = 0; i < arrayFingerprint.length; i++) {
            bts.set(i, arrayFingerprint[i]);
        }
        return bts;
    }

    /**
     *
     * @param fromIndex
     * @param molFingerprint
     * @throws CDKException
     */
    private synchronized void set(int fromIndex, MolFingerprint molFingerprint) throws CDKException {
        for (int i = fromIndex; (i < fromIndex + molFingerprint.length()) && (i < arrayFingerprint.length); i++) {
            arrayFingerprint[i] = molFingerprint.getBooleanArray()[i - fromIndex];
        }
    }

    private synchronized void set(int fromIndex, boolean[] fgprt) throws CDKException {
        for (int i = fromIndex; (i < fromIndex + fgprt.length) && (i < arrayFingerprint.length); i++) {
            arrayFingerprint[i] = fgprt[i - fromIndex];
        }
    }

    @Override
    public synchronized String toString() {
        String strFp = "";
        for (int i = 0; i < arrayFingerprint.length; i++) {
            strFp += (arrayFingerprint[i] ? "1" : "0");
        }
        return strFp;
    }

    /**
     *
     */
    public synchronized void println() {
        out.println(toString());
    }

    /**
     *
     * @return
     */
    public synchronized int length() {
        return arrayFingerprint.length;
    }

    /**
     *
     * @param index
     * @return
     * @throws CDKException
     */
    public synchronized boolean getBit(int index) throws CDKException {
        if ((index >= arrayFingerprint.length) || (index < 0)) {
            throw new CDKException("EBIFingerprint.getBit(int index) failed for index out of bounds.");
        }
        return arrayFingerprint[index];
    }

    /**
     *
     * @param index
     * @param value
     * @throws CDKException
     */
    public synchronized void setBit(int index, boolean value) throws CDKException {
        if ((index >= arrayFingerprint.length) || (index < 0)) {
            throw new CDKException("EBIFingerprint.setBit(int index, boolean value) failed for index out of bounds.");
        }
        arrayFingerprint[index] = value;
    }

    /**
     *
     * @return
     */
    public synchronized boolean[] getBooleanArray() {
        boolean[] bs = new boolean[arrayFingerprint.length];
        arraycopy(arrayFingerprint, 0, bs, 0, arrayFingerprint.length);
        return bs;
    }

    /**
     *
     * @param b
     * @throws CDKException
     */
    public synchronized void append(Byte b) throws CDKException {
        boolean[] bt = new boolean[8];
        for (int i = 0; i < 8; i++) {
            bt[i] = (b & (1 << (7 - i))) != 0;
        }
        append(new MolFingerprint(bt));
    }

    /**
     *
     * @param fp
     * @throws CDKException
     */
    public synchronized void append(MolFingerprint fp) throws CDKException {
        MolFingerprint newFp = new MolFingerprint(arrayFingerprint.length + fp.length());
        newFp.set(0, arrayFingerprint);
        newFp.set(arrayFingerprint.length, fp);
        arrayFingerprint = newFp.getBooleanArray();
    }

    /**
     *
     * @param molFp
     * @return
     * @throws CDKException
     */
    public synchronized MolFingerprint or(MolFingerprint molFp) throws CDKException {
        return or(arrayFingerprint, molFp.getBooleanArray());
    }

    /**
     *
     * @param molFp
     * @return
     * @throws CDKException
     */
    public synchronized MolFingerprint and(MolFingerprint molFp) throws CDKException {
        return and(arrayFingerprint, molFp.getBooleanArray());
    }

    /**
     *
     * @param fingerprint
     * @return
     * @throws Exception
     */
    public synchronized double similarity(MolFingerprint fingerprint) throws Exception {
        double similarity;
        similarity = getTanimotoSimilarity(fingerprint.getBitSet(), bitsetFingerprint);
        return similarity;
    }

    /**
     * Returns 0 if two fingerprints are equal and if they share same labels it
     * returns difference in their weight
     *
     * @param t
     * @return
     */
    @Override
    public synchronized int compareTo(MolFingerprint t) {
        return compare(this, t);
    }

    /**
     * Return true if two Fingerprints are equal
     *
     * @param object
     * @return
     */
    @Override
    public synchronized boolean equals(Object object) {
        if (this == object) {
            return true;
        }
        if (!(object instanceof MolFingerprint)) {
            return false;
        }

        MolFingerprint fpn = (MolFingerprint) object;

        if (this.arrayFingerprint.length != fpn.getBooleanArray().length) {
            return false;
        }

        for (int i = 0; i < arrayFingerprint.length; i++) {
            if (this.arrayFingerprint[i] != arrayFingerprint[i]) {
                return false;
            }
        }
        return true;
    }

    @Override
    public synchronized int hashCode() {
        int hash = 7;
        hash = 19 * hash + Arrays.hashCode(this.arrayFingerprint);
        hash = 19 * hash + (this.bitsetFingerprint != null ? this.bitsetFingerprint.hashCode() : 0);
        return hash;
    }

    /**
     * Returns 0 if two fingerprints are equal and if they share same labels it
     * returns difference in their weight
     *
     * @param o1
     * @param o2
     * @return
     */
    @Override
    public int compare(MolFingerprint o1, MolFingerprint o2) {
        int len1 = o1.getBooleanArray().length;
        int len2 = o2.getBooleanArray().length;
        int n = min(len1, len2);
        if (len1 == len2) {
            if (o1.equals(o2)) {
                return 0;
            } else {
                return -1;
            }
        }
        return max(len1, len2) - n;
    }

}
