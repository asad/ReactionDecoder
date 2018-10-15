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
package uk.ac.ebi.reactionblast.mechanism;

import java.io.Serializable;
import static java.lang.Math.log10;
import java.util.HashMap;
import java.util.Map.Entry;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;

/**
 * This class maintains the references of a set of RMatrices. The method
 * getGlobalEntropy can be called in order to get the global entropy of a stored
 * RMatrix.
 *
 * @author Lorenzo Baldacci {lorenzo@ebi.ac.uk|lbaldacc@csr.unibo.it}
 */
public class GlobalEntropyCalculator implements Serializable {

    private static final long serialVersionUID = 7879978965972591251L;
    private final HashMap<Integer, RMatrix> matrixMap = new HashMap<>();
    private final HashMap<String, HashMap<Integer, Integer>> typeMap = new HashMap<>();
    private final HashMap<String, Integer> freqMap = new HashMap<>();

    /**
     * Class constructor.
     */
    public GlobalEntropyCalculator() {
    }

    /**
     * The method adds a RMatrix to the ArrayList maintaining the set of
     * RMatrices.
     *
     * @param index The unique index associated to the RMatrix.
     * @param m The RMatrix to be added.
     * @throws CDKException
     */
    public void addMatrix(int index, RMatrix m) throws CDKException {
        matrixMap.put(index, m);
        String key = "";
        HashMap<Integer, Integer> val = null;
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = i + 1; j < m.getColumnDimension(); j++) {
                IAtom ai = m.getReactantAtom(i);
                IAtom aj = m.getReactantAtom(j);
                if (ai.getSymbol().compareTo(aj.getSymbol()) <= 0) {
                    key = ai.getSymbol() + "." + aj.getSymbol();
                } else {
                    key = aj.getSymbol() + "." + ai.getSymbol();
                }
                val = typeMap.get(key);
                if (val == null) {
                    val = new HashMap<>();
                    typeMap.put(key, val);
                    freqMap.put(key, 0);
                }
                Integer frequency = val.get(new Double(m.getValue(i, j)).intValue());
                if (frequency == null) {
                    frequency = 0;
                }
                val.put(new Double(m.getValue(i, j)).intValue(), frequency + 1);
                freqMap.put(key, freqMap.get(key) + 1);
            }
        }
    }

    /**
     * The methods returns the global entropy of the RMatrix in idx-th position.
     *
     * @param idx The RMatrix position.
     * @return The global entropy of the RMatrix in idx-th position
     * @throws CDKException
     */
    public synchronized double getGlobalEntropy(int idx) throws CDKException {
        RMatrix m = matrixMap.get(idx);
        double lH = 0;
        String key = "";
        HashMap<Integer, Integer> val = null;
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = i + 1; j < m.getColumnDimension(); j++) {
                IAtom ai = m.getReactantAtom(i);
                IAtom aj = m.getReactantAtom(j);
                if (ai.getSymbol().compareTo(aj.getSymbol()) <= 0) {
                    key = ai.getSymbol() + "." + aj.getSymbol();
                } else {
                    key = aj.getSymbol() + "." + ai.getSymbol();
                }
                int matValue = (int) m.getValue(i, j);
                val = typeMap.get(key);
                int freq = val.get(matValue);
                double px = freq / freqMap.get(key);
                if (px != 0.0f) {
                    lH += px * log10(px) / log10(2);
                }
            }
        }
        return -lH / ((m.getRowDimension() * m.getColumnDimension() - m.getRowDimension()) / 2);
    }

    /**
     *
     * @return
     */
    public synchronized Iterable<Integer> getKeySet() {
        return matrixMap.keySet();
    }

    /**
     * Returns the index of the RMatrix.
     *
     * @param m The RMatrix for which the index is required.
     * @return The index of the RMatrix, -1 if the RMatrix is not found.
     * @throws CDKException
     */
    public synchronized int getMatrixIndex(RMatrix m) throws CDKException {
        int index = -1;
        if (!matrixMap.containsValue(m)) {
            throw new CDKException("GlobalEntropyCalculator.getMatrixIndex: The matrix has not been found.");
        }
        for (Entry<Integer, RMatrix> e : matrixMap.entrySet()) {
            if (e.getValue() == m) {
                index = e.getKey();
            }
        }
        return index;
    }

    /**
     * Returns the number of RMatrices stored in the object.
     *
     * @return The number of RMatrices stored in the object.
     */
    public synchronized int getMatrixCount() {
        return matrixMap.size();
    }

    /**
     * Returns the RMatrix in idx-th position
     *
     * @param index The index of the required RMatrix
     * @return The RMatrix in idx-th position
     * @throws CDKException
     */
    public synchronized RMatrix getMatrix(int index) throws CDKException {
        if (!matrixMap.containsKey(index)) {
            throw new CDKException("The index of the required RMatrix is out of bouns");
        }
        return matrixMap.get(index);
    }
}
