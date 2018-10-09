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

package uk.ac.ebi.reactionblast.fingerprints.interfaces;

import java.util.BitSet;
import java.util.Collection;
import org.openscience.cdk.exception.CDKException;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public interface IPatternFingerprinter extends Comparable<IPatternFingerprinter> {

    /**
     *
     * @return
     */
    public abstract double[] getValuesAsArray();

    /**
     *
     * @return Keys of the fingerprints
     */
    public abstract Collection<IFeature> getFeatures();

    /**
     *
     * @return values of the fingerprints
     */
    public abstract Collection<Double> getValues();

    /**
     *
     * @param index
     * @return
     * @throws CDKException
     */
    public abstract IFeature getFeature(int index) throws CDKException;

    /**
     *
     * @param pattern
     * @return
     */
    public abstract Double getWeight(String pattern);

    /**
     *
     * @param index
     * @return
     */
    public abstract Double getWeight(int index);

    /**
     * Size of the hashed fingerprint
     *
     * @return
     */
    public abstract int getFingerprintSize();

    /**
     *
     * @return
     */
    @Override
    public abstract int hashCode();

    /**
     *
     * @param feature
     * @return
     */
    @Override
    public abstract boolean equals(Object feature);

    /**
     *
     * @return
     */
    @Override
    public abstract String toString();

    /**
     *
     * @return
     */
    public abstract BitSet getHashedFingerPrint();

    /**
     *
     * @param fingerprint
     * @throws CDKException
     */
    public abstract void addBinary(BitSet fingerprint) throws CDKException;

    /**
     *
     * @param fngp
     * @throws CDKException
     */
    public abstract void add(IPatternFingerprinter fngp) throws CDKException;

    /**
     * @return the fingerprintID
     */
    public abstract String getFingerprintID();

    /**
     * @param fingerprintID the fingerprintID to set
     */
    public abstract void setFingerprintID(String fingerprintID);

    /**
     * Number of unique features of this fingerprint
     *
     * @return
     */
    public abstract int getFeatureCount();

    /**
     *
     * @param key
     * @return
     */
    public abstract boolean hasFeature(IFeature key);

    /**
     *
     * @return
     */
    public abstract double[] getWeightedHashedFingerPrint();

    /**
     *
     * @param feature
     * @throws CDKException
     */
    public abstract void add(IFeature feature) throws CDKException;

    /**
     *
     * @return
     * @throws CloneNotSupportedException
     */
    public IPatternFingerprinter clone() throws CloneNotSupportedException;
}
