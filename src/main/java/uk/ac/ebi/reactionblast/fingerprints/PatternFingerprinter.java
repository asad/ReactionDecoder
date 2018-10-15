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

import java.io.Serializable;
import static java.lang.String.valueOf;
import static java.lang.System.getProperty;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import static java.util.Collections.synchronizedSortedSet;
import static java.util.Collections.unmodifiableCollection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator.getFingerprinterSize;
import static uk.ac.ebi.reactionblast.fingerprints.PatternComparators.overallComparator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class PatternFingerprinter implements Cloneable, IPatternFingerprinter,
        Comparable<IPatternFingerprinter>,
        Comparator<IPatternFingerprinter>,
        Serializable {

    private static final long serialVersionUID = 0156306561546552043757L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(PatternFingerprinter.class);

    /**
     *
     * @param map
     * @return
     */
    public static IPatternFingerprinter makePatternFingerprint(Map<String, Double> map) {
        return makePatternFingerprint(map.keySet(), map.values());
    }

    /**
     *
     * @param keyCollection
     * @param valueCollection
     * @return
     */
    public static IPatternFingerprinter makePatternFingerprint(Collection<String> keyCollection, Collection<Double> valueCollection) {
        List<IFeature> features = new ArrayList<>();

        List<String> keyList = new ArrayList<>(keyCollection);
        List<Double> valueList = new ArrayList<>(valueCollection);
        for (int index = 0; index < keyList.size(); index++) {
            String key = keyList.get(index);
            Double value = valueList.get(index);
            features.add(new Feature(key, value));
        }
        return new PatternFingerprinter(features);
    }
    private final Set<IFeature> featureSet;
    private String fingerprintID = "?";
    private int fingerprintSize;

    /**
     *
     */
    public PatternFingerprinter() {
        this(getFingerprinterSize());
    }

    /**
     *
     * @param features
     */
    public PatternFingerprinter(
            Collection<IFeature> features) {
        this(features, getFingerprinterSize());
    }

    /**
     *
     * @param fingerprintSize
     */
    public PatternFingerprinter(int fingerprintSize) {
        this.fingerprintSize = fingerprintSize;
        featureSet = synchronizedSortedSet(new TreeSet<IFeature>());
    }

    /**
     *
     * @param features
     * @param fingerprintSize
     */
    public PatternFingerprinter(Collection<IFeature> features, int fingerprintSize) {
        this(fingerprintSize);
        for (final IFeature feature : features) {
            if (!this.featureSet.contains(feature)) {
                this.featureSet.add(new Feature(feature.getPattern()));
            } else {
                for (IFeature localFeature : featureSet) {
                    if (localFeature.getPattern().equals(feature.getPattern())) {
                        double newWeight = localFeature.getWeight() + feature.getWeight();
                        localFeature.setValue(newWeight);
                        break;
                    }
                }
            }
        }
    }

    /**
     *
     * @param fingerprint
     * @throws CDKException
     */
    @Override
    public synchronized void addBinary(BitSet fingerprint) throws CDKException {
        if (featureSet == null) {
            throw new CDKException("Cannot perform PatternFingerprint.add() as Fingerprint not initialized");
        }
        for (int i = 0; i < fingerprint.size(); i++) {
            if (fingerprint.get(i)) {
                add(new Feature(valueOf(i), 1.0));
            }
        }
    }

    /**
     *
     * @throws CDKException
     */
    @Override
    public synchronized void add(IFeature feature) throws CDKException {
        if (featureSet == null) {
            throw new CDKException("Cannot perform PatternFingerprint.add() as Fingerprint not initialized");
        }

        if (!this.featureSet.contains(feature)) {
            this.featureSet.add(new Feature(feature.getPattern(), feature.getWeight()));
        } else {
            for (IFeature localFeature : featureSet) {
                if (localFeature.getPattern().equals(feature.getPattern())) {
                    double newWeight = localFeature.getWeight() + feature.getWeight();
                    localFeature.setValue(newWeight);
                    break;
                }
            }
        }
    }

    /**
     *
     * @param fngp
     * @throws CDKException
     */
    @Override
    public synchronized void add(IPatternFingerprinter fngp) throws CDKException {
        if (featureSet == null || fngp == null) {
            throw new CDKException("Cannot perform PatternFingerprint.add() as Fingerprint not initialized");
        }
        if (fngp.getFingerprintSize() != this.fingerprintSize) {
            throw new CDKException("Cannot perform PatternFingerprint.add() as Fingerprint size not equal");
        }
        for (IFeature feature : fngp.getFeatures()) {
            if (!this.featureSet.contains(feature)) {
                this.featureSet.add(new Feature(feature.getPattern(), feature.getWeight()));
            } else {
                for (IFeature localFeature : featureSet) {
                    if (localFeature.getPattern().equals(feature.getPattern())) {
                        double newWeight = localFeature.getWeight() + feature.getWeight();
                        localFeature.setValue(newWeight);
                        break;
                    }
                }
            }
        }
    }

    @Override
    public double[] getValuesAsArray() {
        int pos = 0;
        double[] res = new double[featureSet.size()];
        for (IFeature feature : featureSet) {
            res[pos] = feature.getWeight();
            pos += 1;
        }
        return res;
    }

    @Override
    public Collection<IFeature> getFeatures() {
        return unmodifiableCollection(featureSet);
    }

    @Override
    public Collection<Double> getValues() {
        List<Double> collection = new ArrayList<>(featureSet.size());
        int i = 0;
        for (IFeature feature : featureSet) {
            collection.add(i, feature.getWeight());
            i += 1;
        }
        return collection;
    }

    @Override
    public int getFeatureCount() {
        return featureSet.size();
    }

    @Override
    public BitSet getHashedFingerPrint() {
        double[] weightedHashedFingerPrint = getWeightedHashedFingerPrint();
        BitSet binary = new BitSet(this.fingerprintSize);
        for (int i = 0; i < weightedHashedFingerPrint.length; i++) {
            if (weightedHashedFingerPrint[i] > 0.) {
                binary.set(i, true);
            } else {
                binary.set(i, false);
            }
        }
        return binary;
    }

    /**
     *
     * @return
     */
    @Override
    public double[] getWeightedHashedFingerPrint() {
        RandomNumber randomNumberGen = new RandomNumber();
        double[] hashedFingerPrint = new double[this.fingerprintSize];
        for (int i = 0; i < hashedFingerPrint.length; i++) {
            hashedFingerPrint[i] = 0.;
        }
        Collection<IFeature> features = this.getFeatures();
        features.stream().forEach((feature) -> {
            long hashCode = feature.hashCode();
            int randomNumber = randomNumberGen.generateMersenneTwisterRandomNumber(this.fingerprintSize, hashCode);
            hashedFingerPrint[randomNumber] += feature.getWeight();
        });
        return hashedFingerPrint;
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");
        DecimalFormat df = new DecimalFormat();
        result.append(NEW_LINE);
        result.append("ID=").append(this.fingerprintID);
        result.append(" (").append(this.featureSet.size()).append(")");
        result.append(NEW_LINE);
        result.append(this.getFeatures());
        result.append(NEW_LINE);
        return result.toString();
    }

    @Override
    public IFeature getFeature(int index) throws CDKException {
        if (featureSet.size() >= index) {
            int i = 0;
            for (final IFeature key : featureSet) {
                if (i == index) {
                    return key;
                }
                i++;
            }
        }
        return null;
    }

    @Override
    public Double getWeight(String pattern) {
        if (!featureSet.isEmpty()) {
            int i = 0;
            for (final IFeature key : featureSet) {
                if (key.getPattern().equals(pattern)) {
                    return key.getWeight();
                }
                i++;
            }
        }
        return -1.0;
    }

    @Override
    public Double getWeight(int index) {
        if (featureSet.size() >= index) {
            int i = 0;
            for (final IFeature value : featureSet) {
                if (i == index) {
                    return value.getWeight();
                }
                i++;
            }
        }
        return -1.0;
    }

    /**
     * @return the fingerprintID
     */
    @Override
    public String getFingerprintID() {
        return fingerprintID;
    }

    /**
     * @param fingerprintID the fingerprintID to set
     */
    @Override
    public void setFingerprintID(String fingerprintID) {
        this.fingerprintID = fingerprintID;
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
    public synchronized int compare(IPatternFingerprinter o1, IPatternFingerprinter o2) {
        Comparator<IPatternFingerprinter> comparator = overallComparator();
        return comparator.compare(o1, o2);
    }

    /**
     * Returns 0 if two fingerprints are equal and if they share same labels it
     * returns difference in their weight
     *
     * @param t
     * @return
     */
    @Override
    public synchronized int compareTo(IPatternFingerprinter t) {
        return compare(this, t);
    }

    /**
     * Return true if two Fingerprints are equal
     *
     * @param object
     * @return
     */
    @Override
    public boolean equals(Object object) {
        if (object == null) {
            return false;
        }
        if (getClass() != object.getClass()) {
            return false;
        }
        final PatternFingerprinter other = (PatternFingerprinter) object;
        if (this.featureSet != other.featureSet && (this.featureSet == null || !this.featureSet.equals(other.featureSet))) {
            return false;
        }
        if ((this.fingerprintID == null) ? (other.fingerprintID != null) : !this.fingerprintID.equals(other.fingerprintID)) {
            return false;
        }
        return this.fingerprintSize == other.fingerprintSize;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 83 * hash + (this.featureSet != null ? this.featureSet.hashCode() : 0);
        hash = 83 * hash + (this.fingerprintID != null ? this.fingerprintID.hashCode() : 0);
        hash = 83 * hash + this.fingerprintSize;
        return hash;
    }

    @Override
    public int getFingerprintSize() {
        return fingerprintSize;
    }

    /**
     *
     * @param key
     * @return
     */
    @Override
    public boolean hasFeature(IFeature key) {
        return this.featureSet.contains(key);
    }

    @Override
    public IPatternFingerprinter clone() throws CloneNotSupportedException {
        IPatternFingerprinter p = new PatternFingerprinter(this.fingerprintSize);
        try {
            p.add(this);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return p;
    }

}
