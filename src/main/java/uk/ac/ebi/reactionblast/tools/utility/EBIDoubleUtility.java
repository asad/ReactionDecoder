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
package uk.ac.ebi.reactionblast.tools.utility;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class EBIDoubleUtility {

    private static final long serialVersionUID = 7683452581122892189L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(EBIDoubleUtility.class);

    /**
     *
     * @param val1
     * @param val2
     * @return fused array val1+val2
     * @throws CDKException
     */
    public static double[] append(double[] val1, double[] val2) throws CDKException {

        double[] feature = null;

        if (val1.length > 0 && val2.length > 0) {
            feature = new double[val1.length + val2.length];

            int index = 0;
            for (int i = 0; i < val1.length; i++) {
                feature[index++] = val1[i];
            }

            for (int j = 0; j < val2.length; j++) {
                feature[index++] = val2[j];
            }

        } else {
            throw new CDKException("Index < 0: ");
        }

        return feature;

    }

    /**
     *
     * @param val1
     * @param val2
     * @return fused array val1+val2
     * @throws CDKException
     */
    public static IPatternFingerprinter Union(IPatternFingerprinter val1, IPatternFingerprinter val2) throws CDKException {
        PatternFingerprinter patternFingerprinter = new PatternFingerprinter(val1.getFingerprintSize() + val2.getFingerprintSize());

        if (val1.getFingerprintSize() > 0 && val2.getFingerprintSize() > 0) {

            for (int i = 0; i < val1.getFeatureCount(); i++) {
                IFeature feature = val1.getFeature(i);
                patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
            }

            for (int j = 0; j < val2.getFeatureCount(); j++) {
                IFeature feature = val2.getFeature(j);
                patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));

            }

        } else {
            throw new CDKException("Index < 0: ");
        }
        return patternFingerprinter;
    }

    /**
     *
     * @param val1
     * @param val2
     * @param val3
     * @return fused array val1+val2+val3
     * @throws CDKException
     */
    public static double[] append(double[] val1, double[] val2, double[] val3) throws CDKException {

        double[] feature = null;

        if (val1.length > 0 && val2.length > 0 && val3.length > 0) {
            feature = new double[val1.length + val2.length + val3.length];

            int index = 0;
            for (int i = 0; i < val1.length; i++) {
                feature[index++] = val1[i];
            }

            for (int j = 0; j < val2.length; j++) {
                feature[index++] = val2[j];
            }

            for (int k = 0; k < val3.length; k++) {
                feature[index++] = val3[k];
            }

        } else {
            throw new CDKException("Index < 0: ");
        }

        return feature;

    }

    /**
     *
     * @param val1
     * @param val2
     * @param val3
     * @return fused array val1+val2+val3
     * @throws CDKException
     */
    public static IPatternFingerprinter Union(IPatternFingerprinter val1, IPatternFingerprinter val2, IPatternFingerprinter val3) throws CDKException {

        PatternFingerprinter patternFingerprinter = new PatternFingerprinter(val1.getFingerprintSize() + val2.getFingerprintSize() + val3.getFingerprintSize());

        if (val1.getFingerprintSize() > 0 && val2.getFingerprintSize() > 0) {

            for (int i = 0; i < val1.getFeatureCount(); i++) {
                IFeature feature = val1.getFeature(i);
                patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
            }

            for (int j = 0; j < val2.getFeatureCount(); j++) {
                IFeature feature = val2.getFeature(j);
                patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
            }

            for (int k = 0; k < val3.getFeatureCount(); k++) {
                IFeature feature = val3.getFeature(k);
                patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
            }

        } else {
            throw new CDKException("Index < 0: ");
        }
        return patternFingerprinter;

    }

    /**
     *
     * @param val1
     * @param val2
     * @return is val1 contained in val2
     * @throws CDKException
     */
    public static boolean isSubset(double[] val1, double[] val2) throws CDKException {

        boolean flag = true;

        if (val1.length > 0 && val2.length > 0) {

            for (int i = 0; i < val1.length; i++) {
                if (val1[i] > val2[i]) {
                    flag = false;
                    break;
                }
            }

        } else {
            throw new CDKException("Index <0: ");
        }

        return flag;

    }

    /**
     *
     * @param val1
     * @param val2
     * @return is val2 contained in val1
     * @throws CDKException
     */
    public static boolean isSuperset(double[] val1, double[] val2) throws CDKException {

        boolean flag = true;

        if (val1.length > 0 && val2.length > 0) {

            for (int i = 0; i < val1.length; i++) {
                if (val1[i] < val2[i]) {
                    flag = false;
                    break;
                }
            }

        } else {
            throw new CDKException("Index <0: ");
        }

        return flag;
    }

    private EBIDoubleUtility() {
    }
}
