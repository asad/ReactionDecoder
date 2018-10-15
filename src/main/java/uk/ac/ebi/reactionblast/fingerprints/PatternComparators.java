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
import java.util.Comparator;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class PatternComparators {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(PatternComparators.class);

    public static Comparator<IPatternFingerprinter> overallComparator() {
        return (IPatternFingerprinter o1, IPatternFingerprinter o2) -> {
            int len1 = o1.getFeatureCount();
            int len2 = o2.getFeatureCount();
            if (!o1.getFingerprintID().equals(o2.getFingerprintID())) {
                return o1.getFingerprintID().compareTo(o2.getFingerprintID());
            }
            int n = min(len1, len2);
            if (len1 == len2) {
                int pos = 0;
                while (n-- != 0) {
                    try {
                        if (!o1.getFeature(pos).equals(o2.getFeature(pos))) {
                            return o1.getFeature(pos).compareTo(o2.getFeature(pos));
                        } else if (!o1.getFeature(pos).equals(o2.getFeature(pos))) {
                            double v1 = o1.getWeight(pos);
                            double v2 = o2.getWeight(pos);
                            if (v1 != v2) {
                                return (int) (max(v1, v2) - min(v1, v2));
                            }
                        }
                    } catch (CDKException ex) {
                        LOGGER.error(SEVERE, null, ex);
                    }
                    pos++;
                }
            }
            return max(len1, len2) - n;
        };
    }

    public static Comparator<IPatternFingerprinter> dataComparator() {
        return (IPatternFingerprinter o1, IPatternFingerprinter o2) -> {
            int len1 = o1.getFeatureCount();
            int len2 = o2.getFeatureCount();
            if (!o1.getFingerprintID().equals(o2.getFingerprintID())) {
                return o1.getFingerprintID().compareTo(o2.getFingerprintID());
            }
            int n = min(len1, len2);
            if (len1 == len2) {
                int pos = 0;
                while (n-- != 0) {
                    double v1 = o1.getWeight(pos);
                    double v2 = o2.getWeight(pos);
                    if (v1 != v2) {
                        return (int) (max(v1, v2) - min(v1, v2));
                    }
                    pos++;
                }
            }
            return max(len1, len2) - n;
        };
    }

    private PatternComparators() {
    }
}
