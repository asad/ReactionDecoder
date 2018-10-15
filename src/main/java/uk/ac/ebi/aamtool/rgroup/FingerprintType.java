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
package uk.ac.ebi.aamtool.rgroup;

import java.util.Set;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FingerprintType {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(FingerprintType.class);

    private final Set<String> commonCommonFP;
    private final Set<String> commonDifferenceFP;
    private final Set<String> allPatternsFP;

    /**
     *
     * @param commonCommonFP
     * @param commonDifferenceFP
     * @param allPatternsFP
     */
    public FingerprintType(Set<String> commonCommonFP, Set<String> commonDifferenceFP, Set<String> allPatternsFP) {
        this.commonCommonFP = commonCommonFP;
        this.commonDifferenceFP = commonDifferenceFP;
        this.allPatternsFP = allPatternsFP;
    }

    /**
     * @return the commonCommonFP
     */
    public Set<String> getCommonCommonFP() {
        return commonCommonFP;
    }

    /**
     * @return the commonDifferenceFP
     */
    public Set<String> getCommonDifferenceFP() {
        return commonDifferenceFP;
    }

    /**
     * @return the allPatternsFP
     */
    public Set<String> getAllPatternsFP() {
        return allPatternsFP;
    }

}
