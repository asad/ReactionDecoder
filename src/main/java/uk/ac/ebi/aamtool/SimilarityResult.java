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
package uk.ac.ebi.aamtool;

import java.util.Map;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class SimilarityResult {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(SimilarityResult.class);

    private final String query;
    private final String target;
    private final Map<String, String> similarityReactions;

    /**
     *
     * @param query
     * @param target
     * @param similarityReactions
     */
    public SimilarityResult(String query, String target, Map<String, String> similarityReactions) {
        this.query = query;
        this.target = target;
        this.similarityReactions = similarityReactions;
    }

    /**
     * @return the query
     */
    public String getQuery() {
        return query;
    }

    /**
     * @return the target
     */
    public String getTarget() {
        return target;
    }

    /**
     * @return the Similarity Reactions Map
     */
    public Map<String, String> getSimilarityReactions() {
        return similarityReactions;
    }

}
