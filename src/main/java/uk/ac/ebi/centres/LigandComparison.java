/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

package uk.ac.ebi.centres;

import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * Simple holder for a ligand comparison. The comparison holds the value
 * produced in the comparison as well as the type of the comparison.
 *
 * @author John May
 */
public class LigandComparison implements Comparison {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(LigandComparison.class);

    private final Integer         order;
    private final Descriptor.Type type;


    /**
     * Construct a new ligand comparison with order and type.
     *
     * @param order the order of two ligands
     * @param type  the type of the comparison
     *
     * @see PriorityRule
     * @see java.util.Comparator
     */
    public LigandComparison(Integer order, Descriptor.Type type) {
        this.order = order;
        this.type = type;
    }


    /**
     * @inheritDoc
     */
    @Override
    public Integer getOrder() {
        return order;
    }


    /**
     * @inheritDoc
     */
    @Override
    public Descriptor.Type getType() {
        return type;
    }
}
