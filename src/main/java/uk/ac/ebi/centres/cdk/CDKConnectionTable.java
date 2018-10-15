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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.cdk;

import static com.google.common.collect.Maps.newHashMapWithExpectedSize;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN_INVERTED;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_INVERTED;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.graph.BasicConnectionTable;

/**
 * @author John May
 */
public class CDKConnectionTable extends BasicConnectionTable<IAtom> {

    private static final Map<IBond.Order, Integer> ORDERS = newHashMapWithExpectedSize(4);
    private static final Map<IBond.Stereo, Integer> DEPTHS = newHashMapWithExpectedSize(4);
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CDKConnectionTable.class);

    static {
        ORDERS.put(SINGLE, 1);
        ORDERS.put(DOUBLE, 2);
        ORDERS.put(TRIPLE, 3);
        ORDERS.put(QUADRUPLE, 4);
    }

    static {
        DEPTHS.put(UP, -1);
        DEPTHS.put(DOWN, 1);
        DEPTHS.put(UP_INVERTED, 1);
        DEPTHS.put(DOWN_INVERTED, -1);
    }

    /**
     *
     * @param container
     */
    public CDKConnectionTable(IAtomContainer container) {
        for (int i = 0; i < container.getAtomCount(); i++) {
            container.getAtom(i).setProperty("number", i + 1);
        }
        for (IBond bond : container.bonds()) {
            addConnection(bond.getAtom(0),
                    bond.getAtom(1),
                    getOrder(bond.getOrder()), // might need to check for aromatic
                    getDepth(bond.getStereo()));

        }
    }

    private int getOrder(IBond.Order order) {
        Integer value = ORDERS.get(order);
        return value != null ? value : 0;
    }

    private int getDepth(IBond.Stereo stereo) {
        // might need to check for aromatic
        Integer value = DEPTHS.get(stereo);
        return value != null ? value : 0;
    }
}
