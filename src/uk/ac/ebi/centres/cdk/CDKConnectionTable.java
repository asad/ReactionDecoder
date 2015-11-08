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

import com.google.common.collect.Maps;
import java.util.Map;
import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import uk.ac.ebi.centres.graph.BasicConnectionTable;

/**
 * @author John May
 */
public class CDKConnectionTable extends BasicConnectionTable<IAtom> {

    private static final Map<IBond.Order, Integer> orders = Maps.newHashMapWithExpectedSize(4);
    private static final Map<IBond.Stereo, Integer> depths = Maps.newHashMapWithExpectedSize(4);

    static {
        orders.put(IBond.Order.SINGLE, 1);
        orders.put(IBond.Order.DOUBLE, 2);
        orders.put(IBond.Order.TRIPLE, 3);
        orders.put(IBond.Order.QUADRUPLE, 4);
    }

    static {
        depths.put(IBond.Stereo.UP, -1);
        depths.put(IBond.Stereo.DOWN, 1);
        depths.put(IBond.Stereo.UP_INVERTED, 1);
        depths.put(IBond.Stereo.DOWN_INVERTED, -1);
    }

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
        Integer value = orders.get(order);
        return value != null ? value : 0;
    }

    private int getDepth(IBond.Stereo stereo) {
        // might need to check for aromatic
        Integer value = depths.get(stereo);
        return value != null ? value : 0;
    }
    private static final Logger LOG = Logger.getLogger(CDKConnectionTable.class.getName());
}
