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
package uk.ac.ebi.centres.graph;

import java.util.AbstractMap;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.ConnectionTable;

/**
 * @author John May
 * @param <A>
 */
public class BasicConnectionTable<A> implements ConnectionTable<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(BasicConnectionTable.class);

    private final Map<A, Map<A, Map.Entry<Integer, Integer>>> connections = new HashMap<>();
    private final Map<A, Map<A, Map.Entry<Integer, Integer>>> stereo = new HashMap<>();

    /**
     *
     * @param first
     * @param second
     * @param order
     */
    public void addConnection(A first, A second, int order) {
        addConnection(first, second, order, 0);
    }

    /**
     *
     * @param first
     * @param second
     * @param order
     * @param sign
     */
    public void addConnection(A first, A second, int order, int sign) {
        newConnection(first, second, order,
                sign >= 1 ? 1 : sign <= -1 ? -1 : 0);
        newConnection(second, first, order,
                sign >= 1 ? -1 : sign <= -1 ? 1
                                : 0); // note the sign is inverted
    }

    private void newConnection(A first, A second, int order, int sign) {
        if (!connections.containsKey(first)) {
            connections.put(first, new HashMap<>());
        }
        connections.get(first).put(second, new AbstractMap.SimpleEntry<>(order, sign));
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public Collection<A> getConnected(A atom) {
        return connections.get(atom).keySet();
    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    @Override
    public int getOrder(A first, A second) {
        return connections.get(first).get(second).getKey();
    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    @Override
    public Integer getDepth(A first, A second) {
        return connections.get(first).get(second).getValue();
    }

    /**
     *
     * @return
     */
    @Override
    public Integer getAtomCount() {
        return connections.keySet().size();
    }
}
