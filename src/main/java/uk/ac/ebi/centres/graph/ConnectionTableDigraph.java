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

import java.util.Collection;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.ConnectionTable;
import uk.ac.ebi.centres.DescriptorManager;
import uk.ac.ebi.centres.Ligand;

/**
 * @author John May
 * @param <A>
 */
public class ConnectionTableDigraph<A>
        extends AbstractDigraph<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(ConnectionTableDigraph.class);

    private final ConnectionTable<A> table;

    /**
     *
     * @param root
     * @param manager
     * @param table
     */
    public ConnectionTableDigraph(Ligand<A> root,
            DescriptorManager<A> manager,
            ConnectionTable<A> table) {
        super(root, manager);
        this.table = table;
    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    @Override
    public int getOrder(A first, A second) {
        return table.getOrder(first, second);
    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    @Override
    public int getDepth(A first, A second) {
        return table.getDepth(first, second);
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public Collection<A> getConnected(A atom) {
        return table.getConnected(atom);
    }
}
