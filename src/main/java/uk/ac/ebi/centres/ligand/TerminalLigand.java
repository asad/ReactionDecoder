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
package uk.ac.ebi.centres.ligand;

import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import static java.util.Collections.emptyList;
import java.util.List;
import java.util.Set;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.ConnectionProvider;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;

/**
 * @author John May
 * @param <A>
 */
public class TerminalLigand<A> extends NonterminalLigand<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(TerminalLigand.class);

    /**
     *
     * @param provider
     * @param descriptor
     * @param atom
     * @param parent
     * @param distance
     */
    public TerminalLigand(ConnectionProvider<A> provider, MutableDescriptor descriptor, A atom, A parent, int distance) {
        super(provider, descriptor, atom, parent, distance);
        setDuplicate(TRUE);
    }

    /**
     *
     * @param provider
     * @param descriptor
     * @param visited
     * @param atom
     * @param parent
     * @param distance
     */
    public TerminalLigand(ConnectionProvider<A> provider, MutableDescriptor descriptor, Set<A> visited, A atom, A parent, int distance) {
        super(provider, descriptor, visited, atom, parent, distance);
        setDuplicate(TRUE);
    }

    @Override
    public List<Ligand<A>> getLigands() {
        // suppress use of connection provider
        return emptyList();
    }

    @Override
    public String toString() {
        return super.toString() + "'";
    }

    /**
     *
     * @return
     */
    @Override
    public boolean isTerminal() {
        return TRUE;
    }

    /**
     *
     * @return
     */
    @Override
    public boolean isBranching() {
        return FALSE;
    }
}
