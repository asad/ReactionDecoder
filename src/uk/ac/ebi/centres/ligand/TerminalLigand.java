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

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import uk.ac.ebi.centres.ConnectionProvider;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;

/**
 * @author John May
 */
public class TerminalLigand<A> extends NonterminalLigand<A> {

    public TerminalLigand(ConnectionProvider<A> provider, MutableDescriptor descriptor, A atom, A parent, int distance) {
        super(provider, descriptor, atom, parent, distance);
        setDuplicate(Boolean.TRUE);
    }

    public TerminalLigand(ConnectionProvider<A> provider, MutableDescriptor descriptor, Set<A> visited, A atom, A parent, int distance) {
        super(provider, descriptor, visited, atom, parent, distance);
        setDuplicate(Boolean.TRUE);
    }

    @Override
    public List<Ligand<A>> getLigands() {
        // suppress use of connection provider
        return Collections.emptyList();
    }

    @Override
    public String toString() {
        return super.toString() + "'";
    }

    @Override
    public boolean isTerminal() {
        return Boolean.TRUE;
    }

    @Override
    public boolean isBranching() {
        return Boolean.FALSE;
    }
    private static final Logger LOG = Logger.getLogger(TerminalLigand.class.getName());
}
