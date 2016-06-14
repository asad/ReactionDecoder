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

import java.util.List;
import uk.ac.ebi.centres.graph.Arc;

/**
 * @author John May
 * @param <A>
 */
public interface ConnectionProvider<A> {

    /**
     *
     */
    public void build();

    /**
     *
     * @param atom
     * @return
     */
    public List<Ligand<A>> getLigands(A atom);

    /**
     *
     * @param ligand
     * @return
     */
    public List<Ligand<A>> getLigands(Ligand<A> ligand);

    /**
     *
     * @param ligand
     * @return
     */
    public List<Arc<A>> getArcs(Ligand<A> ligand);

    /**
     *
     * @param ligand
     * @return
     */
    public Arc<A> getParentArc(Ligand<A> ligand);

    /**
     *
     * @param ligand
     */
    public void reroot(Ligand<A> ligand);

    /**
     * Clear the digraph ready for GC
     */
    public void dispose();
}
