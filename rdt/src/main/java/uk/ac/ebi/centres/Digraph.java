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

/**
 * Defines a directed graph of ligands. Each graph has an atom which provides the root of the graph. The graph is
 * normally build on demand with the root ligand being the only ligand that always exists. The graph also holds onto
 * centre's of chirality via the {@link Centre}. Centre's can be accessed via and atom for tetrahedral/trigonal centres
 * or via two atoms (indicating a bond) for planar centres. This centre access allows perception algorithms to control
 * the mutation and access of descriptors.
 * <p/>
 * All digraphs are acyclic with ghost atoms being created for cyclic graphs.
 *
 * @author John May
 * @param <A>
 * @see Centre
 */
public interface Digraph<A> {

    /**
     * Access the root ligand of this digraph. The root is normally the chiral atom that is being determined. When
     * determining ligands around double bond, two roots are normally needed. These two roots will provide each other as
     * ligands and thus should be filtered out (currently outside of digraph).
     *
     * @return the root ligand of the directed graph
     */
    public Ligand<A> getRoot();

    /**
     * Access the proximal ligands next to the root. If no root is set this method will throw an
     * {@link IllegalArgumentException}.
     *
     * @return list of proximal ligands
     *
     * @see Ligand
     */
    public List<Ligand<A>> getProximal();

    /**
     * Access all created ligands for the provided atom. In acyclic structure with only single bonds there will only be
     * a single ligand per atom. Cyclic structures and double bonds create 'ghost' ligands in the structure which can
     * also be accessed. It is important to note that the method will only return ligands for atoms that have been
     * built.
     * <p/>
     * As the digraph is normally constructed on a per-ligand basis an atom may have no ligands if they haven't been
     * visited yet. This can be fixed by invoking {@link #build()}.
     *
     * @param atom to access ligands for
     *
     * @return a list of ligands that have been constructed for the provided atom
     *
     * @see #build()
     */
    public List<Ligand<A>> getLigands(A atom);

    /**
     * Exhaustively expands from the root creating all ligands. Normally the graph is constructor on a per-ligand basis
     * starting at the root. When using auxiliary descriptors expanding the whole graph is often required as centres
     * need to be determined on remote regions.
     *
     * @see #getLigands(A)
     */
    public void build();

    /**
     * Reroot the digraph on the provided ligand. This does not recalculate the whole graph (each centre should have
     * it's own graph) but instead transforms the directions of the the edges to point away from the new root. This
     * method is primarily used for generating auxiliary descriptors.
     *
     * @param ligand the ligand which will be the new root
     */
    public void reroot(Ligand<A> ligand);
}
