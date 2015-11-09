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
import java.util.Set;
import uk.ac.ebi.centres.graph.Arc;

/**
 * Defines a ligand which in the digraph is a side chain of a central atom. The ligand is a node in the {@link Digraph}.
 * Each ligand provides the atom that this ligand represents as well as the child ligands. The ligand also allows
 * determination of parent and visited atoms.
 *
 * @author John May
 * @param <A>
 * @see Digraph
 */
public interface Ligand<A> {

    /**
     *
     * @return
     */
    public boolean isTerminal();

    /**
     *
     * @return
     */
    public boolean isBranching();

    /**
     *
     * @return
     */
    public boolean isDuplicate();

    /**
     *
     * @return
     */
    public A getAtom();

    /**
     * Resets any caches
     */
    public void reset();

    /**
     * Sets the descriptor for this centre. This method will throw an illegal argument exception if the descriptor is
     * set to null.
     *
     * @param descriptor the new descriptor for this centre
     *
     * @see uk.ac.ebi.centres.descriptor.General
     * @see uk.ac.ebi.centres.descriptor.Tetrahedral
     * @see uk.ac.ebi.centres.descriptor.Planar
     * @see uk.ac.ebi.centres.descriptor.Trigonal
     */
    public void setDescriptor(Descriptor descriptor);

    /**
     * Access the descriptor for this centre. This descriptor is the primary descriptor for this centre and not an
     * auxiliary descriptor. Auxiliary descriptors should be set on a per ligand basis. This method should not return
     * null and instead return {@link uk.ac.ebi.centres.descriptor.General#UNKNOWN} for unknown/not yet determined
     * centres.
     *
     * @return descriptor for this centre
     *
     * @see uk.ac.ebi.centres.descriptor.General
     * @see uk.ac.ebi.centres.descriptor.Tetrahedral
     * @see uk.ac.ebi.centres.descriptor.Planar
     * @see uk.ac.ebi.centres.descriptor.Trigonal
     */
    public Descriptor getDescriptor();

    /**
     * Access the child ligands. Child ligands are ligands that are further away from the central atom in sphere n+1.
     *
     * @return a list of the child ligands
     */
    public List<Ligand<A>> getLigands();

    /**
     * Access all previously visited atoms. For convenience the {@link
     * #isVisited} method can be used to determine if an atom has already been visited.
     *
     * @return previously visited atoms
     *
     * @see #isVisited(Object)
     */
    public Set<A> getVisited();

    /**
     * Determines if the provided atom is this ligands parent. There is only ever a single parent atom which is closer
     * to the root of the graph. Normally the root ligand will have it's parent set to itself.
     *
     * @param atom a potential parent
     *
     * @return whether the atom is the parent (n-1) of this ligand
     */
    public Boolean isParent(A atom);

    /**
     * Set the parent of ligand for when we are shuffling around for auxiliary descriptors. Not supported on planar
     * ligands
     *
     * @param atom
     */
    public void setParent(A atom);

    /**
     *
     * @return
     */
    public A getParent();

    /**
     * Determine if the provided atom has already been visited
     *
     * @param atom the atom which may be visited
     *
     * @return whether the provided atom has been visited
     */
    public Boolean isVisited(A atom);

    /**
     * Sets the auxiliary descriptor for this ligand
     *
     * @param descriptor the new auxiliary descriptor
     */
    public void setAuxiliary(Descriptor descriptor);

    /**
     * Access the auxiliary descriptor for this ligand.
     *
     * @return the auxiliary descriptor
     */
    public Descriptor getAuxiliary();

    /**
     *
     * @return
     */
    public List<Arc<A>> getArcs();

    /**
     *
     * @return
     */
    public Arc<A> getParentArc();

    /**
     * Access the distance from the root this ligand is
     *
     * @return
     */
    public int getDistanceFromRoot();

    /**
     * Access the depth (z coord)
     *
     * @return
     */
    public int getDepth();
}
