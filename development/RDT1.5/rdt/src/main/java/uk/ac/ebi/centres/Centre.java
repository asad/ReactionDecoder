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

import java.util.Collection;
import java.util.List;
import java.util.Set;

/**
 * Defines a stereo centre (normally on an atom or bond) that provides access and mutation of the centres descriptor.
 * This centre could plug directly into the molecular object (atom or bond) but would normally be a wrapper around the
 * molecular object which can then be transferred when all centres that can be perceived, have been perceived.
 *
 * @author John May
 * @param <A>
 * @see Descriptor
 * @see uk.ac.ebi.centres.descriptor.General
 * @see uk.ac.ebi.centres.descriptor.Tetrahedral
 * @see uk.ac.ebi.centres.descriptor.Planar
 * @see uk.ac.ebi.centres.descriptor.Trigonal
 */
public interface Centre<A> extends Ligand<A> {

    /**
     * Access the centre atoms that define this centre. In tetrahedral and trigonal centres this is a set of length one
     * whilst in planar centres this is a set of length two.
     *
     * @return the atoms of this centre
     */
    public Set<A> getAtoms();

    /**
     * Perceives the descriptor for this centre given a priority rule and a calculator for the sign of the space. The
     * descriptor is not set directly as descriptors are used in rules and should be exhaustively perceived before being
     * assigned. This allows descriptor calculation to be order independent.
     *
     * @param rule the rule by which this centre's proximal ligands are sorted
     * @param calculator the sign calculator to use (normally 2D or 3D).
     *
     * @return a perceived descriptor for this centre.
     */
    public Descriptor perceive(PriorityRule<A> rule, SignCalculator<A> calculator);

    /**
     *
     * @param proximal
     * @param rule
     * @param calculator
     * @return
     */
    public Descriptor perceive(List<Ligand<A>> proximal, PriorityRule<A> rule, SignCalculator<A> calculator);

    /**
     *
     * @param centres
     * @param rule
     * @param calculator
     * @return
     */
    public int perceiveAuxiliary(Collection<Centre<A>> centres,
            PriorityRule<A> rule,
            SignCalculator<A> calculator);

    /**
     * Clean up the digraph
     */
    public void dispose();
}
