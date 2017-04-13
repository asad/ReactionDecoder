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

import java.util.Comparator;
import java.util.List;

/**
 * Defines a comparator for ligands. The default {@link Comparator} can be used to prioritise ligands however the
 * {@link #compareLigands(Ligand, Ligand)} also adds meta data about the type of the descriptor via the {@link
 * Comparison}
 *
 * @author John May
 * @param <A>
 * @see Comparison
 * @see Comparator
 */
public interface PriorityRule<A> extends Comparator<Ligand<A>> {

    /**
     *
     */
    public enum Type {

        /**
         *
         */
        CONSTITUTIONAL,

        /**
         *
         */
        GEOMETRICAL,

        /**
         *
         */
        TOPOGRAPHICAL,

        /**
         *
         */
        COMBINED
    }

    /**
     *
     * @return
     */
    public Type getRuleType();

    /**
     *
     * @param o1
     * @param o2
     * @return
     */
    public int recursiveCompare(Ligand<A> o1, Ligand<A> o2);

    /**
     * Prioritises ligands using the provided sorter and indicates whether the ligands were unique.
     *
     * @param ligands a list of ligands to prioritise
     *
     * @return whether the ligands were unique
     */
    public Priority prioritise(List<Ligand<A>> ligands);

    /**
     * Allows injection of a ligand sorter. The ligand sort is used when two compared ligands are ranked equally. The
     * sorter needs to be injected for as when combination of priority rules is required the sorting must be done all
     * proceeding rules.
     *
     * @param sorter the ligand sorter to use
     */
    public void setSorter(LigandSorter<A> sorter);

    /**
     * Analogous to {@link #compare(Object, Object)} the prioritise method combines the {@link Descriptor.Type} to the
     * order and can indicate what comparison method was used. The single comparison cases the type doesn't change
     * however when using a combined comparator the type may change depending on which comparator was used.
     *
     * @param o1 first ligand
     * @param o2 second ligand
     *
     * @return the order of the two objects
     *
     * @see #compare(Object, Object)
     */
    public Comparison compareLigands(Ligand<A> o1, Ligand<A> o2);

    /**
     * Access the descriptor type this rule indicates. Normally rules will indicate. In rare cases a rule produce a
     * pseudo-asymmetric centre.
     *
     * @return the type of rule
     */
    public Descriptor.Type getType();

    /**
     * Indicates the rule should halt. This allows us to terminate a timed out thread but stopping all comparisons
     *
     * @param value
     */
    public void setHalt(boolean value);
}
