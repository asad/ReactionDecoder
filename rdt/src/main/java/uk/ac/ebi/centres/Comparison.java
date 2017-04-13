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

/**
 * Wrapper for a comparison between two ligands. The comparison hold the value (order) from the
 * {@link java.util.Comparator} and the type of the comparison {@link uk.ac.ebi.centres.Descriptor.Type}
 *
 * @author John May
 */
public interface Comparison {

    /**
     * Access the order of the comparison between two ligands. The order is the value returned by the
     * {@link java.util.Comparator}. Negative value indicates the first ligand ranks lower then second whilst a positive
     * value indicated the second ligand ranks lower then the first. A value of 0 indicates the ligands order equally.
     *
     * @return the order of two compared ligands
     *
     * @see java.util.Comparator
     */
    public Integer getOrder();

    /**
     * Access the type of the comparison. The type is defined by the method used to prioritise. The main reason for the
     * type inclusion is some comparisons indicate pseudo-asymmetry. Normally the only rule used to define
     * pseudo-asymmetry is R proceeds S but integrating the type in the comparison makes for a consistent API and allows
     * the the same rule to be use across multiple threads.
     *
     * @return the inferred type that this comparison produces
     */
    public Descriptor.Type getType();
}
