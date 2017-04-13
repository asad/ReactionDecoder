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
package uk.ac.ebi.centres.priority.access;

/**
 * Wrapper provides access to a given an atom's mass number (of type A) to priority rules that require it. The method
 * can be injected into a priority to rule too allow that rule to act on the any given atom class type.
 *
 * @param <A> the atom class type
 *
 * @author John May
 */
public interface MassNumberAccessor<A> {

    /**
     * Access the mass number for a provided atom. The mass can be 0 but should never be negative.
     *
     * @param atom the atom to access the mass number for
     *
     * @return a positive integer value that is the mass number for the given atom
     */
    public int getMassNumber(A atom);
}
