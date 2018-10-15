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
package uk.ac.ebi.centres.priority.access;

/**
 * Wrapper provides access to a given an atom's atomic number (of type A) to
 * priority rules that require it. The method can be injected into a priority to
 * rule too allow that rule to act on the any given atom class type. This
 * special case flips the hydrogen and pseudo atom priority h = 0 and R = 1 so
 * that implicit and explicit centres give the same label.
 *
 * @param <A> the atom class type
 * @author John May
 */
public class PsuedoAtomicNumberModifier<A> implements AtomicNumberAccessor<A> {

    private final AtomicNumberAccessor<A> accessor;

    /**
     *
     * @param accessor
     */
    public PsuedoAtomicNumberModifier(AtomicNumberAccessor<A> accessor) {
        this.accessor = accessor;
    }

    /**
     * Access the atomic number for a provided atom. The atomic can be 0 to
     * allow for cases such as 'R' or '*' but should never be negative.
     *
     * @param atom the atom to access the atomic number for
     * @return a positive integer value that is the atomic number for the given
     * atom
     */
    @Override
    public int getAtomicNumber(A atom) {
        int value = accessor.getAtomicNumber(atom);

        switch (value) {
            case 1:
                return 0;
            case 0:
                return 1;
            default:
                return value;
        }

    }

}
