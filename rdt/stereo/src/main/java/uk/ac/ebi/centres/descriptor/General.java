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
package uk.ac.ebi.centres.descriptor;

import uk.ac.ebi.centres.Descriptor;

/**
 * Enumeration of general descriptors. These general descriptors are not specific to any particular stereo-centre
 * (tetrahedral, planar or trigonal) and can also be used to indicate non-stereo-centres.
 *
 * @author John May
 * @see Tetrahedral
 * @see Trigonal
 * @see Planar
 */
public enum General
        implements Descriptor {

    /**
     * A centre which is asymmetric but is lacking information to define a dull descriptor. An example could be a 2D
     * tetrahedral centre that does not have a up/down bond.
     */
    UNSPECIFIED(Type.NON_STEREOGENIC),
    /**
     * Indicates an atom is known not to be a stereo-centre
     */
    NONE(Type.NON_STEREOGENIC),
    /**
     * Indicates that it is unknown whether the atom is a stereo-centre
     */
    UNKNOWN(Type.NON_STEREOGENIC);
    private final Type type;

    private General(Type type) {
        this.type = type;
    }

    /**
     * @inheritDoc
     */
    @Override
    public Type getType() {
        return type;
    }
}
