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
package uk.ac.ebi.centres.descriptor;

import uk.ac.ebi.centres.Descriptor;

/**
 * Enumeration of asymmetric and pseudo-asymmetric tetrahedral descriptors. Lower case indicates a descriptor are
 * pseudo-asymmetric.
 *
 * @author John May
 * @see Planar
 * @see Trigonal
 * @see General
 */
public enum Tetrahedral
        implements Descriptor {

    /**
     * Indicates the priority of ligands around the chiral atoms proceeds in a clockwise rotation.
     */
    R(Type.ASYMMETRIC),
    /**
     * Indicates the priority of ligands around the chiral atoms proceeds in a anti-clockwise rotation.
     */
    S(Type.ASYMMETRIC),
    /**
     * Indicates the priority of ligands around the chiral atoms proceeds in a clockwise rotation. The priority is
     * pseudo-asymmetric if the centre is only defined by opposing stereo descriptors.
     *
     * @see <a href="http://goldbook.iupac.org/P04921.html">pseudo-asymmetric carbon atom</a>
     */
    r(Type.PSEUDO_ASYMMETRIC),
    /**
     * Indicates the priority of ligands around the chiral atoms proceeds in a anti-clockwise rotation. The priority is
     * pseudo-asymmetric if the centre is only defined by opposing stereo descriptors.
     *
     * @see <a href="http://goldbook.iupac.org/P04921.html">pseudo-asymmetric carbon atom</a>
     */
    s(Type.PSEUDO_ASYMMETRIC);
    private final Type type;

    private Tetrahedral(Type type) {
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
