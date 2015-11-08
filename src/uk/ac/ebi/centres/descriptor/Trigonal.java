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
 * Enumeration of asymmetric and pseudo-asymmetric trigonal planar descriptors. These descriptors are generally found on
 * prochiral centres. The pseudo- asymmetric centres are indicated in lower case.
 *
 * @author John May
 * @see Planar
 * @see Tetrahedral
 * @see General
 */
public enum Trigonal
        implements Descriptor {

    /**
     * A trigonal prochiral centre whose ligands priority proceeded clockwise.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Prochirality">Prochirality</a>
     */
    Re(Type.ASYMMETRIC),
    /**
     * A trigonal prochiral centre whose ligands priority proceeded anti- clockwise.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Prochirality">Prochirality</a>
     */
    Si(Type.ASYMMETRIC),
    /**
     * A trigonal prochiral centre whose ligands priority proceeded clockwise. The priority of the ligands is only
     * defined by opposing stereo-centres and thus makes this centre pseudo-asymmetric.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Prochirality">Prochirality</a>
     */
    re(Type.PSEUDO_ASYMMETRIC),
    /**
     * A trigonal prochiral centre whose ligands priority proceeded anti- clockwise. The priority of the ligands is only
     * defined by opposing stereo-centres and thus makes this centre pseudo-asymmetric.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Prochirality">Prochirality</a>
     */
    si(Type.PSEUDO_ASYMMETRIC);
    private final Type type;

    private Trigonal(Type type) {
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
