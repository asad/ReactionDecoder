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
 * Enumeration of asymmetric and pseudo-asymmetric planar descriptors. Planar descriptors are generally found around
 * double bonds which cannot rotate. The pseudo-asymmetric centres are indicated in lower case.
 *
 * @author John May
 * @see Tetrahedral
 * @see Trigonal
 * @see General
 */
public enum Planar
        implements Descriptor {

    /**
     * The two highest priority ligands are on the opposite (entgegen) side of a double bond. Also refereed to as
     * 'trans' isomerism.
     */
    E(Type.ASYMMETRIC),
    /**
     * The two highest priority ligands are on the together (zusammen) on the same side of a double bond. Also refereed
     * to as 'cis' isomerism.
     */
    Z(Type.ASYMMETRIC),
    /**
     * The two highest priority ligands are on the opposite (entgegen) side of a double bond. This centres is
     * pseudo-asymmetric and indicates the priority is only defined by opposite tetrahedral centres.
     */
    e(Type.PSEUDO_ASYMMETRIC),
    /**
     * The two highest priority ligands are on the together (zusammen) on the same side of a double bond. This centres
     * is pseudo-asymmetric and indicates the priority is only defined by opposite tetrahedral centres.
     */
    z(Type.PSEUDO_ASYMMETRIC);
    private final Type type;

    private Planar(Type type) {
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
