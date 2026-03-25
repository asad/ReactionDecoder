/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package com.bioinceptionlabs.reactionblast.model;

/**
 * Graph edge representing a chemical bond. Toolkit-agnostic — implementations
 * wrap CDK IBond, RDKit Bond, OpenBabel OBBond, etc.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface BondEdge {

    AtomNode getSource();

    AtomNode getTarget();

    BondOrder getOrder();

    void setOrder(BondOrder order);

    boolean isAromatic();

    void setAromatic(boolean aromatic);

    boolean connects(AtomNode atom);

    enum BondOrder {
        SINGLE(1), DOUBLE(2), TRIPLE(3), QUADRUPLE(4), UNSET(0);

        private final int numeric;

        BondOrder(int numeric) {
            this.numeric = numeric;
        }

        public int numeric() {
            return numeric;
        }
    }
}
