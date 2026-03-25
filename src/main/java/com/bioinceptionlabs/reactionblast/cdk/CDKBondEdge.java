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
package com.bioinceptionlabs.reactionblast.cdk;

import org.openscience.cdk.interfaces.IBond;
import com.bioinceptionlabs.reactionblast.model.AtomNode;
import com.bioinceptionlabs.reactionblast.model.BondEdge;

/**
 * CDK adapter for BondEdge. Wraps a CDK IBond as a graph edge.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CDKBondEdge implements BondEdge {

    private final IBond cdkBond;
    private final CDKAtomNode source;
    private final CDKAtomNode target;

    public CDKBondEdge(IBond cdkBond) {
        if (cdkBond == null) throw new IllegalArgumentException("CDK bond cannot be null");
        this.cdkBond = cdkBond;
        this.source = new CDKAtomNode(cdkBond.getBegin());
        this.target = new CDKAtomNode(cdkBond.getEnd());
    }

    public IBond getCDKBond() {
        return cdkBond;
    }

    @Override public AtomNode getSource() { return source; }
    @Override public AtomNode getTarget() { return target; }

    @Override
    public BondOrder getOrder() {
        if (cdkBond.getOrder() == null) return BondOrder.UNSET;
        switch (cdkBond.getOrder()) {
            case SINGLE: return BondOrder.SINGLE;
            case DOUBLE: return BondOrder.DOUBLE;
            case TRIPLE: return BondOrder.TRIPLE;
            case QUADRUPLE: return BondOrder.QUADRUPLE;
            default: return BondOrder.UNSET;
        }
    }

    @Override
    public void setOrder(BondOrder order) {
        switch (order) {
            case SINGLE: cdkBond.setOrder(IBond.Order.SINGLE); break;
            case DOUBLE: cdkBond.setOrder(IBond.Order.DOUBLE); break;
            case TRIPLE: cdkBond.setOrder(IBond.Order.TRIPLE); break;
            case QUADRUPLE: cdkBond.setOrder(IBond.Order.QUADRUPLE); break;
            default: cdkBond.setOrder(IBond.Order.UNSET); break;
        }
    }

    @Override public boolean isAromatic() { return cdkBond.isAromatic(); }
    @Override public void setAromatic(boolean aromatic) { cdkBond.setIsAromatic(aromatic); }

    @Override
    public boolean connects(AtomNode atom) {
        if (atom instanceof CDKAtomNode) {
            return cdkBond.contains(((CDKAtomNode) atom).getCDKAtom());
        }
        return false;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof CDKBondEdge) return cdkBond == ((CDKBondEdge) o).cdkBond;
        return false;
    }

    @Override
    public int hashCode() {
        return System.identityHashCode(cdkBond);
    }
}
