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

import org.openscience.cdk.interfaces.IAtom;
import com.bioinceptionlabs.reactionblast.model.AtomNode;

/**
 * CDK adapter for AtomNode. Wraps a CDK IAtom as a graph node.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CDKAtomNode implements AtomNode {

    private final IAtom cdkAtom;

    public CDKAtomNode(IAtom cdkAtom) {
        if (cdkAtom == null) throw new IllegalArgumentException("CDK atom cannot be null");
        this.cdkAtom = cdkAtom;
    }

    public IAtom getCDKAtom() {
        return cdkAtom;
    }

    @Override public String getSymbol() { return cdkAtom.getSymbol(); }
    @Override public int getAtomicNumber() { return cdkAtom.getAtomicNumber() != null ? cdkAtom.getAtomicNumber() : 0; }
    @Override public Integer getFormalCharge() { return cdkAtom.getFormalCharge(); }
    @Override public Integer getMassNumber() { return cdkAtom.getMassNumber(); }
    @Override public boolean isAromatic() { return cdkAtom.isAromatic(); }
    @Override public void setAromatic(boolean aromatic) { cdkAtom.setIsAromatic(aromatic); }
    @Override public Integer getImplicitHydrogenCount() { return cdkAtom.getImplicitHydrogenCount(); }
    @Override public String getId() { return cdkAtom.getID(); }
    @Override public void setId(String id) { cdkAtom.setID(id); }
    @Override public Object getProperty(String key) { return cdkAtom.getProperty(key); }
    @Override public void setProperty(String key, Object value) { cdkAtom.setProperty(key, value); }
    @Override public boolean getFlag(int flag) { return cdkAtom.getFlag(flag); }
    @Override public void setFlag(int flag, boolean value) { cdkAtom.setFlag(flag, value); }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof CDKAtomNode) return cdkAtom == ((CDKAtomNode) o).cdkAtom;
        return false;
    }

    @Override
    public int hashCode() {
        return System.identityHashCode(cdkAtom);
    }

    @Override
    public String toString() {
        return getSymbol() + (getId() != null ? ":" + getId() : "");
    }
}
