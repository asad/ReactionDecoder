/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.container.helper;

import java.io.Serializable;

import uk.ac.ebi.reactionblast.mapping.interfaces.IKey;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Key implements IKey, Serializable {

    private static final long serialVersionUID = 92392372979041041L;
    private final int sourceIndex;
    private final int targetIndex;

    /**
     *
     * @param sourceIndex
     * @param targetIndex
     */
    public Key(int sourceIndex, int targetIndex) {
        this.sourceIndex = sourceIndex;
        this.targetIndex = targetIndex;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Key:");
        sb.append(this.sourceIndex);
        sb.append(":");
        sb.append(this.targetIndex);
        return sb.toString();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Key other = (Key) obj;
        if (this.sourceIndex != other.sourceIndex) {
            return false;
        }
        return this.targetIndex == other.targetIndex;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 97 * hash + this.sourceIndex;
        hash = 97 * hash + this.targetIndex;
        return hash;
    }

    @Override
    public int compareTo(Key t) {
        final int BEFORE = -1;
        final int EQUAL = 0;
        final int AFTER = 1;
        String key1 = this.sourceIndex + "_" + this.targetIndex;
        String key2 = t.sourceIndex + "_" + t.targetIndex;

        if (key1.equals(key2)) {
            return EQUAL;
        } else {
            return key1.compareTo(key2);
        }
    }

    /**
     *
     * @param t
     * @param t1
     * @return
     */
    @Override
    public int compare(Key t, Key t1) {
        return t.compareTo(t1);
    }
}
