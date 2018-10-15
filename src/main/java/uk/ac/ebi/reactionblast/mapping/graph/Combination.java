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
package uk.ac.ebi.reactionblast.mapping.graph;

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Combination implements Serializable, Comparable<Combination>, Comparator<Combination> {

    private static final long serialVersionUID = 786786786L;
    private final int row;
    private final int col;

    /**
     *
     * @param row
     * @param col
     */
    public Combination(int row, int col) {
        this.row = row;
        this.col = col;
    }

    /**
     *
     * @return
     */
    public int getRowIndex() {
        return row;
    }

    /**
     *
     * @return
     */
    public int getColIndex() {
        return col;
    }

    @Override
    public String toString() {
        return "Combination{" + "row=" + row + ", col=" + col + '}';
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 13 * hash + this.row;
        hash = 13 * hash + this.col;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Combination other = (Combination) obj;
        if (this.row != other.row) {
            return false;
        }
        return this.col == other.col;
    }

    @Override
    public int compareTo(Combination o) {
        String a = this.row + "_" + this.col;
        String b = o.row + "_" + o.col;
        return a.compareTo(b);
    }

    @Override
    public int compare(Combination o1, Combination o2) {
        String a = o1.row + "_" + o1.col;
        String b = o2.row + "_" + o2.col;
        return a.compareTo(b);
    }
}
