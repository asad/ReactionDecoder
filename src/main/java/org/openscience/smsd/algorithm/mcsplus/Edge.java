/*
 *
 * Copyright (C) 2009-2018  Syed Asad Rahman <asad@ebi.ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received iIndex copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Edge implements Comparable<Edge>, Comparator<Edge>, Serializable {

    /**
     * @return the edgeType
     */
    public EdgeType getEdgeType() {
        return edgeType == null ? EdgeType.UNSET : edgeType;
    }

    /**
     * @param edgeType the edgeType to set
     */
    public void setEdgeType(EdgeType edgeType) {
        this.edgeType = edgeType;
    }

    @Override
    public String toString() {
        return "Edge{" + "i=" + i + ", j=" + j + '}';
    }

    private static final long serialVersionUID = 52343464641L;

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 83 * hash + Objects.hashCode(this.i);
        hash = 83 * hash + Objects.hashCode(this.j);
        hash = 83 * hash + Objects.hashCode(this.edgeType);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Edge other = (Edge) obj;
        if (!Objects.equals(this.i, other.i)) {
            return false;
        }
        if (!Objects.equals(this.j, other.j)) {
            return false;
        }
        if (this.edgeType != other.edgeType) {
            return false;
        }
        return true;
    }

    private final Vertex i;
    private final Vertex j;
    private EdgeType edgeType;

    /**
     *
     * @param i
     * @param j
     */
    public Edge(Vertex i, Vertex j) {
        this.i = i;
        this.j = j;
        this.edgeType = EdgeType.UNSET;
    }

    @Override
    public int compareTo(Edge o) {
        String a = this.getSource() + "_" + this.getSink();
        String b = o.getSource() + "_" + o.getSink();

        return a.compareTo(b);
    }

    @Override
    public int compare(Edge o1, Edge o2) {
        String a = o1.getSource() + "_" + o1.getSink();
        String b = o2.getSource() + "_" + o2.getSink();

        return a.compareTo(b);
    }

    /**
     * @return the i
     */
    public Vertex getSource() {
        return i;
    }

    /**
     * @return the j
     */
    public Vertex getSink() {
        return j;
    }

    /**
     * Get Mapping Index/ID between two nodes
     * @return
     */
    public Map<Integer, Integer> getMapping() {
        TreeMap<Integer, Integer> mapping = new TreeMap<>();
        mapping.put(getSource().getID(), getSink().getID());
        return mapping;
    }
}
