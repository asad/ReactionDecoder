/*
 *
 * Copyright (C) 2009-2020  Syed Asad Rahman <asad@ebi.ebi.ac.uk>
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
package org.openscience.smsd.graph;

import java.io.Serializable;
import java.util.Objects;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Edge implements Serializable {

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 11 * hash + Objects.hashCode(this.source);
        hash = 11 * hash + Objects.hashCode(this.sink);
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
        if (!Objects.equals(this.source, other.source)) {
            return false;
        }
        return Objects.equals(this.sink, other.sink);
    }

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
        return "Edge{" + "i=" + source + ", j=" + sink + '}';
    }

    private static final long serialVersionUID = 52343464641L;

    private final Integer source;
    private final Integer sink;
    private EdgeType edgeType;

    /**
     *
     * @param i
     * @param j
     */
    public Edge(Integer i, Integer j) {
        this.source = i;
        this.sink = j;
        this.edgeType = EdgeType.UNSET;
    }

    /**
     * @return the source
     */
    public Integer getSource() {
        return source;
    }

    /**
     * @return the sink
     */
    public Integer getSink() {
        return sink;
    }

    public boolean isC_Edge() {
        return this.edgeType == EdgeType.C_EDGE;
    }

    public boolean isD_Edge() {
        return this.edgeType == EdgeType.D_EDGE;
    }
}
