/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph;

import java.io.Serializable;
import java.util.Comparator;

/**
 * This class holds Nodes (vertex)
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class Vertex implements Comparable<Vertex>, Comparator<Vertex>, Serializable {

    private Integer query;
    private Integer target;

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 37 * hash + this.id;
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
        final Vertex other = (Vertex) obj;
        return this.id == other.id;
    }

    @Override
    public String toString() {
        if (label.isEmpty()) {
            return "V{" + "id=" + id + "(" + query + "," + target + ")" + '}';
        } else if (query == null) {
            return "V{" + "id=" + id + ", label=" + label + '}';
        }
        return "V{" + "Q=" + query + ", T=" + target + ", id=" + id + ", label=" + label + '}';
    }

    /**
     * @return the label
     */
    public String getLabel() {
        return label;
    }

    /**
     * @param label the label to set
     */
    public void setLabel(String label) {
        this.label = label;
    }

    private final int id;
    private String label;

    public Vertex(int node) {
        this.id = node;
        this.label = "";
        this.query = null;
        this.target = null;
    }

    @Override
    public int compareTo(Vertex o) {
        return this.getID() - o.getID();
    }

    @Override
    public int compare(Vertex o1, Vertex o2) {
        return o1.getID() - o2.getID();
    }

    /**
     * @return the id
     */
    public int getID() {
        return id;
    }

    /**
     *
     * @param a
     * @param b
     */
    public void setCompatibilityBondPair(Integer a, Integer b) {
        this.query = a;
        this.target = b;
    }

    /**
     *
     * @return
     */
    public Integer getTargetBondIndex() {
        return target;
    }

    /**
     *
     * @return
     */
    public Integer getQueryBondIndex() {
        return query;
    }
}
