/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

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
public class Node implements Comparable<Node>, Comparator<Node>, Serializable {

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 43 * hash + this.node;
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
        final Node other = (Node) obj;
        if (this.node != other.node) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        if (id.isEmpty()) {
            return "" + node + "";
        } else {
            return "Node{" + "node=" + node + ", id=" + id + '}';
        }
    }

    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(String id) {
        this.id = id;
    }

    final int node;
    private String id;

    public Node(int node) {
        this.node = node;
        this.id = "";
    }

    @Override
    public int compareTo(Node o) {
        return this.node - o.node;
    }

    @Override
    public int compare(Node o1, Node o2) {
        return o1.node - o2.node;
    }

}
