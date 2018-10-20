/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class stores the results of the compatibility graph between query and
 * target molecule. It also marks edges in the compatibility graph as c-edges or
 * d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
/**
 * Stores the Result (c-edges, d-edges and Graph Nodes)
 */
public class Result {

    final List<Edge> cEdges;
    final List<Edge> dEdges;
    final List<Integer> compGraphNodes;

    public Result() {
        cEdges = Collections.synchronizedList(new ArrayList<>());
        dEdges = Collections.synchronizedList(new ArrayList<>());
        compGraphNodes = Collections.synchronizedList(new ArrayList<>());
    }

    public synchronized List<Edge> getCEgdes() {
        return Collections.synchronizedList(cEdges);
    }

    public synchronized List<Edge> getDEgdes() {
        return Collections.synchronizedList(dEdges);
    }

    /**
     * @return the Compatibility Graph Nodes
     */
    public List<Integer> getCompGraphNodes() {
        return compGraphNodes;
    }
}
