/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph;

import java.util.Collection;
import java.util.Set;
import java.util.Stack;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public interface IClique {

    /**
     *
     * R: set of vertices belonging to the current clique
     *
     * X: set of vertices which are not allowed to be added to R, defined as X
     * in paper
     *
     * P: is a set of vertices which <b>can</b> be added to R, because they are
     * neighbours of vertex u via <i>c-edges</i>
     *
     * Q: is a set of vertices which <b>cannot</b> be added to R, because they
     * are neighbours of vertex u via <i>d-edges</i>
     *
     * Vertex: stored all the vertices for the Graph G Vertex[G]: nodes of
     * vector graph are stored in Vertex
     *
     */
    void findMaximalCliques();

    /**
     *
     * @return Collection of cliques (each of which is represented as a Set of
     * vertices)
     */
    Collection<Set<Vertex>> getCliques();

    /**
     * Finds the largest maximal cliques of the graph.
     *
     * @return the largest cliques
     */
    Stack<Set<Vertex>> getMaxCliquesSet();

}
