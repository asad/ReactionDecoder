/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public class GraphKoch implements IClique {

    private final static boolean DEBUG = false;
    private final Collection<Set<Vertex>> cliques;

    private final Set<Edge> cEdges;
    private final Set<Edge> dEdges;
    private final Graph graph;

    private Integer iterations;

    /**
     *
     * @return Collection of cliques (each of which is represented as a Set of
     * vertices)
     */
    @Override
    public Collection<Set<Vertex>> getCliques() {
        return cliques;
    }

    /**
     * Finds the largest maximal cliques of the graph.
     *
     * @return the largest cliques
     */
    @Override
    public Stack<Set<Vertex>> getMaxCliquesSet() {
        Stack<Set<Vertex>> maxCliquesSet = new Stack<>();
        int best_clique_size = 0;
        for (Set<Vertex> clique : cliques) {
            if (clique.size() >= best_clique_size) {
                if (clique.size() > best_clique_size) {
                    while (!maxCliquesSet.empty()) {
                        maxCliquesSet.pop();
                    }
                    best_clique_size = clique.size();
                }
                if (clique.size() == best_clique_size) {
                    maxCliquesSet.push(clique);
                }
            }
        }
        return maxCliquesSet;
    }

    /**
     *
     * @param comp_graph_nodes
     * @param cEdges
     * @param dEdges
     */
    public GraphKoch(Graph comp_graph_nodes, Set<Edge> cEdges, Set<Edge> dEdges) {
        this.cEdges = cEdges;
        this.dEdges = dEdges;
        this.graph = comp_graph_nodes;
        this.cliques = new HashSet<>();
        this.iterations = 0;
    }

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
    @Override
    public void findMaximalCliques() {
        Set<Vertex> result = new LinkedHashSet<>();
        //set of vertices which have already been used for the initialization of Enumerate_C_Cliques()
        Set<Vertex> T = new LinkedHashSet<>();		// T <- Empty
        Set<Vertex> P, D, N, S;
        int currentmaxresult = 0;
        for (Vertex u : graph.nodes()) {				//for all u ELEMENTOF Vertex

            P = new LinkedHashSet<>();			// P <- Empty
            D = new LinkedHashSet<>();			// D <- Empty
            S = new LinkedHashSet<>();			// S <- Empty
            N = findNeighbors(u);	// N <- {v ELEMENTOF Vertex | {u,v} ELEMENTOF E}
            for (Vertex v : N) {					// for each v ELEMENTOF N
                if (isCEdge(u, v)) {		// if u and v are adjacent via a c-edge
                    if (T.contains(v)) {		// then if v ELEMENTOF T
                        //s.add(v);			// S <- S UNION {v}
                    } else {
                        P.add(v);			// else P <- P UNION {v}
                    }
                } else if (isDEdge(u, v)) {// else if u and v are adjacent via a d-edge
                    D.add(v);				// D <- D UNION {v}
                }
            }
            Set<Vertex> C = new LinkedHashSet<>();
            C.add(u);
            Set<Vertex> subresult = Enumerate_C_Cliques(graph, C, P, D, T, currentmaxresult); //ENUMERATE....
//            Set<Vertex> subresult = Enumerate_C_Cliques_Complex(graph, C, P, D, T, currentmaxresult); //ENUMERATE....
            if (subresult != null && subresult.size() > result.size()) {
                result = subresult;
                currentmaxresult = result.size();
                cliques.add(result);
            }
            T.add(u);						// T <- T UNION {v}
            if (DEBUG) {
                System.out.println("Current Max " + currentmaxresult);
            }
        }
        if (DEBUG) {
            System.out.println("cliques " + cliques);
        }
    }

    /**
     *
     * @param g The graph where the largest clique needs to be found
     * @param C Set of vertices belonging to the current clique
     * @param P Set of vertices which can be added to C, because they are
     * neighbours of vertex u via C-Edges
     * @param D Set of vertices which cannot directly be added to C because they
     * are neighbours of u via D-Edges
     * @return the largest clique in graph g
     */
    private Set<Vertex> Enumerate_C_Cliques(Graph comp_graph_nodes,
            Set<Vertex> C, Set<Vertex> P, Set<Vertex> D, Set<Vertex> T,
            int currentmaxresult) {
        Set<Vertex> result = new LinkedHashSet<>(C);

        if (iterations > 10000) {
            System.out.println("Reached max limit, 10000 itertions. ");
            return result;
        }
        iterations++;

        if (P.isEmpty() || P.size() + C.size() + D.size() <= currentmaxresult) { //if p=EMPTY and s=EMPTY
            cliques.add(result);
            return result;                               //REPORT.CLIQUE
        } else {
            Set<Vertex> P_Copy = new LinkedHashSet<>(P);
            for (Vertex ui : P) {                    			 //for i <- 1 to k
                P_Copy.remove(ui);             			 //P <-P\{ui}
                Set<Vertex> P_Prime = new LinkedHashSet<>(P_Copy);    //P' <- P
                Set<Vertex> D_Prime = new LinkedHashSet<>(D);		 //D' <- D
                Set<Vertex> N = findNeighbors(ui);//N <- { v ELEMENTOF Vertex | {ui,v} ELEMENTOF E }

                for (Vertex v : D) {							 // for all v ELEMENTOF D'
                    //(note that D and D' are the same at this point, to allow concurrent modification we loop over D)
                    if (isCEdge(ui, v)) {	 		 // if v and ui are adjacent via a c-edge
                        P_Prime.add(v);					 //P' <- P' UNION {v}
                        D_Prime.remove(v);				 //D' <- D'\{v}
                    }
                }

                Set<Vertex> C_Copy = new LinkedHashSet<>(C);
                C_Copy.add(ui);               			 //C UNION {ui}
                P_Prime.retainAll(N);                      //P' INTERSECTION N
                D_Prime.retainAll(N);						 //D' INTERSECTION N

                Set<Vertex> clique = Enumerate_C_Cliques(comp_graph_nodes, C_Copy, P_Prime, D_Prime, T, currentmaxresult); //ENUMERATE.C_CLIQUES....

                if (clique != null && clique.size() > result.size()) {
                    result = clique;
                    currentmaxresult = clique.size();
                    cliques.add(result);
                }
            }
        }
        return result;
    }

    /**
     *
     * @param g The graph where the largest clique needs to be found
     * @param C Set of vertices belonging to the current clique
     * @param P Set of vertices which can be added to C, because they are
     * neighbours of vertex u via C-Edges
     * @param D Set of vertices which cannot directly be added to C because they
     * are neighbours of u via D-Edges
     * @param s set of vertices which are not allowed to be added to C
     * @return the largest clique in graph g
     */
    private Set<Vertex> Enumerate_C_Cliques_Complex(Graph comp_graph_nodes,
            Set<Vertex> C, Set<Vertex> P, Set<Vertex> D, Set<Vertex> T, int currentmaxresult) {

        Set<Vertex> result = new LinkedHashSet<>(C);
        if (iterations > 10000) {
            System.out.println("Reached max limit, 10000 itertions. ");
            return result;
        }
        iterations++;

        if (P.isEmpty() || P.size() + C.size() + D.size() <= currentmaxresult) {//if P=EMPTY
            cliques.add(result);
            return result;                               //REPORT.CLIQUE
        } else {
            LinkedHashSet<Vertex> P_Copy = new LinkedHashSet<>(P);
            Vertex ut = P_Copy.iterator().next();
            for (Vertex ui : P) {                    			 //for i <- 1 to k
                Set<Vertex> target = new LinkedHashSet<>(D);
                target.removeAll(findNeighbors(ut)); //target is all vertices from D that are not adjacent to ut
                if (!comp_graph_nodes.hasEdge(ut, ui)
                        || // if ui is not adjacent to ut
                        hasCPath(ui, target, new LinkedHashSet<>())) { //or ui is connected via a C-path to a
                    //vertex form D that is not adjacent to ut

                    P_Copy.remove(ui);             			 //P <-P\{ui}
                    Set<Vertex> P_Prime = new LinkedHashSet<>(P_Copy);    //P' <- P
                    Set<Vertex> D_Prime = new LinkedHashSet<>(D);		 //D' <- D
                    Set<Vertex> N = findNeighbors(ui);//N <- { v ELEMENTOF Vertex | {ui,v} ELEMENTOF E }

                    for (Vertex v : D) {							 // for all v ELEMENTOF D'
                        //(note that D and D' are the same at this point, to allow concurrent modification we loop over D)
                        if (P.contains(v)) {					 //if v ELEMENTOF P
                            P_Prime.add(v);					 // then P' = P' UNION {v}
                        } else if (D.contains(v)) {			 // else if v ELEMENTOF D   \\can v be added to P?
                            if (isCEdge(ui, v)) {		 // then if v and ui are adjacent via a C-edge
                                //  \\is v an initializing vertex
                                if (!T.contains(v)) {			 // if v ELEMENTOF T
                                } else {
                                    P_Prime.add(v);			 // else P'=P' UNION {v}
                                }
                                D_Prime.remove(v);
                            }
                        }
                    }

                    Set<Vertex> C_Copy = new LinkedHashSet<>(C);
                    C_Copy.add(ui);               			 //C UNION {ui}
                    P_Prime.retainAll(N);                      //P' INTERSECTION N
                    D_Prime.retainAll(N);						 //D' INTERSECTION N
                    Set<Vertex> clique = Enumerate_C_Cliques_Complex(comp_graph_nodes, C_Copy, P_Prime, D_Prime, T, currentmaxresult); //ENUMERATE.C_CLIQUES....
                    if (clique.size() > result.size()) {
                        result = clique;
                        currentmaxresult = clique.size();
                        cliques.add(result);
                    }
                }
            }
        }
        return result;
    }

    private boolean isCEdge(Vertex u, Vertex v) {
        for (Edge e : cEdges) {
            if (e.getSource() == u && e.getSink() == v) {
                return true;
            }
            if (e.getSource() == v && e.getSink() == u) {
                return true;
            }
        }
        return false;
    }

    private boolean isDEdge(Vertex u, Vertex v) {
        for (Edge e : dEdges) {
            if (e.getSource() == u && e.getSink() == v) {
                return true;
            }
            if (e.getSource() == v && e.getSink() == u) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether an edge between vertices source and sink exists. whether
     * an edge exists between vertices x and y.
     *
     * @param x the index of the first Graph is in the edge list
     * @param y the index of the second Graph is in the edge list
     *
     * @param source the index of the first Graph is in the c-edge list
     * @param sink the index of the first Graph is in the c-edge list
     * @return true if a contact exists else false
     */
    private Set<Vertex> findNeighbors(Vertex central_node) {
        Set<Vertex> neighbors = new TreeSet<>();
        neighbors.addAll(this.graph.getNeighbours(central_node));

        if (DEBUG) {
            System.out.println("Node:" + central_node.getID() + " => Neighbors: " + neighbors);
        }
        return neighbors;
    }

    private boolean hasCPath(Vertex source, Set<Vertex> target, Set<Vertex> exclude) {
        //first check if there is a C_Edge from source to any element of target
        for (Vertex v : target) {
            if (isCEdge(source, v)) {
                return true;
            }
        }
        boolean result = false;
        //add source to the exclude list (no edge from source to any element from target exists)
        exclude.add(source);
        //check the same for every C_Neighbour of source
        Set<Vertex> neighbours = neighbourCVertices(source);
        //Remove all neighbours that have already been checked
        neighbours.removeAll(exclude);
        for (Vertex neighbour : neighbours) {
            //if there is a C-Path from a C-Neighbour of source to a vertex in target,
            //then there is a C_Path form source to the same vertex in target
            result |= hasCPath(neighbour, target, exclude);
        }
        return result;
    }

    /**
     * return only the vertices that neighbour u with C_Edges
     *
     * @param u a Vertex of g
     * @return {v ELEMENTOF Vertex | {u,v} ELEMENTOF E}
     */
    private Set<Vertex> neighbourCVertices(Vertex u) {
        Set<Vertex> neighbors = new LinkedHashSet<>();
        Set<Vertex> allNeighbours = this.graph.getNeighbours(u);
        for (Edge e : this.graph.edges()) {
//            System.out.println(" e " + e + ", EDGE Type " + e.getEdgeType());
            if (e.getEdgeType() != EdgeType.UNSET) {
                if (e.getSource().equals(u) && allNeighbours.contains(e.getSink())) {
                    if (allNeighbours.contains(e.getSink())) {
                        neighbors.add(e.getSink());
                    }
                }
                if (e.getSink().equals(u) && allNeighbours.contains(e.getSource())) {
                    if (allNeighbours.contains(e.getSource())) {
                        neighbors.add(e.getSource());
                    }
                }
            }
        }
        return neighbors;
    }
}
