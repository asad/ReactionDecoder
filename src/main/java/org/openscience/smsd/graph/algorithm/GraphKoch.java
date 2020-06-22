/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph.algorithm;

import java.util.ArrayList;
import org.openscience.smsd.graph.IClique;
import org.openscience.smsd.graph.Vertex;
import org.openscience.smsd.graph.Graph;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import org.openscience.smsd.graph.Edge;
import org.openscience.smsd.tools.IterationManager;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public class GraphKoch implements IClique {

    private final static boolean DEBUG = false;
    private final static boolean DEBUG2 = false;
    private final Collection<Set<Vertex>> cliques;
    private final Graph graph;
    IterationManager manager;

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
                    maxCliquesSet.push(new TreeSet<>(clique));
                }
            }
        }
//        System.out.println("maxCliquesSet " + maxCliquesSet);
        return maxCliquesSet;
    }

    /**
     *
     * @param compatibilityGraph
     */
    public GraphKoch(Graph compatibilityGraph) {
        this.graph = compatibilityGraph;
        this.cliques = new HashSet<>();
        int interation = this.graph.V() * 100;
        if (interation > 50000) {
            interation = 50000;
        }
        if (DEBUG2) {
            System.out.println("Edges: " + (this.graph.E()) + ", Vertex: " + (this.graph.V()));
            System.out.println(" GraphKoch Iterations " + interation);
        }
        this.manager = new IterationManager(interation);

        if (DEBUG) {
            System.out.println("GraphKoch");
        }
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
        if (DEBUG) {
            System.out.println("Starting koch ");
        }
        Set<Vertex> result = new LinkedHashSet<>();
        int currentmaxresult = 0;
        Set<Vertex> T = new LinkedHashSet<>();		// T <- Empty
        Set<Vertex> P, D, N, S;

        //set of vertices which have already been used for the initialization of Enumerate_C_Cliques()
        for (Vertex u : graph.nodes()) {				//for all u ELEMENTOF Vertex
            if (manager.isMaxIteration()) {
                //System.out.println("Reached max limit, " + manager.getIterationLimit() + " itertions. ");
                return;
            }

            P = new LinkedHashSet<>();			// P <- Empty
            D = new LinkedHashSet<>();			// D <- Empty
            S = new LinkedHashSet<>();			// S <- Empty
            N = findNeighbors(u);	// N <- {v ELEMENTOF Vertex | {u,v} ELEMENTOF E}
            if (DEBUG) {
                System.out.println("\n\n\nfindNeighbors = u => " + N.size());
            }
            for (Vertex v : N) {					// for each v ELEMENTOF N
                if (isCEdge(u, v)) {		// if u and v are adjacent via a c-edge
                    if (DEBUG) {
                        System.out.println("u " + u + ", v " + v);
                    }
                    if (T.contains(v)) {		// then if v ELEMENTOF T
                        S.add(v);			// S <- S UNION {v}
                    } else {
                        P.add(v);			// else P <- P UNION {v}
                    }
                } else if (isDEdge(u, v)) {// else if u and v are adjacent via a d-edge
                    D.add(v);				// D <- D UNION {v}
                }
            }
            Set<Vertex> C = new LinkedHashSet<>();
            C.add(u);

            if (DEBUG) {
                System.out.println("C " + C + ", P " + P + ", D " + D + ", T " + T);
            }
            Set<Vertex> subresult;
            subresult = Enumerate_C_Cliques(C, P, D, currentmaxresult); //ENUMERATE....(small footprint)
            if (subresult != null && subresult.size() >= result.size()) {
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
     * @param comp_graph_nodes The graph where the largest clique needs to be
     * found
     * @param C Set of vertices belonging to the current clique
     * @param P Set of vertices which can be added to C
     * @return the largest clique in graph comp_graph_nodes
     */
    private Set<Vertex> enumerateCliques(Set<Vertex> C, Set<Vertex> P, int currentmaxresult) {
        Set<Vertex> result = new LinkedHashSet<>(C);

        if (manager.isMaxIteration()) {
            //System.out.println("Reached max limit, " + manager.getIterationLimit() + " itertions. ");
            return result;
        }
        manager.increment();

        if (DEBUG2 && manager.getCounter() % 1000 == 0) {
            System.out.print("    Found clique #" + manager.getCounter()
                    + "/" + manager.getIterationLimit()
                    + " of size " + result.size() + ".\n");
        }

        if (P.isEmpty() || P.size() + C.size() <= currentmaxresult) { //if P=EMPTY
            ///or P.size() + C.size() <= currentmaxresult, 
            //if this is true, then the new clique can not be bigger then the clique 
            //that has already been reported/found
            return result;                                   //REPORT {CLIQUE}
        } else {
            List<Vertex> P_Copy = new ArrayList<>(P);
            Vertex ut = P_Copy.get(0);                           	 //Let ut be a vertex from P
            for (Vertex currentVertex : P) {                        //for i <- 1 to k
                if (!graph.hasEdge(ut, currentVertex)) {      //if ui is not adjacent ut ut
                    P_Copy.remove(currentVertex);             //P <-P\{ui}
                    Set<Vertex> P_Prime = new LinkedHashSet<>(P_Copy);    //P' <- P
                    Set<Vertex> N = new LinkedHashSet<>();
                    for (Edge edge : graph.edgesOf(currentVertex)) {
                        Vertex neighbour = graph.getEdgeSource(edge);
                        if (neighbour.equals(currentVertex)) {
                            neighbour = graph.getEdgeTarget(edge);
                        }
                        N.add(neighbour);
                    }                                        //N <- { v ELEMENTOF Vertex | {ui,v} ELEMENTOF Edge }

                    Set<Vertex> C_Copy = new LinkedHashSet<>(C);
                    C_Copy.add(currentVertex);                //C UNION {ui}
                    P_Prime.retainAll(N);                      //P' INTERSECTION N

                    Set<Vertex> clique = enumerateCliques(C_Copy, P_Prime, currentmaxresult); //ENUMERATE.CLIQUES....
                    if (clique.size() > result.size()) {
                        result = clique;
                        currentmaxresult = clique.size();
                    }
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
     * @return the largest clique in graph g
     */
    private Set<Vertex> Enumerate_C_Cliques(
            Set<Vertex> C, Set<Vertex> P, Set<Vertex> D,
            int currentmaxresult) {
        Set<Vertex> result = new LinkedHashSet<>(C);

        if (manager.isMaxIteration()) {
            //System.out.println("Reached max limit, " + manager.getIterationLimit() + " itertions. ");
            return result;
        }
        manager.increment();

        if (DEBUG2 && manager.getCounter() % 1000 == 0) {
            System.out.print("    Found clique #" + manager.getCounter()
                    + "/" + manager.getIterationLimit()
                    + " of size " + result.size() + ".\n");
        }

        if (P.isEmpty() || P.size() + C.size() + D.size() <= currentmaxresult) { //if p=EMPTY and s=EMPTY
            return result;                               //REPORT.CLIQUE
        } else {
            Set<Vertex> P_Copy = new LinkedHashSet<>(P);
            for (Vertex ui : P) {                    			 //for i <- 1 to k
                P_Copy.remove(ui);             			 //P <-P\{ui}
                Set<Vertex> P_Prime = new LinkedHashSet<>(P_Copy);    //P' <- P
                Set<Vertex> D_Prime = new LinkedHashSet<>(D);		 //D' <- D
                Set<Vertex> N = findNeighbors(ui);//N <- { v ELEMENTOF Vertex | {ui,v} ELEMENTOF E }
                if (DEBUG) {
                    System.out.println(ui + " = Nighbours: " + N.size());
                }
                D.stream().filter((v) -> (isCEdge(ui, v))).map((v) -> {
                    // if v and ui are adjacent via a c-edge
                    P_Prime.add(v);					 //P' <- P' UNION {v}
                    return v;
                }).forEachOrdered((v) -> {
                    D_Prime.remove(v);				 //D' <- D'\{v}
                }); // for all v ELEMENTOF D'
                //(note that D and D' are the same at this point, to allow concurrent modification we loop over D)

                Set<Vertex> C_Copy = new LinkedHashSet<>(C);
                C_Copy.add(ui);               			 //C UNION {ui}
                P_Prime.retainAll(N);                      //P' INTERSECTION N
                D_Prime.retainAll(N);						 //D' INTERSECTION N

                Set<Vertex> clique = Enumerate_C_Cliques(C_Copy, P_Prime, D_Prime, currentmaxresult); //ENUMERATE.C_CLIQUES....

                if (clique != null && clique.size() > result.size()) {
                    result = clique;
                    currentmaxresult = clique.size();
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
    private Set<Vertex> Enumerate_C_Cliques_Complex(
            Set<Vertex> C, Set<Vertex> P, Set<Vertex> D, Set<Vertex> T, int currentmaxresult) {

        Set<Vertex> result = new LinkedHashSet<>(C);
        if (manager.isMaxIteration()) {
            //System.out.println("Reached max limit," + manager.getIterationLimit() + " itertions. ");
            return result;
        }
        manager.increment();

        if (DEBUG2 && manager.getCounter() % 1000 == 0) {
            System.out.print("    Found clique #" + manager.getCounter()
                    + "/" + manager.getIterationLimit()
                    + " of size " + result.size() + ".\n");
        }

        if (P.isEmpty() || P.size() + C.size() + D.size() <= currentmaxresult) {//if P=EMPTY
            return result;                               //REPORT.CLIQUE
        } else {
            LinkedHashSet<Vertex> P_Copy = new LinkedHashSet<>(P);
            Vertex ut = P_Copy.iterator().next();
            for (Vertex ui : P) {                    			 //for i <- 1 to k
                Set<Vertex> target = new LinkedHashSet<>(D);
                //target.removeAll(findNeighbors(ut)); //target is all vertices from D that are not adjacent to ut
                target.removeAll(findNeighbors(ut)); //target is all vertices from D that are not adjacent to ut
                if (!graph.hasEdge(ut, ui)
                        || // if ui is not adjacent to ut
                        hasCPath(ui, target, new LinkedHashSet<>())) { //or ui is connected via a C-path to a
                    //vertex from D that is not adjacent to ut

                    P_Copy.remove(ui);             			 //P <-P\{ui}
                    Set<Vertex> P_Prime = new LinkedHashSet<>(P_Copy);    //P' <- P
                    Set<Vertex> D_Prime = new LinkedHashSet<>(D);		 //D' <- D
                    Set<Vertex> N = findNeighbors(ui);//N <- { v ELEMENTOF Vertex | {ui,v} ELEMENTOF E }
                    D.forEach((v) -> {
                        // for all v ELEMENTOF D'
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
                    });

                    Set<Vertex> C_Copy = new LinkedHashSet<>(C);
                    C_Copy.add(ui);               			 //C UNION {ui}
                    P_Prime.retainAll(N);                      //P' INTERSECTION N
                    D_Prime.retainAll(N);						 //D' INTERSECTION N
                    Set<Vertex> clique = Enumerate_C_Cliques_Complex(C_Copy, P_Prime, D_Prime, T, currentmaxresult); //ENUMERATE.C_CLIQUES....
                    if (clique.size() > result.size()) {
                        result = clique;
                        currentmaxresult = clique.size();
                    }
                }
            }
        }
        return result;
    }

    private boolean isCEdge(Vertex u, Vertex v) {
        return graph.isCEdge(u, v);
    }

    private boolean isDEdge(Vertex u, Vertex v) {
        return graph.isDEdge(u, v);
    }

    /**
     * Returns whether an edge between vertices source and sink exists. whether
     * an edge exists between vertices x and y.
     *
     * @param u a Vertex of g
     * @return {v ELEMENTOF V | {u,v} ELEMENTOF E}
     *
     * @return true if a contact exists else false
     */
    private Set<Vertex> findNeighbors(Vertex central_node) {
        Set<Vertex> allNeighbours = this.graph.getNeighbours(central_node);
        if (DEBUG) {
            System.out.println("Vertex:" + central_node.getID() + " => all Neighbours: " + allNeighbours);
        }
        return allNeighbours;
    }

    private boolean hasCPath(Vertex source, Set<Vertex> target, Set<Vertex> exclude) {
        //first check if there is a C_Edge from source to any element of target
        if (target.stream().anyMatch((v) -> (isCEdge(source, v)))) {
            return true;
        }

        boolean result = false;
        //add source to the exclude list (no edge from source to any element from target exists)
        exclude.add(source);
        //check the same for every C_Neighbour of source
        Set<Vertex> neighbours = neighbourCVertices(source);
        //Remove all neighbours that have already been checked
        neighbours.removeAll(exclude);
        result = neighbours.stream().map((neighbour)
                -> hasCPath(neighbour, target, exclude)).reduce(result, (accumulator, _item)
                -> accumulator | _item); //if there is a C-Path from a C-Neighbour of source to a vertex in target,
        //then there is a C_Path from source to the same vertex in target
        return result;
    }

    /**
     * return only the vertices that neighbour u with C_Edges
     *
     * @param u a Vertex of g
     * @return {v ELEMENTOF Vertex | {u,v} ELEMENTOF E}
     */
    private Set<Vertex> neighbourCVertices(Vertex u) {
        Set<Vertex> allNeighbours = this.graph.getCEdgeNeighbours(u);
        return allNeighbours;
    }
}
