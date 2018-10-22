/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

/**
 * This class implements Bron Kerbosch with pivot between query and target
 * molecule. It also marks edges in the compatibility graph as c-edges or
 * d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class GraphBronKerboschPivotal {

    /**
     * Finds the largest maximal cliques of the graph.
     *
     * @return the largest cliques
     */
    public Stack<Set<Node>> getMaxCliquesSet() {
        Stack<Set<Node>> maxCliquesSet = new Stack<>();
        int best_clique_size = 0;
        for (Set<Node> clique : cliques) {
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
     * @param cliques set the cliques
     */
    public void setCliques(Stack<Set<Node>> cliques) {
        this.cliques = cliques;
    }

    /**
     *
     * @return Collection of cliques (each of which is represented as a Set of
     * vertices)
     */
    public Stack<Set<Node>> getCliques() {
        return cliques;
    }

    private final static boolean DEBUG = false;
    private Stack<Set<Node>> cliques;

    private final List<Integer> comp_graph_nodes;
    private final List<Edge> c_edges;
    private final List<Edge> d_edges;

    private Integer iterations;

    public GraphBronKerboschPivotal(List<Integer> comp_graph_nodes,
            List<Edge> cEdges,
            List<Edge> dEdges) {
        this.comp_graph_nodes = comp_graph_nodes;
        this.c_edges = cEdges;
        this.d_edges = dEdges;
        this.cliques = new Stack<>();
        this.iterations = 0;

    }

    /**
     * Finds all maximal cliques of the graph.
     *
     *
     * A clique is maximal if it is impossible to extend it by adding another
     * vertex from the graph.
     *
     * Note: A maximal clique is not necessarily the largest clique in the
     * graph.
     *
     */
    public void findMaximalCliques() {

        Stack<Node> potential_clique_R = new Stack<>();//R, 
        List<Node> candidates_P = new ArrayList<>();//P
        List<Node> already_found_X = new ArrayList<>();//X

        // add all candidate vertices
        int candidates_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < candidates_size; a++) {
            candidates_P.add(new Node(comp_graph_nodes.get(a * 3 + 2)));
        }

        int printDepth = 1;
        BronKerboschWithPivot(potential_clique_R, candidates_P, already_found_X, printDepth);
        if (DEBUG) {
            System.out.println("BK Cliques Found: \n " + getCliques());
        }
    }

    /**
     *
     * With pivoting The basic form of the algorithm, described above, is
     * inefficient in the case of graphs with many non-maximal cliques: it makes
     * a recursive call for every clique, maximal or not. To save time and allow
     * the algorithm to backtrack more quickly in branches of the search that
     * contain no maximal cliques, Bron and Kerbosch introduced a variant of the
     * algorithm involving a "pivot vertex" u, chosen from P (or more generally,
     * as later investigators realized,[4] from P ⋃ X). Any maximal clique must
     * include either u or one of its non-neighbors, for otherwise the clique
     * could be augmented by adding u to it. Therefore, only u and its
     * non-neighbors need to be tested as the choices for the vertex v that is
     * added to R in each recursive call to the algorithm. In pseudocode:
     *
     * https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm#cite_note-1
     *
     * ``` BronKerboschPivoting(R,P,X):
     *
     * ----- if P and X are both empty:
     *
     * ---------- report R as a maximal clique
     *
     * ----- choose a pivot vertex u in P ⋃ X
     *
     * ----- for each vertex v in P \ N(u):
     *
     * ---------- BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
     *
     * ---------- P := P \ {v}
     *
     * ---------- X := X ⋃ {v}'
     *
     * ```
     *
     * If the pivot is chosen to minimize the number of recursive calls made by
     * the algorithm, the savings in running time compared to the non-pivoting
     * version of the algorithm can be significant.[5]
     *
     *
     * R := is the set of nodes of a maximal clique. <potential clique>
     *
     * P := is the set of possible nodes in a maximal clique. <candidates>
     *
     * X := is the set of nodes that are excluded. <already found>
     *
     * @param R R := is the set of nodes of a maximal clique. (potential clique)
     * @param P P := is the set of possible nodes in a maximal clique
     * (candidates)
     * @param X X := is the set of nodes that are excluded. (already found)
     * @param printDepth
     */
    private void BronKerboschWithPivot(
            Stack<Node> R,
            List<Node> P,
            List<Node> X,
            int printDepth) {

        if (DEBUG) {
            System.out.println("BronKerboschWithPivot called: C=" + toText(R, "{", "}")
                    + ", P=" + toText(P, "{", "}") + ", S=" + toText(X, "{", "}"));
        }

        if ((P.isEmpty()) && (X.isEmpty())) {
            Set<Node> pushed = cliques.push(new HashSet<>(R));
            //printClique(R);
            return;
        }
        if (iterations > 5000) {
            System.out.println("Reached max limit, 5000 itertions. ");
            return;
        }

        // System.out.println(); 
        List<Node> P1 = new ArrayList<>(P);
        // Find Pivot 
        Node u = getMaxDegreeVertex(union(P, X));

        iterations++;
        if (this.iterations % 1000 == 0) {
            if (DEBUG) {
                System.out.print("    Found clique #" + this.iterations + " of size " + R.size() + ".\n");
            }
            if (DEBUG) {
                printClique(R);
            }
        }
        if (DEBUG) {
            System.out.println("" + printDepth + " Pivot is " + (u));
        }
        //P = P / Nbrs(u) 
        P = removeNeigbour(P, u);

        for (Node v : P) {
            R.push(v);
            BronKerboschWithPivot(R, intersect(P1, getNeigbours(v)),
                    intersect(X, getNeigbours(v)), printDepth + 1);
            R.pop();//remove v
            P1.remove(v);
            X.add(v);
        }
    }

    /*
     * Returns max degree of a vertex
     */
    private Node getMaxDegreeVertex(List<Node> t) {
        int i = 0, temp = 0;
        Node n = null;
        while (i < t.size()) {
            if (getDegreeVertex(t.get(i)) > temp) {
                temp = getDegreeVertex(t.get(i));
                n = t.get(i);
            }
            i += 1;
        }
        return n;
    }

    /*
     * Returns degree of a vertex
     */
    private int getDegreeVertex(Node node) {
        return findNeighbors(node).size();
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
    private Collection<Node> findNeighbors(Node central_node) {

        Set<Node> neighbors = new TreeSet<>();
        /*
         * Add all the valid nodes via c-edges
         */

        for (Edge e : c_edges) {
            if (e.getSource() == central_node.node) {
                neighbors.add(new Node(e.getSink()));
            }
            if (e.getSink() == central_node.node) {
                neighbors.add(new Node(e.getSource()));
            }
        }

        if (DEBUG) {
            System.out.println("neighbors c-edges " + neighbors);
        }

        /*
         * remove all the nodes via d-edges
         */
        for (Edge e : d_edges) {
            if (e.getSource() == central_node.node) {
                neighbors.remove(new Node(e.getSink()));
            }
            if (e.getSink() == central_node.node) {
                neighbors.remove(new Node(e.getSource()));
            }
        }

        if (DEBUG) {
            System.out.println("neighbors d-edges " + neighbors);
        }
        return neighbors;
    }

    private List<Node> getNeigbours(Node v) {
        List<Node> t = new ArrayList<>();
        t.addAll(findNeighbors(v));
        return t;
    }

    // Intersection of two sets 
    private List<Node> intersect(List<Node> source, List<Node> sink) {
        List<Node> intersection = new ArrayList<>(source);
        intersection.retainAll(sink);
        return intersection;
    }

    // Union of two sets 
    private List<Node> union(List<Node> source, List<Node> sink) {
        List<Node> union = new ArrayList<>(source);
        union.addAll(sink);
        return union;
    }

    // Removes the neigbours 
    private List<Node> removeNeigbour(List<Node> source, Node v) {
        List<Node> remaining = new ArrayList<>(source);
        remaining.removeAll(getNeigbours(v));
        return remaining;
    }

    /**
     * Debug function to get a string representation of a list.
     *
     * @param solution the list
     * @param start
     * @param end
     * @return the string
     */
    private static String toText(List<Node> solution, String start, String end) {
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        solution.forEach((i) -> {
            sb.append(i).append(" ");
        });
        sb.append(end);
        return sb.toString();
    }

    /**
     *
     * @param s
     * @return
     */
    private static String addSpacer(int s) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < s; i++) {
            sb.append("    ");
        }
        return sb.toString();
    }

    private void printClique(List<Node> R) {
        System.out.print("Clique Set R=[");
        R.forEach((v) -> {
            System.out.print(" " + v);
        });
        System.out.print(" ]\n");
    }
}
