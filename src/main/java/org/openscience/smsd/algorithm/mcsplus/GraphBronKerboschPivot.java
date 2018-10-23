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
public class GraphBronKerboschPivot {

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

    private final static boolean DEBUG = true;
    private Stack<Set<Node>> cliques;

    private final List<Integer> comp_graph_nodes;
    private final List<Edge> c_edges;
    private final List<Edge> d_edges;

    private Integer iterations;

    public GraphBronKerboschPivot(List<Integer> comp_graph_nodes,
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

        TreeSet<Node> potential_clique_R = new TreeSet<>();//R, 
        TreeSet<Node> candidates_P = new TreeSet<>();//P
        TreeSet<Node> already_found_X = new TreeSet<>();//X

        // add all candidate vertices
        int candidates_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < candidates_size; a++) {
            candidates_P.add(new Node(comp_graph_nodes.get(a * 3 + 2)));
        }

        int printDepth = 1;

        BronKerboschWithPivot(potential_clique_R, candidates_P, already_found_X, printDepth);
//        BronKerboschWithoutPivot(potential_clique_R, candidates_P, already_found_X, printDepth);
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
            TreeSet<Node> R,
            TreeSet<Node> P,
            TreeSet<Node> X,
            int printDepth) {

        if (DEBUG) {
            System.out.println("BronKerboschWithPivot called: R=" + toText(R, "{", "}")
                    + ", P=" + toText(P, "{", "}") + ", X=" + toText(X, "{", "}"));
        }

        if ((P.isEmpty()) && (X.isEmpty())) {
            Set<Node> pushed = cliques.push(new HashSet<>(R));
            if (DEBUG) {
                printClique(R);
            }
            return;
        }
        if (iterations > 5000) {
            System.out.println("Reached max limit, 5000 itertions. ");
            return;
        }

        List<Node> P1 = new ArrayList<>(P);

        iterations++;
        if (this.iterations % 1000 == 0) {
            if (DEBUG) {
                System.out.print("    Found clique #" + this.iterations + " of size " + R.size() + ".\n");
            }
            if (DEBUG) {
                printClique(R);
            }
        }

        /*
         * Find Pivot 
         */
        Node u = getMaxDegreeVertex(new ArrayList<>(union(P, X)));
        /*
         * P = P / Nbrs(u) 
         */
        P = new TreeSet<>(removeNeigbour(P, u));
        if (DEBUG) {
            System.out.println("P_Prime: " + P + printDepth + " Pivot is " + (u));
        }
        for (Node v : P) {
            //Push the node into selection set
            R.add(v);
            //Find neighbours
            List<Node> neighbors = new ArrayList<>(findNeighbors(v));
            if (DEBUG) {
                System.out.println("Neighbours of v " + v + " are " + neighbors);
            }
            BronKerboschWithPivot(R, new TreeSet<>(intersect(P1, neighbors)),
                    new TreeSet<>(intersect(X, neighbors)), printDepth + 1);
            R.remove(v);
            P1.remove(v);
            X.add(v);
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
    private void BronKerboschWithoutPivot(
            TreeSet<Node> R,
            TreeSet<Node> P,
            TreeSet<Node> X,
            int printDepth) {

        if (DEBUG) {
            System.out.println("BronKerboschWithPivot called: R=" + toText(R, "{", "}")
                    + ", P=" + toText(P, "{", "}") + ", X=" + toText(X, "{", "}"));
        }

        if ((P.isEmpty()) && (X.isEmpty())) {
            Set<Node> pushed = cliques.push(new HashSet<>(R));
            if (DEBUG) {
                printClique(R);
            }
            return;
        }

        /*
         * Find Pivot 
         */
        Node v = null;
        if (!P.isEmpty()) {
            v = P.first();
        }
        while (!P.isEmpty() && v != P.last()) {
            R.add(v);
            //Find neighbours
            List<Node> neighbors = new ArrayList<>(findNeighbors(v));
            BronKerboschWithoutPivot(R, new TreeSet<>(intersect(P, neighbors)),
                    new TreeSet<>(intersect(X, neighbors)), printDepth + 1);
            P.remove(v);
            X.add(v);
            if (!P.isEmpty()) {
                v = P.first();
            }
        }
    }

    /*
     * Returns max degree of a vertex
     */
    private Node getMaxDegreeVertex(List<Node> t) {
        int i = 0, temp = 0;
        Node n = null;
        while (i < t.size()) {
            int degreeVertex = getDegreeVertex(t.get(i));
            if (degreeVertex > temp) {
                temp = degreeVertex;
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
        c_edges.stream().map((e) -> {
            if (e.getSource() == central_node.node) {
                neighbors.add(new Node(e.getSink()));
            }
            return e;
        }).filter((e) -> (e.getSink() == central_node.node)).forEachOrdered((e) -> {
            neighbors.add(new Node(e.getSource()));
        });

//        if (DEBUG) {
//            System.out.println(central_node.node + "=> neighbors c-edges " + neighbors);
//        }

        /*
         * remove all the nodes via d-edges
         */
        d_edges.stream().map((e) -> {
            if (e.getSource() == central_node.node) {
                neighbors.remove(new Node(e.getSink()));
            }
            return e;
        }).filter((e) -> (e.getSink() == central_node.node)).forEachOrdered((e) -> {
            neighbors.remove(new Node(e.getSource()));
        });

//        if (DEBUG) {
//            System.out.println(central_node.node + "=> neighbors d-edges " + neighbors);
//        }
        if (DEBUG) {
            System.out.println("Node:" + central_node.node + " => Neighbors: " + neighbors);
        }
        return neighbors;
    }

    // Intersection of two sets 
    private Collection<Node> intersect(Collection<Node> source, Collection<Node> sink) {
        Set<Node> intersection = new HashSet<>(source);
        intersection.retainAll(sink);
        return intersection;
    }

    // Union of two sets 
    private Collection<Node> union(Collection<Node> source, Collection<Node> sink) {
        Set<Node> union = new HashSet<>(source);
        union.addAll(sink);
        return union;
    }

    // Removes the neigbours 
    private Collection<Node> removeNeigbour(Collection<Node> source, Node v) {
        Set<Node> remaining = new HashSet<>(source);
        remaining.removeAll(findNeighbors(v));
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
    private static String toText(Collection<Node> solution, String start, String end) {
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        solution.forEach((i) -> {
            sb.append(i).append(",");
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

    private void printClique(Collection<Node> R) {
        System.out.print("Clique Set R=[");
        R.forEach((v) -> {
            System.out.print(" " + v + ",");
        });
        System.out.print(" ]\n");
    }
}
