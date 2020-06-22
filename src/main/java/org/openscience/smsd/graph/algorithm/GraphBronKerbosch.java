/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph.algorithm;

import org.openscience.smsd.graph.IClique;
import org.openscience.smsd.graph.Vertex;
import org.openscience.smsd.graph.Graph;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import org.openscience.smsd.tools.IterationManager;

/**
 * This class implements Bron Kerbosch with pivot between query and target
 * molecule. It also marks edges in the compatibility graph as c-edges or
 * d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman @ bioinceptionlabs.com>
 */
public class GraphBronKerbosch implements IClique {

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
     * @return Collection of cliques (each of which is represented as a Set of
     * vertices)
     */
    @Override
    public Collection<Set<Vertex>> getCliques() {
        return cliques;
    }

    private final static boolean DEBUG = false;
    private final Collection<Set<Vertex>> cliques;
    IterationManager manager;
    private final Graph graph;

    public GraphBronKerbosch(
            Graph comp_graph_nodes) {
        this.graph = comp_graph_nodes;
        this.cliques = new HashSet<>();
        int interation = 2 * this.graph.V();
        if (interation > 1000) {
            interation = 1000;
        }
        this.manager = new IterationManager(interation);
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
    @Override
    public void findMaximalCliques() {

        TreeSet<Vertex> potential_clique_R = new TreeSet<>();//R, 
        TreeSet<Vertex> candidates_P = new TreeSet<>();//P
        TreeSet<Vertex> already_found_X = new TreeSet<>();//X
        // add all candidate vertices
        Iterator<Vertex> iterator = graph.iterator();
        while (iterator.hasNext()) {
            candidates_P.add(iterator.next());
        }

        int printDepth = 1;

        BronKerboschWithPivot(potential_clique_R, candidates_P, already_found_X, printDepth);
//        BronKerboschWithoutPivot(potential_clique_R, candidates_P, already_found_X, printDepth);
//        BronKerbosch(potential_clique_R, candidates_P, already_found_X);
        if (DEBUG) {
            System.out.println("BK Cliques Found: \n " + this.cliques);
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
            TreeSet<Vertex> R,
            TreeSet<Vertex> P,
            TreeSet<Vertex> X,
            int printDepth) {

        if (DEBUG) {
            System.out.println("BronKerboschWithPivot called: R=" + toText(R, "{", "}")
                    + ", P=" + toText(P, "{", "}") + ", X=" + toText(X, "{", "}"));
        }

        if ((P.isEmpty()) && (X.isEmpty())) {
            cliques.add(new HashSet<>(R));
            if (DEBUG) {
                printClique(R);
            }
            return;
        }
        if (manager.isMaxIteration()) {
            //System.out.println("Reached max limit, " + manager.getIterationLimit() + " itertions. ");
            return;
        }
        manager.increment();

        if (DEBUG && manager.getCounter() % 1000 == 0) {
            System.out.print("    Found clique #" + manager.getCounter() + " of size " + R.size() + ".\n");
            printClique(R);
        }

        Set<Vertex> P1 = new TreeSet<>(P);
        if (DEBUG) {
            System.out.println("P_PRIME " + P1);
        }

        /*
         * Find Pivot 
         */
        Vertex u = getMaxDegreeVertex(new ArrayList<>(union(P1, X)));
        /*
         * P = P / Nbrs(u) 
         */
        P1 = new TreeSet<>(removeNeigbour(P1, u));

        if (DEBUG) {
            System.out.println("P_Prime: " + P1 + " Depth: " + printDepth + " Pivot is " + (u));
        }
        for (Vertex v : P1) {
            //Push the id into selection set
            R.add(v);
            //Find neighbours
            List<Vertex> neighbors = new ArrayList<>(findNeighbors(v));
            if (DEBUG) {
                System.out.println("Neighbours of v " + v + " are " + neighbors);
            }
            BronKerboschWithPivot(R, new TreeSet<>(intersect(P, neighbors)),
                    new TreeSet<>(intersect(X, neighbors)), printDepth + 1);
            R.remove(v);
            P.remove(v);
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
            TreeSet<Vertex> R,
            TreeSet<Vertex> P,
            TreeSet<Vertex> X,
            int printDepth) {

        if (DEBUG) {
            System.out.println("BronKerboschWithPivot called: R=" + toText(R, "{", "}")
                    + ", P=" + toText(P, "{", "}") + ", X=" + toText(X, "{", "}"));
        }

        if ((P.isEmpty()) && (X.isEmpty())) {
            cliques.add(new HashSet<>(R));
            if (DEBUG) {
                printClique(R);
            }
            return;
        }

        /*
         * Find Pivot 
         */
        Vertex v = null;
        if (!P.isEmpty()) {
            v = P.first();
        }
        while (!P.isEmpty() && v != P.last()) {
            R.add(v);
            //Find neighbours
            List<Vertex> neighbors = new ArrayList<>(findNeighbors(v));
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
    private Vertex getMaxDegreeVertex(List<Vertex> t) {
        int i = 0, temp = 0;
        Vertex n = null;
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

    private void BronKerbosch(
            TreeSet<Vertex> R,
            TreeSet<Vertex> P,
            TreeSet<Vertex> X) {
        TreeSet<Vertex> candidates_array = new TreeSet<>(P);
        if (!end(P, X)) {
            // for each candidate_node in P do
            candidates_array.stream().map((candidate) -> {
                TreeSet<Vertex> new_candidates = new TreeSet<>();
                TreeSet<Vertex> new_already_found = new TreeSet<>();
                // move candidate id to R
                R.add(candidate);
                P.remove(candidate);
                // create new_candidates by removing nodes in P not
                // connected to candidate id
                P.stream().filter((new_candidate) -> (isNeighbor(candidate, new_candidate))).forEachOrdered((new_candidate) -> {
                    new_candidates.add(new_candidate);
                }); // of if
                // of for
                // create new_already_found by removing nodes in X
                // not connected to candidate id
                X.stream().filter((new_found) -> (isNeighbor(candidate, new_found))).forEachOrdered((new_found) -> {
                    new_already_found.add(new_found);
                }); // of if
                // of for
                // if new_candidates and new_already_found are empty
                if (new_candidates.isEmpty() && new_already_found.isEmpty()) {
                    // R is maximal_clique
                    cliques.add(new HashSet<>(R));
                } // of if
                else {
                    // recursive call
                    BronKerbosch(R,
                            new_candidates,
                            new_already_found);
                } // of else
                // move candidate_node from R to X;
                X.add(candidate);
                return candidate;
            }).forEachOrdered((candidate) -> {
                R.remove(candidate);
            }); // of for
        } // of if
    }

    /**
     *
     * @param candidates
     * @param already_found
     * @return
     */
    private boolean end(TreeSet<Vertex> candidates, TreeSet<Vertex> already_found) {
        // if a id in already_found is connected to all nodes in candidates
        boolean end = false;
        int edgecounter;
        for (Vertex found : already_found) {
            edgecounter = 0;
            edgecounter = candidates.stream().filter((candidate)
                    -> (isNeighbor(found, candidate))).map((_item) -> 1).
                    reduce(edgecounter, Integer::sum); // of if
            // of for
            if (edgecounter == candidates.size()) {
                end = true;
            }
        } // of for
        return end;
    }

    /*
     * Returns degree of a vertex
     */
    private int getDegreeVertex(Vertex node) {
        return this.graph.getDegree(node);
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
    // Intersection of two sets 

    private Collection<Vertex> intersect(Collection<Vertex> source, Collection<Vertex> sink) {
        Set<Vertex> intersection = new HashSet<>(source);
        intersection.retainAll(sink);
        return intersection;
    }

    // Union of two sets 
    private Collection<Vertex> union(Collection<Vertex> source, Collection<Vertex> sink) {
        Set<Vertex> union = new HashSet<>(source);
        union.addAll(sink);
        return union;
    }

    // Removes the neigbours 
    private Collection<Vertex> removeNeigbour(Collection<Vertex> source, Vertex v) {
        Set<Vertex> remaining = new HashSet<>(source);
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
    private static String toText(Collection<Vertex> solution, String start, String end) {
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

    private void printClique(Collection<Vertex> R) {
        System.out.print("Clique Set R=[");
        R.forEach((v) -> {
            System.out.print(" " + v + ",");
        });
        System.out.print(" ]\n");
    }

    private boolean isNeighbor(Vertex found, Vertex candidate) {
        return graph.hasEdge(found, candidate);
    }
}
