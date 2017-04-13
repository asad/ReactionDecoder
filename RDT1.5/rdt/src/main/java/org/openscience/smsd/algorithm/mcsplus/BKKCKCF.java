/**
 *
 * Copyright (C) 2009-2015 Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received index copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, C. Karande: An Algorithm for reporting maximal
 * c-cliques; processedVertex.Comp. Sc. (2005); vol 349; pp. 484-490]
 *
 *
 * BronKerboschCazalsKarandeKochCliqueFinder.java
 *
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class BKKCKCF {

    private final Set<List<Integer>> max_Cliques_Set;
    /**
     * *****************************************************************
     */
    /*
     *T: is a set of vertices which have already been used for the
     * initialization of ENUMERATE_CLIQUES
     */
    private final List<Integer> T;
    /*
     * C: set of vertices belonging to the current clique
     */
    private final List<Integer> C;
    /*
     * S: set of vertices which are not allowed to be added
     * to C
     */
    private final List<Integer> S;
    /*
     *P: is a set of vertices which <b>can</b> be added to C, because they are
     * neighbours of vertex u via <i>c-edges</i>
     */
    private final Stack<Integer> P;
    /*
     * D: is a set of vertices which <b>cannot</b> be added to C, because they are
     * neighbours of vertex u via <i>d-edges</i>
     */
    private final Stack<Integer> D;
    /*
     *V: stored all the vertices for the Graph G
     * V[G]: nodes of vector comp_graph_nodes are stored in V
     */
    private final Stack<Integer> V;
    /**
     * ********************************************************************
     */
    private final List<Integer> C_edges;
    private final List<Integer> D_edges;
    private final List<Integer> comp_graph_nodes;

    private int best_clique_size;
    private List<Integer> C_copy;
    private Stack<Integer> P_copy;
    private Stack<Integer> D_copy;
    private List<Integer> S_copy;

    /**
     * Creates a new instance of BKKCKCF
     *
     * @param compGraphNodes
     * @param cEdges
     * @param dEdges
     */
    public BKKCKCF(
            List<Integer> compGraphNodes,
            List<Integer> cEdges,
            List<Integer> dEdges) {

        this.comp_graph_nodes = Collections.unmodifiableList(new ArrayList<>(compGraphNodes));
        this.C_edges = Collections.unmodifiableList(new ArrayList<>(cEdges));
        this.D_edges = Collections.unmodifiableList(new ArrayList<>(dEdges));
        best_clique_size = 0;
        max_Cliques_Set = new HashSet<>();

        T = new ArrayList<>(); //Initialize the T Vector
        C = new ArrayList<>();
        P = new Stack<>();
        D = new Stack<>();
        S = new ArrayList<>();
        V = new Stack<>();
        int V_set_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
        }
        V.add(0);

        /*
         * N[u]: set of neighbours of vertex u in Graph G
         *
         */
        List<Integer> N;
        int b = 0;

        /*
         * Let T be the set of Nodes already been used in the initialization
         *
         */
        T.clear();

        while (V.get(b) != 0) {

            int central_node = V.get(b);

            P.clear();
            D.clear();
            S.clear();
            C.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(V.get(b));
            //get neighbours and order S, P oder D
            for (int c = 0; c < N.size(); c += 2) {
                //Grouping of the neighbors in S,P and D

                /*
                 * u and v are adjacent via a C-edge
                 */
                if (N.get(c + 1) == 1) {
                    if (T.contains(N.get(c))) {
                        S.add(N.get(c));
                    } else {
                        P.push(N.get(c));
                    }

                } else if (N.get(c + 1) == 2) {
                    // u and v are adjacent via a D-edge
                    D.add(N.get(c));
                }
                //find respective neighbor position in P, which is needed for the deletion from V

                int neighbor_position = -1;

                int elementAtC = N.get(c);

                for (int d = 0; d < V.size(); d++) {
                    if (elementAtC == V.elementAt(d)) {
                        neighbor_position = d;
                    }
                }

                //delete neighbor from set V
                if (neighbor_position != -1) {
                    //System.out.println("neighbor_position : " + neighbor_position);
                    for (int e = neighbor_position; e < (V.size() - 1); e++) {
                        V.set(e, V.get(e + 1));
                    }
                    V.pop();
                    if (neighbor_position < b) {
                        b -= 1;
                    }
                }
            }
            P.add(0);
            C.add(central_node);
            Enumerate_Cliques(C, P, D, S);
            T.add(V.get(b));
            b++;
        }

    }

    private int Enumerate_Cliques(List<Integer> C, Stack<Integer> P, Stack<Integer> D, List<Integer> S) {

        List<Integer> N = new ArrayList<>(); ////Initialization Vector N
        Stack<Integer> P_Prime = new Stack<>();//Defined as P' in the paper

        C_copy = new ArrayList<>();
        P_copy = new Stack<>();
        D_copy = new Stack<>();
        S_copy = new ArrayList<>();

        for (Integer I : P) {
            P_Prime.add(I);
        }

        if (P.size() == 1) {
            if (S.isEmpty()) {
                //store best solutions in stack max_Cliques_Set
                int clique_size = C.size();

                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {
                        max_Cliques_Set.clear();
                        best_clique_size = clique_size;
//                        System.out.println("Best Cliques Size: " + best_clique_size + " " + clique_size);
                    }
                    if (clique_size == best_clique_size) {
                        max_Cliques_Set.add(C);
                    }
                }

                return 0;
            }
        }
        int a = 0;

        while (P_Prime.elementAt(a) != 0) {
            int ui = P_Prime.get(a);
            //remove P_Prime[a] from P
            //find position of P_Prime node in P
            int P_size = P.size();
            Integer ut_node_pos = 100000;
            for (int counter = 0; counter < P_size - 1; counter++) {
                if (P.elementAt(counter) == ui) {
                    ut_node_pos = counter;
                }
            }
            if (ut_node_pos == 100000) {
                System.out.println("ut_node_pos = 100000");
            }
            //delete P_Prime node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++) {
                P.set(counter, P.get(counter + 1));
            }

            P.pop();

            C_copy.clear();
            P_copy.clear();
            D_copy.clear();
            S_copy.clear();
            N.clear();

            for (Integer obj : C) {
                C_copy.add(obj);
            }

            for (Integer obj : P) {
                P_copy.add(obj);
            }
            for (Integer obj : D) {
                D_copy.add(obj);
            }
            for (Integer obj : S) {
                S_copy.add(obj);
            }

            P_copy.pop(); //Entferne Endekennung bei P_copy POP method() in CPP

            N = find_neighbors(P_Prime.get(a));

            int N_size = N.size();

            for (int b = 0; b < N_size; b += 2) {

                int D_set_size = D.size();
                int Nelement_at_b = N.get(b);

                for (int c = 0; c < D_set_size; c++) {

                    if (Nelement_at_b == D.elementAt(c)) {
                        if (N.get(b + 1) == 1) {
                            //u and v are adjacent via a C-edge

                            if (T.contains(Nelement_at_b)) {
                                S_copy.add(N.get(b));
                            } else {
                                P_copy.push(N.get(b));
                            }

                            int D_copy_size = D_copy.size();
                            int Nb_position = 10000;
                            for (int e = 0; e < D_copy_size; e++) {
                                if (Nelement_at_b == D_copy.elementAt(e)) {
                                    Nb_position = e;
                                }
                            }
                            for (int e = Nb_position; e < D_copy_size - 1; e++) {
                                D_copy.set(e, D_copy.get(e + 1));
                            }

                            D_copy.pop();
                        }
                        /*//Abschnitt sinnlos, denn wenn etwas in S war ist, es nach S' kopiert worden
                         if(N[b+1] == 2){     //u and v are adjacent via a D-edge
                         if().....
                         }*/
                    }
                }
                int ut_set_size = P_Prime.size();
                int neighbor_position = -1;
                for (int e = 0; e < ut_set_size; e++) {
                    if (Nelement_at_b == P_Prime.elementAt(e)) {
                        neighbor_position = e;
                    }
                }
                if (neighbor_position != -1) {
                    for (int e = neighbor_position; e < ut_set_size - 1; e++) {
                        P_Prime.set(e, P_Prime.get(e + 1));
                    }
                    P_Prime.pop(); //TODO:Check removeElementsAt to see whether size returns number of elements or index value
                    if (neighbor_position < a) {
                        a -= 1;
                    }
                }
            }
            Stack<Integer> P_copy_N_intersec = new Stack<>();
            Stack<Integer> D_copy_N_intersec = new Stack<>();
            List<Integer> S_copy_N_intersec = new ArrayList<>();

            int nElement;

            for (int sec = 0; sec < N_size; sec += 2) {

                nElement = N.get(sec);

                if (P_copy.contains(nElement)) {
                    P_copy_N_intersec.push(nElement);
                }
                if (D_copy.contains(nElement)) {
                    D_copy_N_intersec.add(nElement);
                }
                if (S_copy.contains(nElement)) {
                    S_copy_N_intersec.add(nElement);
                }

            }
            P_copy_N_intersec.add(0);
            C_copy.add(ui);

            Enumerate_Cliques(C_copy, P_copy_N_intersec, D_copy_N_intersec, S_copy_N_intersec);
            S.add(ui);
            a++;
        }
        return 0;
    }

    private List<Integer> find_neighbors(int central_node) {

        List<Integer> neighbor_vec = new ArrayList<>();

        int C_edge_number = C_edges.size() / 2;

        for (int a = 0; a < C_edge_number; a++) {
            if (C_edges.get(a * 2 + 0) == central_node) {
                neighbor_vec.add(C_edges.get(a * 2 + 1));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
            if (C_edges.get(a * 2 + 1) == central_node) {
                neighbor_vec.add(C_edges.get(a * 2 + 0));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
        }

        int D_edge_number = D_edges.size() / 2;
        for (int a = 0; a < D_edge_number; a++) {
            if (D_edges.get(a * 2 + 0) == central_node) {
                neighbor_vec.add(D_edges.get(a * 2 + 1));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }
            if (D_edges.get(a * 2 + 1) == central_node) {
                neighbor_vec.add(D_edges.get(a * 2 + 0));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }
        }

        return neighbor_vec;
    }

    public synchronized int getBestCliqueSize() {
        return best_clique_size;
    }

    /**
     *
     * @return
     */
    public synchronized Collection<List<Integer>> getMaxCliqueSet() {
        //System.out.println("max_Cliques_Set: " + max_Cliques_Set.size());
        return Collections.unmodifiableCollection(max_Cliques_Set);
    }
}
