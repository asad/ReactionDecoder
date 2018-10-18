/* Copyright (R) 2009-2018  Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus.mcsplus1;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.graphics.direct.MoleculeLabelDrawer;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class BKKCKCF {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MoleculeLabelDrawer.class);

    private final List<Integer> comp_graph_nodes;

    private final List<Integer> c_edges;

    private final List<Integer> d_edges;
    private final Stack<List<Integer>> max_Cliques_Set;
    private int best_clique_size;

    /*
     *T: is a set of vertices which have already been used for the
     * initialization of ENUMERATE_CLIQUES
     */
    protected final List<Integer> T;

    /**
     *
     * @param comp_graph_nodes
     * @param C_edges
     * @param D_edges
     */
    public BKKCKCF(List<Integer> comp_graph_nodes,
            List<Integer> C_edges,
            List<Integer> D_edges) {
        this.comp_graph_nodes = comp_graph_nodes;
        this.c_edges = C_edges;
        this.d_edges = D_edges;
        this.best_clique_size = 0;
        this.max_Cliques_Set = new Stack<>();
        this.T = new Stack<>();
    }

    /*
   
     * R: set of vertices belonging to the current clique
     
     * X: set of vertices which are not allowed to be added
     * to R, defined as X in paper
    
     * P: is a set of vertices which <b>can</b> be added to R, because they are
     * neighbours of vertex u via <i>c-edges</i>
    
     * Q: is a set of vertices which <b>cannot</b> be added to R, because they are
     * neighbours of vertex u via <i>d-edges</i>
     
     * V: stored all the vertices for the Graph G
     * V[G]: nodes of vector comp_graph_nodes are stored in V
     
     */
    int init_Algorithm() {

//        System.out.println( "init_Algorithm " + comp_graph_nodes.size() );
        List<Integer> R = new ArrayList<>();
        Stack<Integer> Q = new Stack<>();
        List<Integer> X = new ArrayList<>();
        List<Integer> N = new ArrayList<>();
        Stack<Integer> P = new Stack<>();

        //nodes of vector comp_graph_nodes are stored in V
        Stack<Integer> V = new Stack<>();//Initialization of Stack V

        int V_set_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < V_set_size; a++) {
            V.push(comp_graph_nodes.get(a * 3 + 2));
        }

        V.push(0);

        int b = 0;

        while (V.get(b) != 0) { // V[b] is node u
            int central_node = V.get(b);

            P.clear();
            Q.clear();
            X.clear();
            N.clear();
            R.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(central_node);

//              System.out.println("N-Neigh: " + N);
            for (int c = 0; c < N.size(); c = c + 2) { // N[c] is node v
                //System.out.println("N[" + c  +  "]= " + N.get(c) + " ");
                /*
                     * Grouping of the neighbors in X,P and Q
                     * u and v are adjacent via a R-edge
                 */
                if (N.get(c + 1) == 1) {
                    if (T.contains(N.get(c))) {
                        X.add(N.get(c));
                    } else {
                        P.push(N.get(c));
                    }

                } else if (N.get(c + 1) == 2) {   // u and v are adjacent via a Q-edge
                    //  System.out.println("u and v are adjacent via a Q-edge");
                    Q.push(N.get(c));
                }
                //find respective neighbor position in P, which is needed for the deletion from V
                int V_size = V.size();
                int neighbor_position = -1;

                //System.out.println("V Size: "+ V.size());
                for (int d = 0; d < V_size; d++) {
                    //System.out.println(" N[c]: " + N.get(c)+ " , V[" +  d + "]: " + V.get(d));
                    if (N.get(c).intValue() == (V.get(d))) {
                        neighbor_position = d;
                    }
                }
                //delete neighbor from set V
                if (neighbor_position != -1) {
                    //  System.out.println("neighbor_position : " + neighbor_position);
                    for (int e = neighbor_position; e < V_size - 1; e++) {
                        V.set(e, V.get(e + 1));
                    }
                    V.pop();
                    if (neighbor_position < b) {
                        b = b - 1;
                    }
                }
            }
            P.add(0);
            R.add(central_node);
            Enumerate_Cliques(R, P, Q, X);
            T.add(central_node);
            b++;
        }

        return 0;
    }

    private int Enumerate_Cliques(List<Integer> R, Stack<Integer> P, Stack<Integer> Q, List<Integer> X) {

        List<Integer> N = new ArrayList<>();////Initialization Vector N
        Stack<Integer> P_Prime = new Stack<>();//Defined as P' in the paper

        P.stream().forEach((I) -> {
            P_Prime.add(I);
        });

        List<Integer> R_copy = new ArrayList<>();
        Stack<Integer> P_copy = new Stack<>();
        Stack<Integer> Q_copy = new Stack<>();
        List<Integer> X_copy = new ArrayList<>();

        if (P.size() == 1) {
            if (X.isEmpty()) {

                //store best solutions in stack max_Cliques_Set
                int clique_size = R.size();
                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {
                        while (!max_Cliques_Set.empty()) {
                            getMax_Cliques_Set().pop();
                        }
                        best_clique_size = clique_size;
                    }
                    if (clique_size == best_clique_size) {
                        getMax_Cliques_Set().push(R);
                    }
                }
                // System.out.println("max_Cliques_Set: " + max_Cliques_Set.size());
                return 0;
            }
        }
        int a = 0;
        while (P_Prime.get(a) != 0) { // P[a] is node ut

            int ui = P_Prime.get(a);
            //remove P_Prime[a] from P
            //find position of P_Prime node in P
            int P_size = P.size();
            int ut_node_pos = Integer.MAX_VALUE;
            for (int counter = 0; counter < P_size - 1; counter++) {  //-1 wegen Endekennung
                if (P.get(counter).intValue() == P_Prime.get(a)) {
                    ut_node_pos = counter;
                }
            }
            if (ut_node_pos == Integer.MAX_VALUE) {
                LOGGER.debug("ut_node_pos = " + Integer.MAX_VALUE);
            }
            //delete P_Prime node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++) {
                P.setElementAt(P.get(counter + 1), counter);
            }
            P.pop();//TO DO

            R_copy.clear();
            P_copy.clear();
            Q_copy.clear();
            X_copy.clear();
            N.clear();

            R.stream().forEach((obj) -> {
                R_copy.add(obj);
            });

            P.stream().forEach((obj) -> {
                P_copy.add(obj);
            });
            Q.stream().forEach((obj) -> {
                Q_copy.add(obj);
            });
            X.stream().forEach((obj) -> {
                X_copy.add(obj);
            });
            P_copy.pop();

            //find the neighbors of the central node from P
            // System.out.println("P_Prime.get(a)" + P_Prime.get(a));
            N = find_neighbors(P_Prime.get(a));

            int N_size = N.size();

            for (int b = 0; b < N_size; b = b + 2) { // N[b] is node v
                int D_set_size = Q.size();
                int n_element_at_b = N.get(b);

                for (int c = 0; c < D_set_size; c++) {

                    if (n_element_at_b == Q.elementAt(c)) {
                        if (N.get(b + 1) == 1) {
                            //u and v are adjacent via a R-edge

                            if (T.contains(n_element_at_b)) {
                                X_copy.add(N.get(b));
                            } else {
                                P_copy.push(N.get(b));
                            }

                            int D_copy_size = Q_copy.size();
                            int n_b_position = Integer.MAX_VALUE;
                            for (int e = 0; e < D_copy_size; e++) {
                                if (n_element_at_b == Q_copy.elementAt(e)) {
                                    n_b_position = e;
                                }
                            }
                            for (int e = n_b_position; e < D_copy_size - 1; e++) {
                                Q_copy.set(e, Q_copy.get(e + 1));
                            }

                            Q_copy.pop();
                        }
                        /*//Abschnitt sinnlos, denn wenn etwas in X war ist, es nach X' kopiert worden
                         if(N[b+1] == 2){     //u and v are adjacent via a Q-edge
                         if().....
                         }*/
                    }
                }
                //find respective neighbor position in P_Prime, which is needed for the deletion from P_Prime
                int ut_set_size = P_Prime.size();
                int neighbor_position = -1;
                for (int e = 0; e < ut_set_size; e++) {
                    if (N.get(b).equals(P_Prime.get(e))) {
                        neighbor_position = e;
                    }
                }
                if (neighbor_position != -1) {
                    //delete neighbor from set P
                    for (int e = neighbor_position; e < ut_set_size - 1; e++) {
                        P_Prime.setElementAt(P_Prime.get(e + 1), e);
                    }
                    P_Prime.pop(); //TODO:Check whether size returns number of elements or index value
                    if (neighbor_position < a) {
                        a = a - 1;
                    }
                }
            }

            Stack<Integer> P_copy_N_intersec = new Stack<>();
            Stack<Integer> Q_copy_N_intersec = new Stack<>();
            List<Integer> X_copy_N_intersec = new ArrayList<>();

            for (int sec = 0; sec < N_size; sec += 2) {

                int nElement = N.get(sec);

                if (P_copy.contains(nElement)) {
                    P_copy_N_intersec.push(nElement);
                }
                if (Q_copy.contains(nElement)) {
                    Q_copy_N_intersec.push(nElement);
                }
                if (X_copy.contains(nElement)) {
                    X_copy_N_intersec.add(nElement);
                }

            }
            P_copy_N_intersec.push(0);
            R_copy.add(ui);
            Enumerate_Cliques(R_copy, P_copy_N_intersec, Q_copy_N_intersec, X_copy_N_intersec);
            X.add(ui);
            a++;
        }

        return 0;
    }

    private List<Integer> find_neighbors(int central_node) {

        List<Integer> neighbor_vec = new ArrayList<>();

        int C_edge_number = c_edges.size() / 2;
        for (int a = 0; a < C_edge_number; a++) {
            if (c_edges.get(a * 2 + 0) == central_node) {
                neighbor_vec.add(c_edges.get(a * 2 + 1));
                neighbor_vec.add(1);       // 1 means: is connected via R-edge
            }
            if (c_edges.get(a * 2 + 1) == central_node) {
                neighbor_vec.add(c_edges.get(a * 2 + 0));
                neighbor_vec.add(1);       // 1 means: is connected via R-edge
            }
        }

        int D_edge_number = d_edges.size() / 2;
        for (int a = 0; a < D_edge_number; a++) {
            if (d_edges.get(a * 2 + 0) == central_node) {
                neighbor_vec.add(d_edges.get(a * 2 + 1));
                neighbor_vec.add(2);       // 2 means: is connected via Q-edge
            }
            if (d_edges.get(a * 2 + 1) == central_node) {
                neighbor_vec.add(d_edges.get(a * 2 + 0));
                neighbor_vec.add(2);       // 2 means: is connected via Q-edge
            }
        }

        return neighbor_vec;
    }

    /**
     * @return the max_Cliques_Set
     */
    public Stack<List<Integer>> getMax_Cliques_Set() {
        return max_Cliques_Set;
    }

}
