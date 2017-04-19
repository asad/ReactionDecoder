/* Copyright (C) 2009-2017  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus1;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;

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
public class MCSPlus extends Filter {

    final Map<String, Integer> SYMBOL_VALUE = new TreeMap<>();

    /**
     * Creates a new instance of SearchCliques
     *
     *
     * @param f1
     * @param f2
     */
    public MCSPlus(MoleculeHandler f1, MoleculeHandler f2) {

        super(f1, f2);

    }

    private List<List<Integer>> label_atoms(List<Integer> basic_atom_vector, int bond_num, List<IAtom> atomstr, List<Integer> i_tab, List<String> c_tab) {

        ArrayList<List<Integer>> label_list = new ArrayList<>();

//        System.out.println("Vector Atom Str: " + atomstr);
//        System.err.println();
//        System.err.println("basic_atom_vector");
//        for (int b = 0; b < basic_atom_vector.size(); b++) {
//            System.err.print(basic_atom_vector.get(b));
//        }
//        System.err.println();
//        System.err.println("i_tab");
//        for (int b = 0; b < i_tab.size(); b++) {
//            System.err.print(i_tab.get(b));
//        }
//        System.err.println();
        for (int a = 0; a < basic_atom_vector.size(); a++) {

            List<Integer> label = new ArrayList<>(7);
            /*
             * Initialize the vector
             */
            for (int i = 0; i < 7; i++) {
                label.add(0);
            }

            String atom1_type = atomstr.get(a).getSymbol();
            if (SYMBOL_VALUE.containsKey(atom1_type)) {
                label.set(0, SYMBOL_VALUE.get(atom1_type));
            } else {
                int index = SYMBOL_VALUE.size() + 1;
                SYMBOL_VALUE.put(atom1_type, index);
                label.set(0, SYMBOL_VALUE.get(atom1_type));
            }
            int count_neighbors = 1;
            for (int b = 0; b < bond_num; b++) {
                if (basic_atom_vector.get(a).equals(i_tab.get(b * 3 + 0))) {
                    String atom2_type = c_tab.get(b * 2 + 1);
                    if (SYMBOL_VALUE.containsKey(atom2_type)) {
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    } else {
                        int index = SYMBOL_VALUE.size() + 1;
                        SYMBOL_VALUE.put(atom2_type, index);
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    }
                    count_neighbors++;
                }

                if (basic_atom_vector.get(a).equals(i_tab.get(b * 3 + 1))) {
                    String atom2_type = c_tab.get(b * 2 + 0);
                    if (SYMBOL_VALUE.containsKey(atom2_type)) {
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    } else {
                        int index = SYMBOL_VALUE.size() + 1;
                        SYMBOL_VALUE.put(atom2_type, index);
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    }
                    count_neighbors++;
                }
            }
//            System.out.println("SYMBOL_VALUE " + SYMBOL_VALUE);
//            System.out.println("label " + label);
            List<Integer> bubbleSort = Utility.getBubbleSort(label);
            label_list.add(bubbleSort);

        }

        return label_list;
    }

    private List<Integer> reduce_atomset(
            int atom_num,
            int bond_numb,
            List<IAtom> a_str,
            List<Integer> i_table,
            List<String> c_table) {

        List<Integer> phosphate_O_atoms = new ArrayList<>();
        List<Integer> h_atoms = new ArrayList<>();

        for (int a = 0; a < atom_num; a++) {
            if ("O".equals(a_str.get(a).getSymbol())) {
                int O_neighbor_num = 0;
                boolean P_neighbor = false;

                for (int b = 0; b < bond_numb; b++) {
                    if (a + 1 == i_table.get(b * 3 + 0)) {
                        O_neighbor_num++;
                        if (("P".equals(a_str.get(i_table.get(b * 3 + 1) - 1).getSymbol())) && (i_table.get(b * 3 + 2) != 2)) {
                            P_neighbor = true;
                        }
                    }
                    if (a + 1 == i_table.get(b * 3 + 1)) {
                        O_neighbor_num++;
                        if (("P".equals(a_str.get(i_table.get(b * 3 + 0) - 1).getSymbol())) && (i_table.get(b * 3 + 2) != 2)) {
                            P_neighbor = true;
                        }
                    }
                }
                if ((O_neighbor_num == 1) && (P_neighbor)) {
                    phosphate_O_atoms.add(a + 1);
                }
            }
            if ("H".equals(a_str.get(a).getSymbol())) {
                h_atoms.add(a + 1);
            }
        }

        List<Integer> basic_atoms = new ArrayList<>();
        int phosphate_O_atoms_size = phosphate_O_atoms.size();
        int H_atoms_size = h_atoms.size();

        for (int a = 0; a < atom_num; a++) {
            boolean no_P_O_atom = true;
            for (int b = 0; b < phosphate_O_atoms_size; b++) {
                if (a + 1 == phosphate_O_atoms.get(b)) {
                    no_P_O_atom = false;
                }
            }

            boolean no_H_atom = true;
            for (int b = 0; b < H_atoms_size; b++) {
                if (a + 1 == h_atoms.get(b)) {
                    no_H_atom = false;
                }
            }

            if ((no_P_O_atom) && (no_H_atom)) {
                basic_atoms.add(a + 1);
            }
        }
        return basic_atoms;
    }

    private int generate_compatibility_graph_nodes() {

        List<Integer> basic_atom_vec_A = reduce_atomset(atom_number1, bond_number1, atomstr1, i_tab1, c_tab1);
        List<Integer> basic_atom_vec_B = reduce_atomset(atom_number2, bond_number2, atomstr2, i_tab2, c_tab2);

        List<List<Integer>> label_list_molA = label_atoms(basic_atom_vec_A, bond_number1, atomstr1, i_tab1, c_tab1);
        List<List<Integer>> label_list_molB = label_atoms(basic_atom_vec_B, bond_number2, atomstr2, i_tab2, c_tab2);

        int molA_nodes = 0;
        int count_nodes = 1;

        for (List<Integer> labelA : label_list_molA) {
            int molB_nodes = 0;

            for (List<Integer> labelB : label_list_molB) {
                if (labelA.equals(labelB)) {
                    comp_graph_nodes.add(basic_atom_vec_A.get(molA_nodes));
                    comp_graph_nodes.add(basic_atom_vec_B.get(molB_nodes));
                    comp_graph_nodes.add(count_nodes++);
                }
                molB_nodes++;
            }
            molA_nodes++;
        }
        return 0;
    }

    private int generate_compatibility_graph() {

        int vector_size = comp_graph_nodes.size();

        for (int a = 0; a < vector_size; a = a + 3) {
            for (int b = a + 3; b < vector_size; b = b + 3) {
                if ((a != b) && (!comp_graph_nodes.get(a).equals(comp_graph_nodes.get(b)))
                        && (!comp_graph_nodes.get(a + 1).equals(comp_graph_nodes.get(b + 1)))) {
                    boolean molecule1_pair_connected = false;
                    boolean molecule2_pair_connected = false;

                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    for (int x = 0; x < bond_number1; x++) {
                        if ((comp_graph_nodes.get(a).equals(i_tab1.get(x * 3 + 0))
                                && comp_graph_nodes.get(b).equals(i_tab1.get(x * 3 + 1)))
                                || (comp_graph_nodes.get(a).equals(i_tab1.get(x * 3 + 1))
                                && comp_graph_nodes.get(b).equals(i_tab1.get(x * 3 + 0)))) {
                            molecule1_pair_connected = true;
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < bond_number2; y++) {
                        if ((comp_graph_nodes.get(a + 1).equals(i_tab2.get(y * 3 + 0))
                                && comp_graph_nodes.get(b + 1).equals(i_tab2.get(y * 3 + 1)))
                                || (comp_graph_nodes.get(a + 1).equals(i_tab2.get(y * 3 + 1))
                                && comp_graph_nodes.get(b + 1).equals(i_tab2.get(y * 3 + 0)))) {
                            molecule2_pair_connected = true;
                        }
                    }
                    //in case that both molecule pairs are connected a c-edge is generated
                    if (molecule1_pair_connected && molecule2_pair_connected) {

                        C_edges.add(((a / 3) + 1));
                        C_edges.add(((b / 3) + 1));
                    }
                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (!molecule1_pair_connected && !molecule2_pair_connected) {
                        D_edges.add(a / 3 + 1);
                        D_edges.add(b / 3 + 1);
                    }
                }
            }
        }

        //print C and D edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();

        return 0;
    }

//comp_graph_nodes_C_zero is used to build up of the edges of the compatibility graph
    private int generate_compatibility_graph_nodes_if_C_edge_number_is_zero() {

        int count_nodes = 1;

        for (int a = 0; a < atom_num_H_1; a++) {
            String atom1_type = atomstr1.get(a).getSymbol();

            for (int b = 0; b < atom_num_H_2; b++) {
                String atom2_type = atomstr2.get(b).getSymbol();

                if ((atom1_type.equals(atom2_type))) {
                    comp_graph_nodes_C_zero.add(a + 1);
                    comp_graph_nodes_C_zero.add(b + 1);
                    comp_graph_nodes_C_zero.add(SYMBOL_VALUE.get(atom1_type)); //C is label 1
                    comp_graph_nodes_C_zero.add(count_nodes);

                    comp_graph_nodes.add(a + 1);
                    comp_graph_nodes.add(b + 1);
                    comp_graph_nodes.add(count_nodes++);
                }
            }
        }

        return 0;
    }

    private int generate_compatibility_graph_if_C_edge_number_is_zero() {

        int vector_size = comp_graph_nodes_C_zero.size();

        for (int a = 0; a < vector_size; a = a + 4) {
            for (int b = a; b < vector_size; b = b + 4) {
                if (a != b
                        && !comp_graph_nodes_C_zero.get(a).equals(comp_graph_nodes_C_zero.get(b))
                        && !comp_graph_nodes_C_zero.get(a + 1).equals(comp_graph_nodes_C_zero.get(b + 1))) {

                    boolean molecule1_pair_connected = false;
                    boolean molecule2_pair_connected = false;

                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    for (int x = 0; x < bond_number1; x++) {
                        if ((comp_graph_nodes_C_zero.get(a).equals(i_tab1.get(x * 3 + 0))
                                && comp_graph_nodes_C_zero.get(b).equals(i_tab1.get(x * 3 + 1)))
                                || (comp_graph_nodes_C_zero.get(a).equals(i_tab1.get(x * 3 + 1))
                                && comp_graph_nodes_C_zero.get(b).equals(i_tab1.get(x * 3 + 0)))) {
                            molecule1_pair_connected = true;
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < bond_number2; y++) {

                        if ((comp_graph_nodes_C_zero.get(a + 1).equals(i_tab2.get(y * 3 + 0))
                                && comp_graph_nodes_C_zero.get(b + 1).equals(i_tab2.get(y * 3 + 1)))
                                || (comp_graph_nodes_C_zero.get(a + 1).equals(i_tab2.get(y * 3 + 1))
                                && comp_graph_nodes_C_zero.get(b + 1).equals(i_tab2.get(y * 3 + 0)))) {
                            molecule2_pair_connected = true;
                        }
                    }
                    //in case that both molecule pairs are connected a c-edge is generated
                    if ((molecule1_pair_connected && molecule2_pair_connected)) {
                        //System.out.println("C-edge " + ((a/4)+1) + " " + ((b/4)+1) << endl;
                        C_edges.add(((a / 4) + 1));
                        C_edges.add(((b / 4) + 1));
                    }
                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (!molecule1_pair_connected && !molecule2_pair_connected) {
                        D_edges.add(((a / 4) + 1));
                        D_edges.add(((b / 4) + 1));
                    }
                }
            }
        }

        //print C and D edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();

        return 0;
    }

    /*
   
     * C: set of vertices belonging to the current clique
     
     * S: set of vertices which are not allowed to be added
     * to C
    
     * P: is a set of vertices which <b>can</b> be added to C, because they are
     * neighbours of vertex u via <i>c-edges</i>
    
     * D: is a set of vertices which <b>cannot</b> be added to C, because they are
     * neighbours of vertex u via <i>d-edges</i>
     
     * V: stored all the vertices for the Graph G
     * V[G]: nodes of vector comp_graph_nodes are stored in V
     
     */
    private int Init_Algorithm() {

//        System.out.println( "Init_Algorithm " + comp_graph_nodes.size() );
        List<Integer> C = new ArrayList<>();
        List<Integer> D = new ArrayList<>();
        List<Integer> S = new ArrayList<>();
        List<Integer> N = new ArrayList<>();
        Stack<Integer> P = new Stack<>();

        //nodes of vector comp_graph_nodes are stored in V
        List<Integer> V = new ArrayList<>();//Initialization of Vector V

        int V_set_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
        }

        V.add(0);

        int b = 0;

        while (V.get(b) != 0) { // V[b] is node u
            int central_node = V.get(b);

            P.clear();
            D.clear();
            S.clear();
            N.clear();
            C.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(V.get(b));

            //  System.out.println("N-Neigh: " + N);
            for (int c = 0; c < N.size(); c = c + 2) { // N[c] is node v
                //System.out.println("N[" + c  +  "]= " + N.get(c) + " ");

                //Grouping of the neighbors in S,P and D
                if (N.get(c + 1) == 1) {    //u and v are adjacent via a C-edge
                    boolean Nc_belongs_to_T = false;
                    int T_size = T.size();
                    // System.out.println("T_size " +T_size);
                    for (int e = 0; e < T_size; e++) {
                        //    System.out.println("T[" + e  +  "]= " +T.get(e));
                        if (N.get(c).intValue() == (T.get(e))) {

                            S.add(N.get(c));
                            Nc_belongs_to_T = true;
                        }
                    }
                    if (Nc_belongs_to_T == false) {
                        P.add(N.get(c));
                    }
                } else if (N.get(c + 1) == 2) {   // u and v are adjacent via a D-edge
                    //  System.out.println("u and v are adjacent via a D-edge");
                    D.add(N.get(c));
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
                    V.remove(V.size() - 1);
                    if (neighbor_position < b) {
                        b = b - 1;
                    }
                }
            }
            P.add(0);
            C.add(central_node);
            Enumerate_Cliques(C, P, D, S);
            T.add(V.get(b));
            b++;
        }

        return 0;
    }

    private int Enumerate_Cliques(List<Integer> C, Stack<Integer> P, List<Integer> D, List<Integer> S) {

        List<Integer> N = new ArrayList<>();////Initialization Vector N
        Stack<Integer> ut_set = new Stack<>();

        P.stream().forEach((I) -> {
            ut_set.add(I);
        });

        List<Integer> C_copy = new ArrayList<>();
        Stack<Integer> P_copy = new Stack<>();
        Stack<Integer> D_copy = new Stack<>();
        List<Integer> S_copy = new ArrayList<>();

        if (P.size() == 1) {
            if (S.isEmpty()) {

                //store best solutions in stack Max_Cliques_Set
                int clique_size = C.size();
                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {
                        while (!Max_Cliques_Set.empty()) {
                            Max_Cliques_Set.pop();
                        }
                        best_clique_size = clique_size;
                    }
                    if (clique_size == best_clique_size) {
                        Max_Cliques_Set.push(C);
                    }
                }
                // System.out.println("Max_Cliques_Set: " + Max_Cliques_Set.size());
                return 0;
            }
        }
        int a = 0;
        while (!ut_set.get(a).equals(0)) { // P[a] is node ut

            int ui = ut_set.get(a);
            //remove ut_set[a] from P
            //find position of ut_set node in P
            int P_size = P.size();
            int ut_node_pos = 100000;
            for (int counter = 0; counter < P_size - 1; counter++) {  //-1 wegen Endekennung
                if (P.get(counter).equals(ut_set.get(a))) {
                    ut_node_pos = counter;
                }
            }
            if (ut_node_pos == 100000) {
                System.out.println("ut_node_pos = 100000");
            }
            //delete ut_set node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++) {
                P.setElementAt(P.get(counter + 1), counter);
            }
            P.pop();//TO DO

            C_copy.clear();
            P_copy.clear();
            D_copy.clear();
            S_copy.clear();
            N.clear();

            C.stream().forEach((obj) -> {
                C_copy.add(obj);
            });

            P.stream().forEach((obj) -> {
                P_copy.add(obj);
            });
            D.stream().forEach((obj) -> {
                D_copy.add(obj);
            });
            S.stream().forEach((obj) -> {
                S_copy.add(obj);
            });
            P_copy.pop();

            //find the neighbors of the central node from P
            // System.out.println("ut_set.get(a)" + ut_set.get(a));
            N = find_neighbors(ut_set.get(a));

            int N_size = N.size();

            for (int b = 0; b < N_size; b = b + 2) { // N[b] is node v

                int D_set_size = D.size();  //gehe Nachbarn durch, die ber D-edge mit central_node verbunden sind
                for (int c = 0; c < D_set_size; c++) {

                    if (N.get(b).equals(D.get(c))) {
                        if (N.get(b + 1) == 1) {     //u and v are adjacent via a C-edge
                            boolean Nb_belongs_to_T = false;
                            int T_size = T.size();
                            for (int d = 0; d < T_size; d++) {
                                if (N.get(b).intValue() == T.get(d)) {
                                    S_copy.add(N.get(b));
                                    Nb_belongs_to_T = true;
                                }
                            }
                            //store N[b] in P
                            if (Nb_belongs_to_T == false) {
                                P_copy.add(N.get(b));
                            }
                            //delete N[b] bzw. D[c] from set D_copy
                            int D_copy_size = D_copy.size();
                            int Nb_position = 10000;
                            for (int e = 0; e < D_copy_size - 1; e++) {
                                if (N.get(b).equals(D_copy.get(e))) {
                                    Nb_position = e;
                                }
                            }
                            for (int e = Nb_position; e < D_copy_size - 1; e++) {
                                D_copy.setElementAt(D_copy.get(e + 1), e);
                            }
                            D_copy.pop();
                        }
                        /*//Abschnitt sinnlos, denn wenn etwas in S war ist, es nach S' kopiert worden
                     if(N[b+1] == 2){     //u and v are adjacent via a D-edge
                     if().....
                     }*/
                    }
                }
                //find respective neighbor position in ut_set, which is needed for the deletion from ut_set
                int ut_set_size = ut_set.size();
                int neighbor_position = -1;
                for (int e = 0; e < ut_set_size; e++) {
                    if (N.get(b).equals(ut_set.get(e))) {
                        neighbor_position = e;
                    }
                }
                if (neighbor_position != -1) {
                    //delete neighbor from set P
                    for (int e = neighbor_position; e < ut_set_size - 1; e++) {
                        ut_set.setElementAt(ut_set.get(e + 1), e);
                    }
                    ut_set.pop(); //TODO:Check whether size returns number of elements or index value
                    if (neighbor_position < a) {
                        a = a - 1;
                    }
                }
            }

            Stack<Integer> P_copy_N_intersec = new Stack<>();
            List<Integer> D_copy_N_intersec = new ArrayList<>();
            List<Integer> S_copy_N_intersec = new ArrayList<>();

            int P_copy_length = P_copy.size();
            int D_copy_length = D_copy.size();
            int S_copy_length = S_copy.size();

            for (int sec = 0; sec < N_size; sec = sec + 2) {
                for (int pc = 0; pc < P_copy_length; pc++) {
                    if (P_copy.get(pc).equals(N.get(sec))) {
                        P_copy_N_intersec.add(P_copy.get(pc));
                    }
                }
                for (int pd = 0; pd < D_copy_length; pd++) {
                    if (D_copy.get(pd).equals(N.get(sec))) {
                        D_copy_N_intersec.add(D_copy.get(pd));
                    }
                }
                for (int ps = 0; ps < S_copy_length; ps++) {
                    if (S_copy.get(ps).equals(N.get(sec))) {
                        S_copy_N_intersec.add(S_copy.get(ps));
                    }
                }
            }
            P_copy_N_intersec.push(0);
            C_copy.add(ui);
            Enumerate_Cliques(C_copy, P_copy_N_intersec, D_copy_N_intersec, S_copy_N_intersec);
            S.add(ui);
            a++;
        }

        return 0;
    }

    private List<Integer> find_neighbors(int central_node) {

        List<Integer> neighbor_vec = new ArrayList<>();

        // System.out.println("C_edge.size: " + C_edges.size());
        int C_edge_number = C_edges.size() / 2;

        //    System.out.println("C_edge.size/2: " + C_edge_number);
        //    System.out.println("");
        //    System.out.println("C_edges: ");
        for (int a = 0; a < C_edge_number; a++) {
            if (C_edges.get(a * 2 + 0).equals(central_node)) {
                //          System.out.println( C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 1));
                neighbor_vec.add(1);       // 1 means: is connected via C-edge
            }
            if (C_edges.get(a * 2 + 1).equals(central_node)) {
                //           System.out.println(C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 0));
                neighbor_vec.add(1);       // 1 means: is connected via C-edge
            }
        }

        int D_edge_number = D_edges.size() / 2;

        //System.out.println("");
        //System.out.println("D_edges: "+ D_edges.size());
        for (int a = 0; a < D_edge_number; a++) {
            if (D_edges.get(a * 2 + 0).equals(central_node)) {
                //       System.out.println( D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 1));
                neighbor_vec.add(2);       // 2 means: is connected via D-edge
            }
            if (D_edges.get(a * 2 + 1).equals(central_node)) {
                //        System.out.println(D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 0));
                neighbor_vec.add(2);       // 2 means: is connected via D-edge
            }
        }

        return neighbor_vec;
    }

    //extract atom mapping from the clique vector and print it on the screen
    int extract_mapping(List<Integer> clique_vector) {

        List<Integer> temp_vector = new ArrayList<>();
        temp_vector.clear();

        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b = b + 3) {
                if (clique_vector.get(a).equals(comp_graph_nodes.get(b + 2))) {
                    temp_vector.add(comp_graph_nodes.get(b));
                    temp_vector.add(comp_graph_nodes.get(b + 1));
                }
            }
        }

        getFinalMappings().add(temp_vector);

        return 0;
    }

//extract atom mapping from the clique vector and store it in vector clique_MAPPING_Local
    private List<Integer> extract_clique_MAPPING(List<Integer> clique_vector) {

        List<Integer> clique_MAPPING_Local = new ArrayList<>();

        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b = b + 3) {
                if (clique_vector.get(a).equals(comp_graph_nodes.get(b + 2))) {
                    clique_MAPPING_Local.add(comp_graph_nodes.get(b));
                    clique_MAPPING_Local.add(comp_graph_nodes.get(b + 1));
                }
            }
        }

        return clique_MAPPING_Local;
    }

//Function is called by the main program and serves as a starting point for the comparision procedure.
    public int search_cliques() {

        generate_compatibility_graph_nodes();
        generate_compatibility_graph();
        if (C_edges_size == 0) {
            comp_graph_nodes.clear();
            C_edges.clear();
            D_edges.clear();
            C_edges_size = 0;
            D_edges_size = 0;
            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            generate_compatibility_graph_if_C_edge_number_is_zero();
            comp_graph_nodes_C_zero.clear();
        }

        best_clique_size = 0;
        T.clear();
        Init_Algorithm();
        best_MAPPING_size = 0;

//        int clique_number = 1;
        while (!Max_Cliques_Set.empty()) {
//            System.out.println ("Clique number " + clique_number + " :" );
            List<Integer> clique_vector = Max_Cliques_Set.peek();
            int clique_size = clique_vector.size();
            //Is the number of mappings smaller than the number of atoms of molecule A and B?
            //In this case the clique is given to the McGregor algorithm
            if ((clique_size < atom_number1) && (clique_size < atom_number2)) {
//                System.out.print( "clique_size: "+clique_size+" atom_number1: "+atom_number1+" atom_number2: "+atom_number2);
//                System.out.println(" -> McGregor");
                McGregor_IterationStart(clique_vector);

            } else {
                List<Integer> clique_MAPPING = extract_clique_MAPPING(clique_vector);
                extract_mapping(clique_vector);
            }
            Max_Cliques_Set.pop();
//            clique_number++;
        }

        postfilter();
        this.Max_Cliques_Set.clear();
        this.comp_graph_nodes.clear();
        this.comp_graph_nodes_C_zero.clear();
        this.c_tab1.clear();
        this.c_tab2.clear();
        this.C_edges.clear();
        this.D_edges.clear();
        this.T.clear();
        this.C_edges_size = 0;
        this.D_edges_size = 0;
        return 0;
    }

}
