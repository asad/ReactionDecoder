/* Copyright (R) 2009-2018  Syed Asad Rahman <asad@ebi.ac.uk>
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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, R. Karande: An Algorithm for reporting maximal
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

    private final boolean DEBUG = false;

    /**
     * Creates a new instance of SearchCliques
     *
     *
     * @param f1
     * @param f2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public MCSPlus(IAtomContainer f1, IAtomContainer f2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {

        super(f1, f2, shouldMatchBonds, shouldMatchRings, matchAtomType);

    }

    private List<List<Integer>> label_atoms(List<Integer> basic_atom_vector, int bond_num, List<IAtom> atoms, List<Integer> i_tab, List<String> c_tab) {

        ArrayList<List<Integer>> label_list = new ArrayList<>();

//        if (DEBUG) {
//            System.out.println("Vector Atom Str: ");
//            for (int b = 0; b < atoms.size(); b++) {
//                System.err.print(atoms.get(b).getSymbol() + ",");
//            }
//            System.err.println();
//            System.err.println("basic_atom_vector");
//            for (int b = 0; b < basic_atom_vector.size(); b++) {
//                System.err.print(basic_atom_vector.get(b) + ",");
//            }
//            System.err.println();
//            System.err.println("i_tab");
//            for (int b = 0; b < i_tab.size(); b++) {
//                System.err.print(i_tab.get(b) + ",");
//            }
//            System.err.println();
//            System.err.println("c_tab");
//            for (int b = 0; b < c_tab.size(); b++) {
//                System.err.print(c_tab.get(b) + ",");
//            }
//            System.err.println();
//        }
        for (int a = 0; a < basic_atom_vector.size(); a++) {

            List<Integer> label = new ArrayList<>(7);
            /*
             * Initialize the vector
             */
            for (int i = 0; i < 7; i++) {
                label.add(0);
            }

            IAtom atom1 = atoms.get(a);
            String atom1_type = atom1.getSymbol();
            if (matchAtomType && atom1.getAtomTypeName() != null) {
                atom1_type = atoms.get(a).getAtomTypeName();
            }
            if (SYMBOL_VALUE.containsKey(atom1_type)) {
                label.set(0, SYMBOL_VALUE.get(atom1_type));
            } else {
                int value = atom1.getAtomicNumber() == null ? atom1.hashCode() + 1000 : atom1.getAtomicNumber() + 1000;
                SYMBOL_VALUE.put(atom1_type, value);
                label.set(0, SYMBOL_VALUE.get(atom1_type));
            }
            int count_neighbors = 1;
            for (int b = 0; b < bond_num; b++) {
                if (basic_atom_vector.get(a).equals(i_tab.get(b * 3 + 0))) {
                    String atom2_type = c_tab.get(b * 2 + 1);
                    /*Get neighbour Atom*/
                    IAtom atom2 = atoms.get(i_tab.get(b * 3 + 1) - 1);
                    //System.out.println("atom2_type " + atom2_type + ", atom2 " + atom2.getSymbol());

                    if (matchAtomType && atom2.getAtomTypeName() != null) {
                        atom2_type = atom2.getAtomTypeName();
                    }

                    if (SYMBOL_VALUE.containsKey(atom2_type)) {
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    } else {
                        int value = atom2.getAtomicNumber() == null ? atom2.hashCode() + 1000 : atom2.getAtomicNumber() + 1000;
                        SYMBOL_VALUE.put(atom2_type, value);
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    }
                    count_neighbors++;
                }

                if (basic_atom_vector.get(a).equals(i_tab.get(b * 3 + 1))) {
                    String atom2_type = c_tab.get(b * 2 + 0);
                    /*Get neighbour Atom*/
                    IAtom atom2 = atoms.get(i_tab.get(b * 3 + 0) - 1);

                    if (matchAtomType && atom2.getAtomTypeName() != null) {
                        atom2_type = atom2.getAtomTypeName();
                    }

                    if (SYMBOL_VALUE.containsKey(atom2_type)) {
                        label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                    } else {
                        int value = atom2.getAtomicNumber() == null ? atom2.hashCode() + 1000 : atom2.getAtomicNumber() + 1000;
                        SYMBOL_VALUE.put(atom2_type, value);
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

        if (DEBUG) {
            System.out.println("label_list of Atoms: " + label_list.size());
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

        List<Integer> basic_atom_vec_A = reduce_atomset(atom_num_H_1, bond_number1, atomstr1, i_tab1, c_tab1);
        List<Integer> basic_atom_vec_B = reduce_atomset(atom_num_H_2, bond_number2, atomstr2, i_tab2, c_tab2);

        List<List<Integer>> label_list_molA = label_atoms(basic_atom_vec_A, bond_number1, atomstr1, i_tab1, c_tab1);
        List<List<Integer>> label_list_molB = label_atoms(basic_atom_vec_B, bond_number2, atomstr2, i_tab2, c_tab2);

        int molA_nodes = 0;
        int count_nodes = 1;

        for (List<Integer> labelA : label_list_molA) {
            int molB_nodes = 0;
            for (List<Integer> labelB : label_list_molB) {
                if (labelA.equals(labelB)) {
//                    System.out.println("labelA " + labelA + ", labelB " + labelB + "\n");
                    comp_graph_nodes.add(basic_atom_vec_A.get(molA_nodes));
                    comp_graph_nodes.add(basic_atom_vec_B.get(molB_nodes));
                    comp_graph_nodes.add(count_nodes++);
                }
                molB_nodes++;
            }
            molA_nodes++;
        }

        if (DEBUG) {
            System.out.println("comp_graph_nodes: " + comp_graph_nodes.size());
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

                    IBond bond1 = null;
                    IBond bond2 = null;

                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    for (int x = 0; x < bond_number1; x++) {
                        if ((comp_graph_nodes.get(a).equals(i_tab1.get(x * 3 + 0))
                                && comp_graph_nodes.get(b).equals(i_tab1.get(x * 3 + 1)))) {

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(a) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(b) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                                System.out.println("BOND " + i_tab1.get(x * 3 + 2));
//                            }
                            IAtom a1 = this.ac1.getAtom(comp_graph_nodes.get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(comp_graph_nodes.get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);
                            molecule1_pair_connected = true;
                            break;
                        } else if ((comp_graph_nodes.get(a).equals(i_tab1.get(x * 3 + 1))
                                && comp_graph_nodes.get(b).equals(i_tab1.get(x * 3 + 0)))) {

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(a) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(b) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                                System.out.println("BOND " + i_tab1.get(x * 3 + 2));
//                            }
                            IAtom a1 = this.ac1.getAtom(comp_graph_nodes.get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(comp_graph_nodes.get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);
                            molecule1_pair_connected = true;
                            break;
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < bond_number2; y++) {
                        if ((comp_graph_nodes.get(a + 1).equals(i_tab2.get(y * 3 + 0))
                                && comp_graph_nodes.get(b + 1).equals(i_tab2.get(y * 3 + 1)))) {
//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(a + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(b + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                                System.out.println("BOND " + i_tab2.get(y * 3 + 2));
//                            }
                            IAtom a1 = this.ac2.getAtom(comp_graph_nodes.get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(comp_graph_nodes.get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);
                            molecule2_pair_connected = true;
                            break;

                        } else if ((comp_graph_nodes.get(a + 1).equals(i_tab2.get(y * 3 + 1))
                                && comp_graph_nodes.get(b + 1).equals(i_tab2.get(y * 3 + 0)))) {
//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(a + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(b + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                                System.out.println("BOND " + i_tab2.get(y * 3 + 2));
//                            }
                            IAtom a1 = this.ac2.getAtom(comp_graph_nodes.get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(comp_graph_nodes.get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);
                            molecule2_pair_connected = true;
                            break;

                        }
                    }

                    boolean connectedFlag = false;
                    boolean disConnectedFlag = false;
                    boolean matchBondFlag = false;

                    if (molecule1_pair_connected
                            && molecule2_pair_connected) {
                        connectedFlag = true;
                    }

                    if (!molecule1_pair_connected
                            && !molecule2_pair_connected) {
                        disConnectedFlag = true;
                    }

                    if (connectedFlag
                            && isMatchFeasible(bond1, bond2, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                        matchBondFlag = true;
                    }

//                    if (DEBUG) {
//                        System.out.println("matchbondFlag " + connectedFlag);
//                    }
                    //in case that both molecule pairs are connected a c-edge is generated
                    if (connectedFlag && matchBondFlag) {

                        C_edges.add(((a / 3) + 1));
                        C_edges.add(((b / 3) + 1));
                    }
//
                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (disConnectedFlag) {
                        D_edges.add((a / 3) + 1);
                        D_edges.add((b / 3) + 1);
                    }

                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (connectedFlag && !matchBondFlag) {
                        D_edges.add((a / 3) + 1);
                        D_edges.add((b / 3) + 1);
                    }
                }
            }
        }

        //print R and Q edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();

        if (DEBUG) {
            System.out.println("C_edges_size " + C_edges_size);
            System.out.println("D_edges_size " + D_edges_size);
        }

        return 0;
    }

//comp_graph_nodes_C_zero is used to build up of the edges of the compatibility graph
    private int generate_compatibility_graph_nodes_if_C_edge_number_is_zero() {

        int count_nodes = 1;

        for (int a = 0; a < atom_num_H_1; a++) {
            String atom1_type = atomstr1.get(a).getSymbol();
            int value = atomstr1.get(a).getAtomicNumber() == null ? atomstr1.get(a).hashCode() + 1000 : atomstr1.get(a).getAtomicNumber() + 1000;
            SYMBOL_VALUE.put(atom1_type, value);
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

                    IBond bond1 = null;
                    IBond bond2 = null;
                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    for (int x = 0; x < bond_number1; x++) {
                        if ((comp_graph_nodes_C_zero.get(a).equals(i_tab1.get(x * 3 + 0))
                                && comp_graph_nodes_C_zero.get(b).equals(i_tab1.get(x * 3 + 1)))) {
                            molecule1_pair_connected = true;
                            IAtom a1 = this.ac1.getAtom(comp_graph_nodes_C_zero.get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(comp_graph_nodes_C_zero.get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes_C_zero.get(a) " + comp_graph_nodes_C_zero.get(a) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                                System.out.println("comp_graph_nodes_C_zero.get(a) " + comp_graph_nodes_C_zero.get(b) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                            }
                            break;

                        } else if ((comp_graph_nodes_C_zero.get(a).equals(i_tab1.get(x * 3 + 1))
                                && comp_graph_nodes_C_zero.get(b).equals(i_tab1.get(x * 3 + 0)))) {
                            molecule1_pair_connected = true;
                            IAtom a1 = this.ac1.getAtom(comp_graph_nodes_C_zero.get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(comp_graph_nodes_C_zero.get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes_C_zero.get(a) " + comp_graph_nodes_C_zero.get(a) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                                System.out.println("comp_graph_nodes_C_zero.get(a) " + comp_graph_nodes_C_zero.get(b) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                            }
                            break;
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < bond_number2; y++) {

                        if ((comp_graph_nodes_C_zero.get(a + 1).equals(i_tab2.get(y * 3 + 0))
                                && comp_graph_nodes_C_zero.get(b + 1).equals(i_tab2.get(y * 3 + 1)))) {
                            molecule2_pair_connected = true;
                            IAtom a1 = this.ac2.getAtom(comp_graph_nodes_C_zero.get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(comp_graph_nodes_C_zero.get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes_C_zero.get(a+1) " + comp_graph_nodes_C_zero.get(a + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                                System.out.println("comp_graph_nodes_C_zero.get(a+1) " + comp_graph_nodes_C_zero.get(b + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                            }
                            break;
                        } else if ((comp_graph_nodes_C_zero.get(a + 1).equals(i_tab2.get(y * 3 + 1))
                                && comp_graph_nodes_C_zero.get(b + 1).equals(i_tab2.get(y * 3 + 0)))) {
                            molecule2_pair_connected = true;
                            IAtom a1 = this.ac2.getAtom(comp_graph_nodes_C_zero.get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(comp_graph_nodes_C_zero.get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes_C_zero.get(a+1) " + comp_graph_nodes_C_zero.get(a + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                                System.out.println("comp_graph_nodes_C_zero.get(a+1) " + comp_graph_nodes_C_zero.get(b + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                            }
                            break;
                        }
                    }

                    boolean connectedFlag = false;
                    boolean disConnectedFlag = false;
                    boolean matchBondFlag = false;

                    if (molecule1_pair_connected
                            && molecule2_pair_connected) {
                        connectedFlag = true;
                    }

                    if (!molecule1_pair_connected
                            && !molecule2_pair_connected) {
                        disConnectedFlag = true;
                    }

                    if (connectedFlag
                            && isMatchFeasible(bond1, bond2, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                        matchBondFlag = true;
                    }

//                    if (DEBUG) {
//                        System.out.println("matchbondFlag " + connectedFlag);
//                    }
                    //in case that both molecule pairs are connected a c-edge is generated
                    if (connectedFlag && matchBondFlag) {

                        C_edges.add(((a / 4) + 1));
                        C_edges.add(((b / 4) + 1));
                    }
//
                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (disConnectedFlag) {
                        D_edges.add((a / 4) + 1);
                        D_edges.add((b / 4) + 1);
                    }

                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (connectedFlag && !matchBondFlag) {
                        D_edges.add((a / 4) + 1);
                        D_edges.add((b / 4) + 1);
                    }
                }
            }
        }

        //print R and Q edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();
        if (DEBUG) {
            System.out.println("C_edges_size " + C_edges_size);
            System.out.println("D_edges_size " + D_edges_size);
        }
        return 0;
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
//        System.out.println("C_edges_size " + C_edges_size);
//        System.out.println("bond cound " + ac1.getBondCount());
//        System.out.println("bond cound " + ac2.getBondCount());
        if (C_edges_size == 0
                || ((C_edges_size < this.ac1.getAtomCount() / 2
                && C_edges_size < this.ac2.getAtomCount() / 2))
                && (this.ac1.getAtomCount() / 2 < 30
                && this.ac2.getAtomCount() / 2 < 30)) {

//            System.out.println("Switching to complex mode ");
            comp_graph_nodes.clear();
            C_edges.clear();
            D_edges.clear();
            C_edges_size = 0;
            D_edges_size = 0;
            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            generate_compatibility_graph_if_C_edge_number_is_zero();
            comp_graph_nodes_C_zero.clear();
        }

        BKKCKCF cliqueFinder = new BKKCKCF(comp_graph_nodes, C_edges, D_edges);
        cliqueFinder.init_Algorithm();
        this.max_Cliques_Set = cliqueFinder.getMax_Cliques_Set();

        best_MAPPING_size = 0;

//        int clique_number = 1;
        while (!max_Cliques_Set.empty()) {
//            System.out.println ("Clique number " + clique_number + " :" );
            List<Integer> clique_vector = max_Cliques_Set.peek();
            int clique_size = clique_vector.size();
            //Is the number of mappings smaller than the number of atoms of molecule A and B?
            //In this case the clique is given to the McGregor algorithm
            if ((clique_size < atom_number1) && (clique_size < atom_number2)) {
//                System.out.print("clique_size: " + clique_vector + " atom_number1: " + atom_number1 + " atom_number2: " + atom_number2);
//                System.out.println(" -> McGregor");
                try {
                    McGregor_IterationStart(clique_vector);
                } catch (Exception e) {
                    e.printStackTrace();
                }

            } else {
                List<Integer> clique_MAPPING = extract_clique_MAPPING(clique_vector);
                extract_mapping(clique_vector);
            }
            max_Cliques_Set.pop();
//            clique_number++;
        }

        postfilter();

        return 0;
    }

    private void clear() {
        this.max_Cliques_Set.clear();
        this.comp_graph_nodes.clear();
        this.comp_graph_nodes_C_zero.clear();
        this.c_tab1.clear();
        this.c_tab2.clear();
        this.C_edges.clear();
        this.D_edges.clear();
        this.C_edges_size = 0;
        this.D_edges_size = 0;
    }

}
