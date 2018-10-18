
/* Copyright (C) 2005-2006 Markus Leber
 *               2006-2018 Syed Asad Rahman <asad at ebi.ac.uk>
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
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.smsd.algorithm.mcsplus.mcsplus1.BinaryTree.remove_tree_structure;

/**
 * Class which reports MCS solutions based on the McGregor algorithm published
 * in 1982.
 *
 * <p>
 * The SMSD algorithm is described in this paper. <font color="#FF0000">please
 * refer Rahman <i>et.al. 2009</i></font> {
 *
 * @cdk.cite SMSD2009}. </p>
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class McGregor extends Utility {

    protected final List<Integer> c_edges;
    protected final List<Integer> d_edges;

    private List<String> c_tab1_copy;
    private List<String> c_tab2_copy;

    private BinaryTree last, first;

    protected int c_edges_size;
    protected int d_edges_size;

    protected int best_clique_size;

    protected final List<Integer> comp_graph_nodes;
    protected final List<Integer> comp_graph_nodes_C_zero;

    protected final List<IAtom> atomstr1;
    protected final List<IAtom> atomstr2;

    private int nNum_globalA;
    private int nNum_globalB;

    private final List<Integer> i_globalA;
    private final List<Integer> i_globalB;
    private final List<String> c_globalA;
    private final List<String> c_globalB;

    private List<Integer> MARCS;
    private List<Integer> FIXARCS;
    private final Stack<List<Integer>> BESTARCS;
    private int bestarcsleft;

    protected int atom_number1;
    protected int atom_number2;

    protected int atom_num_H_1;
    protected int atom_num_H_2;

    protected int bond_number1;
    protected int bond_number2;

    protected List<Integer> i_tab1;
    protected List<Integer> i_tab2;

    protected List<String> c_tab1;
    protected List<String> c_tab2;

    private boolean new_matrix;
    protected int best_MAPPING_size;
    protected Stack<List<Integer>> max_Cliques_Set;

    private final List<List<Integer>> final_MAPPINGS;

    private final List<String> SignROW;
    protected final IAtomContainer ac1;
    protected final IAtomContainer ac2;
    protected final boolean shouldMatchBonds;
    protected final boolean shouldMatchRings;
    protected final boolean matchAtomType;
    protected final Map<String, Integer> SYMBOL_VALUE;

    /**
     *
     * @param f1
     * @param f2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public McGregor(IAtomContainer f1, IAtomContainer f2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomType = matchAtomType;

        this.SYMBOL_VALUE = new TreeMap<>();

        MoleculeHandler file1 = new MoleculeHandler(f1, shouldMatchBonds);
        MoleculeHandler file2 = new MoleculeHandler(f2, shouldMatchBonds);

        this.atom_number1 = file1.indexOf();
        this.atom_number2 = file2.indexOf();
        this.atom_num_H_1 = file1.getStartHatom_num();
        this.atom_num_H_2 = file2.getStartHatom_num();
        this.bond_number1 = file1.getBondNumber();
        this.bond_number2 = file2.getBondNumber();
        this.atomstr1 = file1.getAtomString();
        this.atomstr2 = file2.getAtomString();
        this.i_tab1 = file1.intTable;
        this.i_tab2 = file2.intTable;
        this.c_tab1 = file1.charTable;
        this.c_tab2 = file2.charTable;
        this.ac1 = file1.getAtomContainer();
        this.ac2 = file2.getAtomContainer();

        this.comp_graph_nodes = new ArrayList<>();
        this.c_edges = new ArrayList<>();//Initialize the c_edges Vector
        this.d_edges = new ArrayList<>();//Initialize the d_edges Vector
        this.comp_graph_nodes_C_zero = new ArrayList<>();//Initialize the comp_graph_nodes_C_zero Vector

        this.c_tab1_copy = new ArrayList<>();
        this.c_tab2_copy = new ArrayList<>();

        this.nNum_globalA = 0;
        this.nNum_globalB = 0;
        this.i_globalA = new ArrayList<>();
        this.i_globalB = new ArrayList<>();
        this.c_globalA = new ArrayList<>();
        this.c_globalB = new ArrayList<>();

        this.MARCS = new ArrayList<>();
        this.FIXARCS = new ArrayList<>();
        this.BESTARCS = new Stack<>();
        this.final_MAPPINGS = new ArrayList<>(); //Initialization of Vector final_MAPPINGS

        this.max_Cliques_Set = new Stack<>(); //Initialization max_Cliques_Set

        //String SignROW = "DGJTVWYZ$%&*#?!~^<>=()[]";
        String[] characters = {"D", "G", "J", "T", "V", "W", "Y", "Z", "$", "%", "&", "*", "#", "?", "!", "~", "^", "<", ">", "=", "(", ")", "[", "]"};
        this.SignROW = Arrays.asList(characters);

    }

    /**
     * ******************************************************************************************************************
     */
    /* MacGregor Methods */
    //TO DO Asad
    /**
     * @param MARCS_vector
     *
     * @param mapped_atoms_num
     * @param current_MAPPING
     * @param bondnum_A
     * @param i_bonds_A
     * @param bondnum_B
     * @param i_bonds_B
     * @return
     * *****************************************************************************************************************
     */
    /*get atom mappings from the McGregor solution matrices
     */
    protected List<Integer> find_mcgregor_MAPPING(List<Integer> MARCS_vector, int mapped_atoms_num, List<Integer> current_MAPPING, int bondnum_A, List<Integer> i_bonds_A, int bondnum_B, List<Integer> i_bonds_B) {
        List<Integer> additional_mapping = new ArrayList<>();
        additional_mapping.clear();
        int pos = 0;
        int number_of_ones = 0;
        for (int x = 0; x < bondnum_A; x++) {
            for (int z = 0; z < bondnum_B; z++) {
                if (MARCS_vector.get(x * bondnum_B + z) == 1) {
                    int cur_pos = x * nNum_globalB + z;
                    int Atom1_moleculeA = i_bonds_A.get(x * 3 + 0);
                    int Atom2_moleculeA = i_bonds_A.get(x * 3 + 1);
                    int Atom1_moleculeB = i_bonds_B.get(z * 3 + 0);
                    int Atom2_moleculeB = i_bonds_B.get(z * 3 + 1);
                    for (int a = 0; a < mapped_atoms_num; a++) {
                        if ((current_MAPPING.get(a * 2 + 0) == (Atom1_moleculeA)) && (current_MAPPING.get(a * 2 + 1) == (Atom1_moleculeB))) {
                            additional_mapping.add(Atom2_moleculeA);
                            additional_mapping.add(Atom2_moleculeB);
                        }
                        if ((current_MAPPING.get(a * 2 + 0) == (Atom1_moleculeA)) && (current_MAPPING.get(a * 2 + 1) == (Atom2_moleculeB))) {
                            additional_mapping.add(Atom2_moleculeA);
                            additional_mapping.add(Atom1_moleculeB);
                        }
                        if ((current_MAPPING.get(a * 2 + 0) == (Atom2_moleculeA)) && (current_MAPPING.get(a * 2 + 1) == (Atom1_moleculeB))) {
                            additional_mapping.add(Atom1_moleculeA);
                            additional_mapping.add(Atom2_moleculeB);
                        }
                        if ((current_MAPPING.get(a * 2 + 0) == (Atom2_moleculeA)) && (current_MAPPING.get(a * 2 + 1) == (Atom2_moleculeB))) {
                            additional_mapping.add(Atom1_moleculeA);
                            additional_mapping.add(Atom1_moleculeB);
                        }
                    }
                }
            }
        }
        int additional_mapping_size = additional_mapping.size();
        //add McGregor mapping to the Clique mapping
        for (int a = 0; a < additional_mapping_size; a = a + 2) {
            current_MAPPING.add(additional_mapping.get(a));
            current_MAPPING.add(additional_mapping.get(a + 1));
        }
        //remove recurring mappings from current_MAPPING
        List<Integer> unique_MAPPING = remove_recurring_mappings(current_MAPPING);
        return unique_MAPPING;
    }

    protected int McGregor_IterationStart(List<Integer> clique_vector) {
        c_tab1_copy.clear();
        generate_c_tab1_copy();
        c_tab2_copy.clear();
        generate_c_tab2_copy();
        List<Integer> mapped_atoms = new ArrayList<>();
        int mapped_atoms_number = 0;
        int neighbor_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS
        int set_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
        int neighbor_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS
        int set_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
        List<Integer> i_bond_neighborsA = new ArrayList<>();
        List<Integer> i_bond_setA = new ArrayList<>();
        List<String> c_bond_neighborsA = new ArrayList<>();
        List<String> c_bond_setA = new ArrayList<>();
        List<Integer> i_bond_neighborsB = new ArrayList<>();
        List<Integer> i_bond_setB = new ArrayList<>();
        List<String> c_bond_neighborsB = new ArrayList<>();
        List<String> c_bond_setB = new ArrayList<>();
        //clear vectors
        i_bond_neighborsA.clear();
        i_bond_setA.clear();
        c_bond_neighborsA.clear();
        c_bond_setA.clear();
        i_bond_neighborsB.clear();
        i_bond_setB.clear();
        c_bond_neighborsB.clear();
        c_bond_setB.clear();

//find mapped atoms of both molecules and store these in mapped_atoms
        mapped_atoms.clear();
        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b = b + 3) {
                if (clique_vector.get(a).intValue() == comp_graph_nodes.get(b + 2)) {
                    mapped_atoms.add(comp_graph_nodes.get(b));
                    mapped_atoms.add(comp_graph_nodes.get(b + 1));
                    mapped_atoms_number++;
                }
            }
        }
//        System.out.print("MoleculeA: " + "\n" + "Already mapped: ");
//        for (int a = 0; a < clique_siz; a++) {
//            System.out.print(mapped_atoms.get(a * 2 + 0) + " ");
//        }
//        System.out.println("");
//
//        System.out.print("MoleculeB: " + "\n" + "Already mapped: ");
//        for (int a = 0; a < clique_siz; a++) {
//            System.out.print(mapped_atoms.get(a * 2 + 1) + " ");
//        }
//        System.out.println();
        //find unmapped atoms of molecule A
        List<Integer> unmapped_atoms_molA = new ArrayList<>();
        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;
        for (int a = 1; a <= atom_num_H_1; a++) {
            //Atomliste sind nur Zahlen von 1 bis atom_number1
            for (int b = 0; b < clique_siz; b++) {
                //da Knotenanzahl == Anzahl zugeordneter Paare
                //cout << mapped_atoms[b*2] <<" ";
                if (a == mapped_atoms.get(b * 2)) {
                    atomA_is_unmapped = false;
                }
            }
            if (atomA_is_unmapped == true) {
                unmapped_atoms_molA.add(a);
                unmapped_numA++;
            }
            atomA_is_unmapped = true;
        }
        //Extract bonds which are related with unmapped atoms of molecule A.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molA, which contain those
        //bonds of molecule A, which are relevant for the McGregor algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule B
        int SR_count = 0;
        boolean bond_considered = false;
        boolean normal_bond = true;
        for (int a = 0; a < bond_number1; a++) {
            for (int b = 0; b < unmapped_numA; b++) {
                if (unmapped_atoms_molA.get(b).intValue() == i_tab1.get(a * 3 + 0)) {
                    for (int c = 0; c < clique_siz; c++) {
                        if (mapped_atoms.get(c * 2).intValue() == i_tab1.get(a * 3 + 1)) {
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 0));
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 1));
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 2));
                            if (c_tab1_copy.get(a * 4 + 3).equals("X")) {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(SignROW.get(SR_count));
                                c_bond_neighborsA.add("X");
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_tab1_copy = change_char_bonds(i_tab1.get(a * 3 + 1), SignROW.get(SR_count), bond_number1, i_tab1, c_tab1_copy);
                                int cor_atom = search_corresponding_atom(clique_siz, i_tab1.get(a * 3 + 1), 1, mapped_atoms);
                                c_tab2_copy = change_char_bonds(cor_atom, SignROW.get(SR_count), bond_number2, i_tab2, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_A++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setA.add(i_tab1.get(a * 3 + 0));
                        i_bond_setA.add(i_tab1.get(a * 3 + 1));
                        i_bond_setA.add(i_tab1.get(a * 3 + 2));
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 0));
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 1));
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
                        set_bondnum_A++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmapped_atoms_molA.get(b).intValue() == i_tab1.get(a * 3 + 1)) {
                    for (int c = 0; c < clique_siz; c++) {
                        if (mapped_atoms.get(c * 2 + 0).equals(i_tab1.get(a * 3 + 0))) {
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 0));
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 1));
                            i_bond_neighborsA.add(i_tab1.get(a * 3 + 2));
                            if (c_tab1_copy.get(a * 4 + 2).equals("X")) {
                                c_bond_neighborsA.add(SignROW.get(SR_count));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add("X");
                                c_tab1_copy = change_char_bonds(i_tab1.get(a * 3 + 0), SignROW.get(SR_count), bond_number1, i_tab1, c_tab1_copy);
                                int cor_atom = search_corresponding_atom(clique_siz, i_tab1.get(a * 3 + 0), 1, mapped_atoms);
                                c_tab2_copy = change_char_bonds(cor_atom, SignROW.get(SR_count), bond_number2, i_tab2, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_A++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setA.add(i_tab1.get(a * 3 + 0));
                        i_bond_setA.add(i_tab1.get(a * 3 + 1));
                        i_bond_setA.add(i_tab1.get(a * 3 + 2));
                        c_bond_setA.add(c_tab1.get(a * 2 + 0));
                        c_bond_setA.add(c_tab1.get(a * 2 + 1));
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
                        set_bondnum_A++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (bond_considered) {
                    break;
                }
            }
            bond_considered = false;
        }
//
//        System.out.println("Bonds of the neighbor set A: ");
//        for (int a = 0; a < neighbor_bondnum_A; a++) {
//            System.out.println(i_bond_neighborsA.get(a * 3 + 0) + " " + i_bond_neighborsA.get(a * 3 + 1) + " " + i_bond_neighborsA.get(a * 3 + 2));
//            System.out.print(c_bond_neighborsA.get(a * 4 + 0) + " " + c_bond_neighborsA.get(a * 4 + 1) + " : " + c_bond_neighborsA.get(a * 4 + 2) + " ");
//            System.out.println(c_bond_neighborsA.get(a * 4 + 3));
//        }
//        System.out.println("");
//        System.out.println("remaining Bond set A: ");
//        for (int a = 0; a < set_bondnum_A; a++) {
//            System.out.println(i_bond_setA.get(a * 3 + 0) + " " + i_bond_setA.get(a * 3 + 1) + " " + i_bond_setA.get(a * 3 + 2));
//            System.out.print(c_bond_setA.get(a * 4 + 0) + " " + c_bond_setA.get(a * 4 + 1) + " : " + c_bond_setA.get(a * 4 + 2) + " ");
//            System.out.println(c_bond_setA.get(a * 4 + 3));
//        }
//        System.out.println("");
        //find unmapped atoms of molecule B
        List<Integer> unmapped_atoms_molB = new ArrayList<>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;
        for (int a = 1; a <= atom_num_H_2; a++) {
            for (int b = 0; b < clique_siz; b++) {
                if (a == mapped_atoms.get(b * 2 + 1)) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped == true) {
                unmapped_atoms_molB.add(a);
                unmapped_numB++;
            }
            atomB_is_unmapped = true;
        }
        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregor algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A
        bond_considered = false;
        normal_bond = true;
        for (int a = 0; a < bond_number2; a++) {
            for (int b = 0; b < unmapped_numB; b++) {
                if (unmapped_atoms_molB.get(b).intValue() == i_tab2.get(a * 3 + 0)) {
                    for (int c = 0; c < clique_siz; c++) {
                        if (mapped_atoms.get(c * 2 + 1).intValue() == i_tab2.get(a * 3 + 1)) {
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 0));
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 1));
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 2));
                            if (c_tab2_copy.get(a * 4 + 3).equals("X")) {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(SignROW.get(SR_count));
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_tab2_copy = change_char_bonds(i_tab2.get(a * 3 + 1), SignROW.get(SR_count), bond_number2, i_tab2, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(clique_siz, i_tab2.get(a * 3 + 1), 2, mapped_atoms);
                                c_bond_neighborsA = change_char_bonds(cor_atom, SignROW.get(SR_count), neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setB.add(i_tab2.get(a * 3 + 0));
                        i_bond_setB.add(i_tab2.get(a * 3 + 1));
                        i_bond_setB.add(i_tab2.get(a * 3 + 2));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
                        set_bondnum_B++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmapped_atoms_molB.get(b).intValue() == i_tab2.get(a * 3 + 1)) {
                    for (int c = 0; c < clique_siz; c++) {
                        if (mapped_atoms.get(c * 2 + 1).intValue() == i_tab2.get(a * 3 + 0)) {
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 0));
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 1));
                            i_bond_neighborsB.add(i_tab2.get(a * 3 + 2));
                            if (c_tab2_copy.get(a * 4 + 2).equals("X")) {
                                c_bond_neighborsB.add(SignROW.get(SR_count));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add("X");
                                c_tab2_copy = change_char_bonds(i_tab2.get(a * 3 + 0), SignROW.get(SR_count), bond_number2, i_tab2, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(clique_siz, i_tab2.get(a * 3 + 0), 2, mapped_atoms);
                                c_bond_neighborsA = change_char_bonds(cor_atom, SignROW.get(SR_count), neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 2));
                                c_bond_neighborsB.add("X");
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setB.add(i_tab2.get(a * 3 + 0));
                        i_bond_setB.add(i_tab2.get(a * 3 + 1));
                        i_bond_setB.add(i_tab2.get(a * 3 + 2));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
                        set_bondnum_B++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (bond_considered) {
                    break;
                }
            }
            bond_considered = false;
        }

//        System.out.println("Bonds of the neighbor set B: ");
//        for (int a = 0; a < neighbor_bondnum_B; a++) {
//            System.out.println(i_bond_neighborsB.get(a * 3 + 0) + " " + i_bond_neighborsB.get(a * 3 + 1) + " " + i_bond_neighborsB.get(a * 3 + 2));
//            System.out.print(c_bond_neighborsB.get(a * 4 + 0) + " " + c_bond_neighborsB.get(a * 4 + 1) + " : " + c_bond_neighborsB.get(a * 4 + 2) + " ");
//            System.out.println(c_bond_neighborsB.get(a * 4 + 3));
//        }
//        System.out.println("");
//        System.out.println("remaining Bond set B: ");
//        for (int a = 0; a < set_bondnum_B; a++) {
//            System.out.println(i_bond_setB.get(a * 3 + 0) + " " + i_bond_setB.get(a * 3 + 1) + " " + i_bond_setB.get(a * 3 + 2));
//            System.out.print(c_bond_setB.get(a * 4 + 0) + " " + c_bond_setB.get(a * 4 + 1) + " : " + c_bond_setB.get(a * 4 + 2) + " ");
//            System.out.println(c_bond_setB.get(a * 4 + 3));
//        }
//        System.out.println("");
//        System.out.println("\"Calling Iterator \"");
        boolean dummy = false;
        Iterator(dummy, mapped_atoms_number, mapped_atoms, neighbor_bondnum_A, neighbor_bondnum_B, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, set_bondnum_A, set_bondnum_B, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);
//        System.out.println("\"DONE Calling Iterator \"");
        return 0;
    }

    private int generate_c_tab1_copy() {

        for (int a = 0; a < bond_number1; a++) {
            c_tab1_copy.add(c_tab1.get(a * 2 + 0));
            c_tab1_copy.add(c_tab1.get(a * 2 + 1));
            c_tab1_copy.add("X");
            c_tab1_copy.add("X");
        }
        return 0;
    }

    private int generate_c_tab2_copy() {

        for (int a = 0; a < bond_number2; a++) {
            c_tab2_copy.add(c_tab2.get(a * 2 + 0));
            c_tab2_copy.add(c_tab2.get(a * 2 + 1));
            c_tab2_copy.add("X");
            c_tab2_copy.add("X");
        }
        return 0;
    }

//Function compaires a structure array with itself. Sometimes a mapping occurs several times within the array.
//The function eliminates these recurring mappings. Function is called in function best_solution.
//The function is called by itself as long as the last list element is processed.
    private List<Integer> remove_recurring_mappings(List<Integer> atom_mapping) {

        boolean exist = true;
        List<Integer> temp_map = new ArrayList<>();
        int temp_counter = 0;
        int atom_mapping_size = atom_mapping.size();
        for (int x = 0; x < atom_mapping_size; x = x + 2) {
            int atom = atom_mapping.get(x);
            for (int y = x + 2; y < atom_mapping_size; y = y + 2) {
                if (atom == atom_mapping.get(y)) {
                    exist = false;
                }
            }
            if (exist == true) {
                temp_map.add(atom_mapping.get(x));
                temp_map.add(atom_mapping.get(x + 1));
                temp_counter = temp_counter + 2;
            }
            exist = true;
        }
        return temp_map;
    }

//In a char field the char sign of corresponding_atom is replaced by new_symbol
    private List<String> change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, List<Integer> i_bond_neighbors, List<String> c_bond_neighbors) {

        List<String> c_bond_neighbors_Local = new ArrayList<>(c_bond_neighbors);

        for (int a = 0; a < neighbor_bondnum; a++) {
            if ((i_bond_neighbors.get(a * 3 + 0).equals(corresponding_atom))
                    && (c_bond_neighbors_Local.get(a * 4 + 2).equals("X"))) {
                c_bond_neighbors_Local.set(a * 4 + 2, c_bond_neighbors_Local.get(a * 4 + 0));
                c_bond_neighbors_Local.set(a * 4 + 0, new_symbol);
            }
            if ((i_bond_neighbors.get(a * 3 + 1).equals(corresponding_atom))
                    && (c_bond_neighbors_Local.get(a * 4 + 3).equals("X"))) {
                c_bond_neighbors_Local.set(a * 4 + 3, c_bond_neighbors_Local.get(a * 4 + 1));
                c_bond_neighbors_Local.set(a * 4 + 1, new_symbol);
            }
        }

        return c_bond_neighbors_Local;
    }

    private int search_corresponding_atom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, List<Integer> mapped_atoms) {

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if (molecule == 1) {
                if (mapped_atoms.get(a * 2 + 0) == (atom_from_other_molecule)) {
                    corresponding_atom = mapped_atoms.get(a * 2 + 1);
                }
            }
            if (molecule == 2) {
                if (mapped_atoms.get(a * 2 + 1) == (atom_from_other_molecule)) {
                    corresponding_atom = mapped_atoms.get(a * 2 + 0);
                }
            }
        }
        return corresponding_atom;
    }

    private int Iterator(boolean MAPPING_check, int mapped_atoms_num,
            List<Integer> mapped_atoms, int neighbor_bondnum_A, int neighbor_bondnum_B,
            List<Integer> i_bond_neighborsA, List<Integer> i_bond_neighborsB,
            List<String> c_bond_neighborsA, List<String> c_bond_neighborsB,
            int set_num_A, int set_num_B, List<Integer> i_bond_setA, List<Integer> i_bond_setB,
            List<String> c_bond_setA, List<String> c_bond_setB) {

        //check possible mappings:
        boolean no_Map = true;
        for (int row = 0; row < neighbor_bondnum_A; row++) {
            String G1A = c_bond_neighborsA.get(row * 4 + 0);
            String G2A = c_bond_neighborsA.get(row * 4 + 1);
            IAtom a1 = this.ac1.getAtom(i_bond_neighborsA.get(row * 3 + 0) - 1);
            IAtom a2 = this.ac1.getAtom(i_bond_neighborsA.get(row * 3 + 1) - 1);
            IBond bond1 = this.ac1.getBond(a1, a2);

            for (int column = 0; column < neighbor_bondnum_B; column++) {
//                System.out.println("c_bond_neighborsA  " + c_bond_neighborsA);
//                System.out.println("i_bond_neighborsA  " + i_bond_neighborsA);
//                System.out.println("neighbor_bondnum_A  " + neighbor_bondnum_A);
//                
//                System.out.println("c_bond_neighborsB  " + c_bond_neighborsB);
//                System.out.println("i_bond_neighborsB  " + i_bond_neighborsB);
//                System.out.println("neighbor_bondnum_B  " + neighbor_bondnum_B);

                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);

                IAtom b1 = this.ac2.getAtom(i_bond_neighborsB.get(column * 3 + 0) - 1);
                IAtom b2 = this.ac2.getAtom(i_bond_neighborsB.get(column * 3 + 1) - 1);
                IBond bond2 = this.ac2.getBond(b1, b2);

                /*
                 * Check if bond matching also possible
                 */
                boolean flag = isMatchFeasible(bond1, bond2, shouldMatchBonds, shouldMatchRings, matchAtomType);

                if ((G1A.equals(G1B)) && (G2A.equals(G2B)) && flag) {
                    no_Map = false;
                    break;
                } else if ((G1A.equals(G2B)) && (G2A.equals(G1B)) && flag) {
                    no_Map = false;
                    break;
                }

//                if (flag) {
//                    System.out.println("flag " + flag + ", no_Map " + no_Map);
//                    System.out.println("bond1 " + bond1.getAtom(0).getSymbol());
//                    System.out.println("bond1 " + bond1.getAtom(1).getSymbol());
//
//                    System.out.println("bond2 " + bond1.getAtom(0).getSymbol());
//                    System.out.println("bond2 " + bond1.getAtom(1).getSymbol());
//
//                    System.out.println(bond1.getOrder().numeric() + "," + bond1.getOrder().numeric());
//                }
            }
            if (!no_Map) {
                break;
            }
        }

//        System.out.println("c_bond_neighborsA - check before " + c_bond_neighborsA.size());
        if ((neighbor_bondnum_A == 0) || (neighbor_bondnum_B == 0) || (MAPPING_check) || (no_Map)) { //MAPPING_check <=> no_further_MAPPINGS

            //solution mappings are pushed in list final_MAPPINGS
            if ((mapped_atoms_num) >= best_MAPPING_size) {
                if ((mapped_atoms_num) > best_MAPPING_size) {
                    getFinalMappings().clear();
                    best_MAPPING_size = mapped_atoms_num;
                }
                getFinalMappings().add(mapped_atoms);
            }
            return 0;
        }

//deletion of global vectors
        i_globalA.clear();
        i_globalB.clear();
        c_globalA.clear();
        c_globalB.clear();

        //redefining of global vectors and variables
        nNum_globalA = neighbor_bondnum_A;
        nNum_globalB = neighbor_bondnum_B;
        i_globalA.addAll(i_bond_neighborsA);
        i_globalB.addAll(i_bond_neighborsB);
        c_globalA.addAll(c_bond_neighborsA);
        c_globalB.addAll(c_bond_neighborsB);
        this.MARCS.clear();
        this.MARCS = new ArrayList<>(neighbor_bondnum_A * neighbor_bondnum_B);
        for (int i = 0; i < neighbor_bondnum_A * neighbor_bondnum_B; i++) {
            MARCS.add(i, 0);
        }
        for (int row = 0; row < neighbor_bondnum_A; row++) {
            for (int column = 0; column < neighbor_bondnum_B; column++) {

                String G1A = c_bond_neighborsA.get(row * 4 + 0);
                String G2A = c_bond_neighborsA.get(row * 4 + 1);
                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);

                if (((G1A.equals(G1B)) && (G2A.equals(G2B))) || ((G1A.equals(G2B)) && (G2A.equals(G1B)))) {
                    MARCS.set(row * neighbor_bondnum_B + column, 1);
//                    System.out.println("Atoms: " + G1A + " " + G2A + " " + G1B + " " + G2B);
                }
            }
        }

        //Initialization of the tree structure which is needed
        // for the identification of redundant matrices
        first = last = new BinaryTree();
        last.setValue(-1);
        last.equal = null;
        last.not_equal = null;

        bestarcsleft = 0;
        startsearch();

        Stack<List<Integer>> BESTARCS_copy = (Stack<List<Integer>>) BESTARCS.clone();

        while (!BESTARCS.empty()) {
            BESTARCS.pop();
        }

        while (!BESTARCS_copy.empty()) {
            List<Integer> MARCS_vector = BESTARCS_copy.peek();
//            print_matrix(MARCS_vector, neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA, neighbor_bondnum_B, i_bond_neighborsB, c_bond_neighborsB);
            List<Integer> new_MAPPING = find_mcgregor_MAPPING(MARCS_vector, mapped_atoms_num, mapped_atoms, neighbor_bondnum_A, i_bond_neighborsA, neighbor_bondnum_B, i_bond_neighborsB);

            int new_MAPPING_size = new_MAPPING.size();
            boolean no_further_MAPPINGS = false;
            if (mapped_atoms_num == new_MAPPING_size / 2) {
                no_further_MAPPINGS = true;
            }
            //new values for neighbor_bondnum_A + neighbor_bondnum_B
            //new arrays for i_bond_neighborsA + i_bond_neighborsB + c_bond_neighborsA + c_bond_neighborsB
            int new_neighbor_numA = 0; //instead of neighbor_bondnum_A
            int new_neighbor_numB = 0; //instead of neighbor_bondnum_B
            List<Integer> new_i_neighborsA = new ArrayList<>(); //instead of i_bond_neighborsA
            List<Integer> new_i_neighborsB = new ArrayList<>(); //instead of i_bond_neighborsB
            List<String> new_c_neighborsA = new ArrayList<>(); //instead of c_bond_neighborsA
            List<String> new_c_neighborsB = new ArrayList<>(); //instead of c_bond_neighborsB
            new_i_neighborsA.clear();
            new_i_neighborsB.clear();
            new_c_neighborsA.clear();
            new_c_neighborsB.clear();

            //new values for set_num_A + set_num_B
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB
            int set_bondnum_A = 0; //instead of set_num_A
            int set_bondnum_B = 0; //instead of set_num_B
            List<Integer> new_i_bond_setA = new ArrayList<>(); //instead of i_bond_setA
            List<Integer> new_i_bond_setB = new ArrayList<>(); //instead of i_bond_setB
            List<String> new_c_bond_setA = new ArrayList<>(); //instead of c_bond_setA
            List<String> new_c_bond_setB = new ArrayList<>(); //instead of c_bond_setB
            new_i_bond_setA.clear();
            new_i_bond_setB.clear();
            new_c_bond_setA.clear();
            new_c_bond_setB.clear();

            List<String> c_setB_copy = generate_c_setB_copy(set_num_B, c_bond_setB);
            List<String> c_setA_copy = new ArrayList<>(c_bond_setA);
            //find unmapped atoms of molecule A
            List<Integer> unmapped_atoms_molA = new ArrayList<>();
            unmapped_atoms_molA.clear();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;
            for (int a = 1; a <= atom_num_H_1; a++) {
                for (int b = 0; b < (new_MAPPING_size / 2); b++) {
                    if (a == new_MAPPING.get(b * 2 + 0)) {
                        atomA_is_unmapped = false;
                    }
                }
                if (atomA_is_unmapped == true) {
                    unmapped_atoms_molA.add(a);
                    unmapped_numA++;
                }
                atomA_is_unmapped = true;
            }
            //The special signs must be transfered to the corresponding atoms of molecule B
            int SR_count = 0;
            boolean bond_considered = false;
            boolean normal_bond = true;
            for (int a = 0; a < set_num_A; a++) {
                for (int b = 0; b < unmapped_numA; b++) {
                    if (unmapped_atoms_molA.get(b).intValue() == i_bond_setA.get(a * 3 + 0)) {
                        for (int c = 0; c < (new_MAPPING_size / 2); c++) {
                            if (new_MAPPING.get(c * 2 + 0).intValue() == i_bond_setA.get(a * 3 + 1)) {
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                if (c_setA_copy.get(a * 4 + 3).equals("X")) {
                                    new_c_neighborsA.add(SignROW.get(SR_count));
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    c_setA_copy = change_char_bonds(i_bond_setA.get(a * 3 + 1), SignROW.get(SR_count), set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = search_corresponding_atom((new_MAPPING_size / 2), i_bond_setA.get(a * 3 + 1), 1, new_MAPPING);
                                    c_setB_copy = change_char_bonds(cor_atom, SignROW.get(SR_count), set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;
                                } else {
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 3));
                                }
                                normal_bond = false;
                                new_neighbor_numA++;
                            }
                        }
                        if (normal_bond) {
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;
                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unmapped_atoms_molA.get(b).intValue() == i_bond_setA.get(a * 3 + 1)) {
                        for (int c = 0; c < (new_MAPPING_size / 2); c++) {
                            if (new_MAPPING.get(c * 2 + 0).intValue() == i_bond_setA.get(a * 3 + 0)) {
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                if (c_setA_copy.get(a * 4 + 2).equals("X")) {
                                    new_c_neighborsA.add(SignROW.get(SR_count));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add("X");
                                    c_setA_copy = change_char_bonds(i_bond_setA.get(a * 3 + 0), SignROW.get(SR_count), set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = search_corresponding_atom((new_MAPPING_size / 2), i_bond_setA.get(a * 3 + 0), 1, new_MAPPING);
                                    c_setB_copy = change_char_bonds(cor_atom, SignROW.get(SR_count), set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;
                                } else {
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 2));
                                    new_c_neighborsA.add("X");
                                }
                                normal_bond = false;
                                new_neighbor_numA++;
                            }
                        }
                        if (normal_bond) {
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;
                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (bond_considered) {
                        break;
                    }
                }
                bond_considered = false;
            }
            //find unmapped atoms of molecule B
            List<Integer> unmapped_atoms_molB = new ArrayList<>();
            unmapped_atoms_molB.clear();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;
            for (int a = 1; a <= atom_num_H_2; a++) {
                for (int b = 0; b < (new_MAPPING_size / 2); b++) {
                    if (a == new_MAPPING.get(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }
                }
                if (atomB_is_unmapped == true) {
                    unmapped_atoms_molB.add(a);
                    unmapped_numB++;
                }
                atomB_is_unmapped = true;
            }

            //The special signs must be transfered to the corresponding atoms of molecule A
            bond_considered = false;
            normal_bond = true;
            for (int a = 0; a < set_num_B; a++) {
                for (int b = 0; b < unmapped_numB; b++) {
                    if (unmapped_atoms_molB.get(b).intValue() == i_bond_setB.get(a * 3 + 0)) {
                        for (int c = 0; c < (new_MAPPING_size / 2); c++) {
                            if (new_MAPPING.get(c * 2 + 1).intValue() == i_bond_setB.get(a * 3 + 1)) {
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));
                                if (c_setB_copy.get(a * 4 + 3).equals("X")) {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(SignROW.get(SR_count));
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    c_setB_copy = change_char_bonds(i_bond_setB.get(a * 3 + 1), SignROW.get(SR_count), set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = search_corresponding_atom((new_MAPPING_size / 2), i_bond_setB.get(a * 3 + 1), 2, new_MAPPING);
                                    new_c_neighborsA = change_char_bonds(cor_atom, SignROW.get(SR_count), new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;
                                } else {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 3));
                                }
                                normal_bond = false;
                                new_neighbor_numB++;
                            }
                        }
                        if (normal_bond) {
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;
                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unmapped_atoms_molB.get(b).intValue() == i_bond_setB.get(a * 3 + 1)) {
                        for (int c = 0; c < (new_MAPPING_size / 2); c++) {
                            if (new_MAPPING.get(c * 2 + 1).equals(i_bond_setB.get(a * 3 + 0))) {
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));
                                if (c_setB_copy.get(a * 4 + 2).equals("X")) {
                                    new_c_neighborsB.add(SignROW.get(SR_count));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add("X");
                                    c_setB_copy = change_char_bonds(i_bond_setB.get(a * 3 + 0), SignROW.get(SR_count), set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = search_corresponding_atom((new_MAPPING_size / 2), i_bond_setB.get(a * 3 + 0), 2, new_MAPPING);
                                    new_c_neighborsA = change_char_bonds(cor_atom, SignROW.get(SR_count), new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;
                                } else {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 2));
                                    new_c_neighborsB.add("X");
                                }
                                normal_bond = false;
                                new_neighbor_numB++;
                            }
                        }
                        if (normal_bond) {
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;
                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (bond_considered) {
                        break;
                    }
                }
                bond_considered = false;
            }
            Iterator(no_further_MAPPINGS, (new_MAPPING_size / 2), new_MAPPING, new_neighbor_numA, new_neighbor_numB,
                    new_i_neighborsA, new_i_neighborsB, new_c_neighborsA, new_c_neighborsB,
                    set_bondnum_A, set_bondnum_B, new_i_bond_setA, new_i_bond_setB, new_c_bond_setA, new_c_bond_setB);

            if (!BESTARCS_copy.isEmpty()) {
                BESTARCS_copy.pop();
            }
        }
        return 0;
    }

    private void startsearch() {

//        System.out.println("FIXARCS neighbor_bondnum_A " + nNum_globalA);
//        System.out.println("FIXARCS neighbor_bondnum_B " + nNum_globalB);
        this.FIXARCS = new ArrayList<>(nNum_globalA * nNum_globalB);
        for (int i = 0; i < nNum_globalA * nNum_globalB; i++) {
            FIXARCS.add(i, 0);
        }

        int x = 0;
        int y = 0;
        while ((x < nNum_globalA) && (MARCS.get(x * nNum_globalB + y) != 1)) {
            y++;
            if (y == nNum_globalB) {
                y = 0;
                x++;
            }
        }
        if (x == nNum_globalA) {
            y = nNum_globalB - 1;
            x = x - 1;
        }

        if (MARCS.get(x * nNum_globalB + y) == 0) {
            partsearch(x, y, MARCS);
        }
        if (MARCS.get(x * nNum_globalB + y) != 0) {
            partsearch(x, y, MARCS);
            MARCS.set(x * nNum_globalB + y, 0);
            partsearch(x, y, MARCS);
        }
    }

    private void partsearch(int xstart, int ystart, List<Integer> TEMPMARCS) {

//        System.out.println("partsearch TEMPMARCS " + TEMPMARCS);
        int x = xstart;
        int y = ystart;

//        System.out.println("X " + x + ", Y " + y);
//        System.out.println("nNum_globalA " + nNum_globalA + ", nNum_globalB " + nNum_globalB);
        if (TEMPMARCS.get(xstart * nNum_globalB + ystart) == 1) {
            TEMPMARCS = remove_redundant_arcs(xstart, ystart, TEMPMARCS);

            int arcsleft = 0;
            for (int a = 0; a < nNum_globalA; a++) {
                for (int b = 0; b < nNum_globalB; b++) {
                    if (TEMPMARCS.get(a * nNum_globalB + b) == 1) {
                        arcsleft++;
                    }
                }
            }

            //test Bestarcsleft and skip rest if needed
//            System.out.println("arcsleft " + arcsleft + " bestarcsleft " + bestarcsleft);
            if (arcsleft >= bestarcsleft) {
                do {
                    y++;
                    if (y == nNum_globalB) {
                        y = 0;
                        x++;
                    }
//                    System.out.println("x * nNum_globalB + y " + (x * nNum_globalB + y));

                } while (x < nNum_globalA && TEMPMARCS.get(x * nNum_globalB + y) != 1);
                if (x < nNum_globalA) {
                    partsearch(x, y, TEMPMARCS);
                    TEMPMARCS.set(x * nNum_globalB + y, 0);
                    partsearch(x, y, TEMPMARCS);
                } else {
                    if (arcsleft > bestarcsleft) {
                        BinaryTree.remove_tree_structure(first);
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (check_MARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }
                }
            }
        } else {
            do {
                y++;
                if (y == nNum_globalB) {
                    y = 0;
                    x++;
                }
            } while (x < nNum_globalA && TEMPMARCS.get(x * nNum_globalB + y) != 1);
            if (x < nNum_globalA) {
                partsearch(x, y, TEMPMARCS);
                TEMPMARCS.set(x * nNum_globalB + y, 0);
                partsearch(x, y, TEMPMARCS);
            } else {
                int arcsleft = 0;
                for (int a = 0; a < nNum_globalA; a++) {
                    for (int b = 0; b < nNum_globalB; b++) {
                        if (TEMPMARCS.get(a * nNum_globalB + b) == 1) {
                            arcsleft++;
                        }
                    }
                }
                if (arcsleft >= bestarcsleft) {
                    if (arcsleft > bestarcsleft) {
                        remove_tree_structure(first);
                        /* TO DO Asad*/
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (check_MARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }
                }
            }
        }
    }

    private List<String> generate_c_setB_copy(int bond_number, List<String> c_setB) {

        List<String> c_setB_copy = new ArrayList<>();
        for (int a = 0; a < bond_number; a++) {
            c_setB_copy.add(c_setB.get(a * 4 + 0));
            c_setB_copy.add(c_setB.get(a * 4 + 1));
            c_setB_copy.add("X");
            c_setB_copy.add("X");
        }
        return c_setB_copy;
    }

//The function is called in function partsearch. The function is given a temporary matrix and a position (row/column)
//within this matrix. First the function sets all entries to zero, which can be exlcuded in respect to the current
//atom by atom matching. After this the function replaces all entries in the same row and column of the current
//position by zeros. Only the entry of the current position is set to one.
//Return value "count_arcsleft" counts the number of arcs, which are still in the matrix.
    private List<Integer> remove_redundant_arcs(int row, int column, List<Integer> MARCS) {

        List<Integer> MARCS_LOCAL = new ArrayList<>(MARCS);
        int G1_atom = i_globalA.get(row * 3 + 0);
        int G2_atom = i_globalA.get(row * 3 + 1);
        int G3_atom = i_globalB.get(column * 3 + 0);
        int G4_atom = i_globalB.get(column * 3 + 1);
        for (int x = 0; x < nNum_globalA; x++) {
            int row_atom1 = i_globalA.get(x * 3 + 0);
            int row_atom2 = i_globalA.get(x * 3 + 1);
            for (int y = 0; y < nNum_globalB; y++) {
                int column_atom3 = i_globalB.get(y * 3 + 0);
                int column_atom4 = i_globalB.get(y * 3 + 1);
                if (((G1_atom == row_atom1) || (G1_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
                    MARCS_LOCAL.set(x * nNum_globalB + y, 0);
                }
                if (((G2_atom == row_atom1) || (G2_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
                    MARCS_LOCAL.set(x * nNum_globalB + y, 0);
                }
                if (((G3_atom == column_atom3) || (G3_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS_LOCAL.set(x * nNum_globalB + y, 0);
                }
                if (((G4_atom == column_atom3) || (G4_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS_LOCAL.set(x * nNum_globalB + y, 0);
                }
            }
        }

        for (int v = 0; v < nNum_globalA; v++) {
            MARCS_LOCAL.set(v * nNum_globalB + column, 0);
        }
        for (int w = 0; w < nNum_globalB; w++) {
            MARCS_LOCAL.set(row * nNum_globalB + w, 0);
        }
        MARCS_LOCAL.set(row * nNum_globalB + column, 1);
        return MARCS_LOCAL;
    }

    /*
     * The function is called in function partsearch. The function is given a temporary matrix.
     * The function checks whether the temporary matrix is already found by calling the function
     *"verify_nodes". If the matrix already exists the function returns false which means that
     * the matrix will not be stored. Otherwise the function returns true which means that the
     * matrix will be stored in function partsearch.
     */
    private boolean check_MARCS(List<Integer> MARCS) {

        List<Integer> posnum_list = new ArrayList<>(nNum_globalA * nNum_globalA);
        for (int i = 0; i < nNum_globalA * nNum_globalA; i++) {
            posnum_list.add(i, 0);
        }
        int y = 0;
        int count_entries = 0;
        for (int x = 0; x < (nNum_globalA * nNum_globalB); x++) {
            if (MARCS.get(x) == 1) {
                posnum_list.set(y++, x);
                count_entries++;
            }
        }
        verify_nodes(posnum_list, first, 0, count_entries);
        return new_matrix;
    }

    private boolean verify_nodes(List<Integer> matrix, BinaryTree cur_struc, int x, int field_length) {

        if ((matrix.get(x).equals(cur_struc.getValue()) && (x < field_length))) {
            if (cur_struc.equal != null) {
                new_matrix = false;
                verify_nodes(matrix, cur_struc.equal, x + 1, field_length);
            }
        }
        if (matrix.get(x) != cur_struc.getValue()) {
            if (cur_struc.not_equal != null) {
                verify_nodes(matrix, cur_struc.not_equal, x, field_length);
            }
            if (cur_struc.not_equal == null) {
                cur_struc.not_equal = new BinaryTree();
                cur_struc.not_equal.setValue(matrix.get(x));
                cur_struc.not_equal.not_equal = null;
                int y = 0;
                BinaryTree last_one = cur_struc.not_equal;

                while ((y + x + 1) < field_length) {
                    last_one.equal = new BinaryTree();
                    last_one = last_one.equal;
                    last_one.setValue(matrix.get(y + x + 1));
                    last_one.not_equal = null;
                    y++;
                }
                last_one.equal = null;
                new_matrix = true;
            }
        }
        return true;
    }

    /**
     * @return the final_MAPPINGS
     */
    public List<List<Integer>> getFinalMappings() {
        return final_MAPPINGS;
    }

}
