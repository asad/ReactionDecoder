/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus1;

import org.openscience.smsd.tools.Utility;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.graph.Edge;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CompatibilityGraph {

    private final AtomMatcher atomMatcher;
    private final BondMatcher bondMatcher;

    /**
     * @return the c_edges
     */
    public List<Edge> getCEdges() {
        return global_c_edges;
    }

    /**
     * @return the d_edges
     */
    public List<Edge> getDEdges() {
        return global_d_edges;
    }

    /**
     * @return the comp_graph_nodes
     */
    public List<Integer> getCompGraphNodes() {
        return comp_graph_nodes;
    }

    private final List<Edge> global_c_edges;
    private final List<Edge> global_d_edges;
    private final List<Integer> comp_graph_nodes;
    private final List<Integer> comp_graph_nodes_C_zero;

    private final int atom_num_H_1;
    private final int atom_num_H_2;

    private final int bond_number1;
    private final int bond_number2;

    private final Map<String, Integer> SYMBOL_VALUE;

    private final List<Integer> i_tab1;
    private final List<Integer> i_tab2;

    private final List<String> c_tab1;
    private final List<String> c_tab2;

    private final IAtomContainer ac1;
    private final IAtomContainer ac2;

    private final List<IAtom> atomstr1;
    private final List<IAtom> atomstr2;

    private final boolean DEBUG = false;

    /**
     * Creates a new instance of SearchCliques
     *
     *
     * @param f1
     * @param f2
     */
    public CompatibilityGraph(IAtomContainer f1, IAtomContainer f2,
            AtomMatcher am, BondMatcher bm) {

        this.global_c_edges = new ArrayList<>();//Initialize the c_edges Vector
        this.global_d_edges = new ArrayList<>();//Initialize the d_edges Vector
        this.atomMatcher = am;
        this.bondMatcher = bm;

        this.SYMBOL_VALUE = new TreeMap<>();

        MoleculeHandler file1 = new MoleculeHandler(f1, false);
        MoleculeHandler file2 = new MoleculeHandler(f2, false);

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
        this.comp_graph_nodes_C_zero = new ArrayList<>();//Initialize the comp_graph_nodes_C_zero Vector
    }

    public int searchCliques() {

        generate_compatibility_graph_nodes();
        generate_compatibility_graph();
        if (DEBUG) {
            System.out.println("c_edges_size " + getCEdges().size());
            System.out.println("bond count: " + ac1.getBondCount());
            System.out.println("bond count: " + ac2.getBondCount());
        }
        if (getCEdges().isEmpty()) {

            if (DEBUG) {
                System.out.println("Switching to complex mode ");
            }
            getCompGraphNodes().clear();
            getCEdges().clear();
            getDEdges().clear();
            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            generate_compatibility_graph_if_C_edge_number_is_zero();
            comp_graph_nodes_C_zero.clear();
        }

        /*
         * Transfor C and D edges from Edge to Integer
         */
        List<Edge> unique_global_c_edges = new ArrayList<>(new HashSet<>(getCEdges()));//remove any duplicates;
        List<Edge> unique_global_d_edges = new ArrayList<>(new HashSet<>(getDEdges()));//remove any duplicates;

        if (DEBUG) {
            System.out.println("**************************************************");
            System.out.println("--MCS PLUS--");
            System.out.println("C_edges: " + unique_global_c_edges.size());
            System.out.println("D_edges: " + unique_global_d_edges.size());
            System.out.println("comp_graph_nodes: " + getCompGraphNodes().size());
        }
        return getCompGraphNodes().size();
    }

    private List<List<Integer>> label_atoms(List<Integer> basic_atom_vector, int bond_num, List<IAtom> atoms, List<Integer> i_tab, List<String> c_tab) {

        ArrayList<List<Integer>> label_list = new ArrayList<>();

//        if (DEBUG) {
//            System.out.println("Vector Atom Str: ");
//            for (int b = 0; b < atoms.size(); b++) {
//                System.err.print(atoms.get(b).getSymbol() + ",");
//            }
//            System.LOGGER.debug();
//            System.LOGGER.debug("basic_atom_vector");
//            for (int b = 0; b < basic_atom_vector.size(); b++) {
//                System.err.print(basic_atom_vector.get(b) + ",");
//            }
//            System.LOGGER.debug();
//            System.LOGGER.debug("i_tab");
//            for (int b = 0; b < i_tab.size(); b++) {
//                System.err.print(i_tab.get(b) + ",");
//            }
//            System.LOGGER.debug();
//            System.LOGGER.debug("c_tab");
//            for (int b = 0; b < c_tab.size(); b++) {
//                System.err.print(c_tab.get(b) + ",");
//            }
//            System.LOGGER.debug();
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
            String atom1_type = atom1.getSymbol(); //+ atom1.getAtomicNumber();

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
                    /*Get neighbour Atom*/
                    IAtom atom2 = atoms.get(i_tab.get(b * 3 + 1) - 1);
                    String atom2_type = c_tab.get(b * 2 + 1);// + atom2.getAtomicNumber();

                    //System.out.println("atom2_type " + atom2_type + ", atom2 " + atom2.getSymbol());
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
                    /*Get neighbour Atom*/
                    IAtom atom2 = atoms.get(i_tab.get(b * 3 + 0) - 1);
                    String atom2_type = c_tab.get(b * 2 + 0);// + atom2.getAtomicNumber();

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

//        System.out.println("atomstr1 " + atomstr1);
//        System.out.println("atomstr2 " + atomstr2);
        List<Integer> basic_atom_vec_A = reduce_atomset(atom_num_H_1, bond_number1, atomstr1, i_tab1, c_tab1);
        List<Integer> basic_atom_vec_B = reduce_atomset(atom_num_H_2, bond_number2, atomstr2, i_tab2, c_tab2);

        List<List<Integer>> label_list_molA = label_atoms(basic_atom_vec_A, bond_number1, atomstr1, i_tab1, c_tab1);
        List<List<Integer>> label_list_molB = label_atoms(basic_atom_vec_B, bond_number2, atomstr2, i_tab2, c_tab2);

        int molA_nodes = 0;
        int count_nodes = 1;

        if (DEBUG) {
            System.out.println("basic_atom_vec_A " + basic_atom_vec_A);
            System.out.println("basic_atom_vec_B " + basic_atom_vec_B);

            System.out.println("label_list_molA " + label_list_molA);
            System.out.println("label_list_molB " + label_list_molB);
        }
        for (List<Integer> labelA : label_list_molA) {
            int molB_nodes = 0;
            for (List<Integer> labelB : label_list_molB) {
                if (labelA.equals(labelB)) {
                    getCompGraphNodes().add(basic_atom_vec_A.get(molA_nodes));
                    getCompGraphNodes().add(basic_atom_vec_B.get(molB_nodes));
                    getCompGraphNodes().add(count_nodes++);
                    if (DEBUG) {
                        System.out.println("labelA " + labelA + ", labelB " + labelB + ", count_nodes " + count_nodes + "\n");
                    }

                }
                molB_nodes++;
            }
            molA_nodes++;
        }

        if (DEBUG) {
            System.out.println("generate_compatibility_graph_nodes comp_graph_nodes: " + getCompGraphNodes().size());
        }

        return 0;
    }

    private int generate_compatibility_graph() {

        int vector_size = getCompGraphNodes().size();

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
                        if ((getCompGraphNodes().get(a).equals(i_tab1.get(x * 3 + 0))
                                && getCompGraphNodes().get(b).equals(i_tab1.get(x * 3 + 1)))) {

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(a) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(b) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                                System.out.println("BOND " + i_tab1.get(x * 3 + 2));
//                            }
                            IAtom a1 = this.ac1.getAtom(getCompGraphNodes().get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(getCompGraphNodes().get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);
                            molecule1_pair_connected = true;
                            if (bond1 != null) {
                                break;
                            }
                        } else if ((getCompGraphNodes().get(a).equals(i_tab1.get(x * 3 + 1))
                                && getCompGraphNodes().get(b).equals(i_tab1.get(x * 3 + 0)))) {

//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(a) + ", i_tab1.get(x * 3 + 1) " + i_tab1.get(x * 3 + 1));
//                                System.out.println("comp_graph_nodes.get(a) " + comp_graph_nodes.get(b) + ", i_tab1.get(x * 3 + 0) " + i_tab1.get(x * 3 + 0));
//                                System.out.println("BOND " + i_tab1.get(x * 3 + 2));
//                            }
                            IAtom a1 = this.ac1.getAtom(getCompGraphNodes().get(a) - 1);
                            IAtom a2 = this.ac1.getAtom(getCompGraphNodes().get(b) - 1);
                            bond1 = this.ac1.getBond(a1, a2);
                            molecule1_pair_connected = true;
                            if (bond1 != null) {
                                break;
                            }
                        }
                    }
                    //exists a bond in molecule 2, so that molecule 2 pair is connected?
                    for (int y = 0; y < bond_number2; y++) {
                        if ((getCompGraphNodes().get(a + 1).equals(i_tab2.get(y * 3 + 0))
                                && getCompGraphNodes().get(b + 1).equals(i_tab2.get(y * 3 + 1)))) {
//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(a + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(b + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                                System.out.println("BOND " + i_tab2.get(y * 3 + 2));
//                            }
                            IAtom a1 = this.ac2.getAtom(getCompGraphNodes().get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(getCompGraphNodes().get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);
                            molecule2_pair_connected = true;
                            if (bond2 != null) {
                                break;
                            }

                        } else if ((getCompGraphNodes().get(a + 1).equals(i_tab2.get(y * 3 + 1))
                                && getCompGraphNodes().get(b + 1).equals(i_tab2.get(y * 3 + 0)))) {
//                            if (DEBUG) {
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(a + 1) + ", i_tab2.get(x * 3 + 1) " + i_tab2.get(y * 3 + 1));
//                                System.out.println("comp_graph_nodes.get(a+1) " + comp_graph_nodes.get(b + 1) + ", i_tab2.get(x * 3 + 0) " + i_tab2.get(y * 3 + 0));
//                                System.out.println("BOND " + i_tab2.get(y * 3 + 2));
//                            }
                            IAtom a1 = this.ac2.getAtom(getCompGraphNodes().get(a + 1) - 1);
                            IAtom a2 = this.ac2.getAtom(getCompGraphNodes().get(b + 1) - 1);
                            bond2 = this.ac2.getBond(a1, a2);
                            molecule2_pair_connected = true;
                            if (bond2 != null) {
                                break;
                            }

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
                            && AtomBondMatcher.matchAtomAndBond(bond1, bond2, atomMatcher, bondMatcher, true)) {
                        matchBondFlag = true;
                    }

                    //in case that both molecule pairs are connected a c-edge is generated
                    if (connectedFlag && matchBondFlag) {
                        Edge edge = new Edge(((a / 3) + 1), ((b / 3) + 1));
                        getCEdges().add(edge);
                    }

                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (disConnectedFlag) {
                        Edge edge = new Edge(((a / 3) + 1), ((b / 3) + 1));
                        getDEdges().add(edge);
                    }

//                    //in case that both molecule pairs are not connected a d-edge is generated
//                    if (connectedFlag && !matchBondFlag) {
//                        Edge edge = new Edge(((a / 3) + 1), ((b / 3) + 1));
//                        getDEdges().add(edge);
//                    }
                }
            }
        }

        if (DEBUG) {
            //print R and Q edges of the compatibility graph
            int c_edges_size = getCEdges().size();
            int d_edges_size = getDEdges().size();

            System.out.println("generate_compatibility_graph C_edges_size " + c_edges_size);
            System.out.println("generate_compatibility_graph D_edges_size " + d_edges_size);
        }

        return 0;
    }

//comp_graph_nodes_C_zero is used to build up of the edges of the compatibility graph
    private int generate_compatibility_graph_nodes_if_C_edge_number_is_zero() {

        if (DEBUG) {
            for (int a = 0; a < atom_num_H_1; a++) {
                String atom1_type = atomstr1.get(a).getSymbol();
                System.out.println("atom1_type " + atom1_type + "(" + (a + 1) + ")");
            }

            for (int b = 0; b < atom_num_H_2; b++) {
                String atom2_type = atomstr2.get(b).getSymbol();
                System.out.println("atom2_type " + atom2_type + "(" + (b + 1) + ")");
            }
        }

        int count_nodes = 1;

        for (int a = 0; a < atom_num_H_1; a++) {
            String atom1_type = atomstr1.get(a).getSymbol();
            int value = atomstr1.get(a).getAtomicNumber() == null
                    ? atomstr1.get(a).hashCode() + 1000 : atomstr1.get(a).getAtomicNumber() + 1000;

            SYMBOL_VALUE.put(atom1_type, value);

            for (int b = 0; b < atom_num_H_2; b++) {
                String atom2_type = atomstr2.get(b).getSymbol();

                if ((atom1_type.equals(atom2_type))) {
                    comp_graph_nodes_C_zero.add(a + 1);
                    comp_graph_nodes_C_zero.add(b + 1);
                    comp_graph_nodes_C_zero.add(SYMBOL_VALUE.get(atom1_type)); //C is label 1
                    comp_graph_nodes_C_zero.add(count_nodes);

                    getCompGraphNodes().add(a + 1);
                    getCompGraphNodes().add(b + 1);
                    getCompGraphNodes().add(count_nodes);

                    if (DEBUG) {
                        System.out.println("a + 1 " + (a + 1));
                        System.out.println("b + 1 " + (b + 1));
                        System.out.println("atoms " + (atom1_type) + "=" + atom2_type);
                        System.out.println("count_nodes " + (count_nodes));
                    }

                    count_nodes++;
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
                            && AtomBondMatcher.matchAtomAndBond(bond1, bond2, atomMatcher, bondMatcher, true)) {
                        matchBondFlag = true;
                    }

//                    if (DEBUG) {
//                        System.out.println("matchbondFlag " + connectedFlag);
//                    }
                    //in case that both molecule pairs are connected a c-edge is generated
                    if (connectedFlag && matchBondFlag) {
                        Edge edge = new Edge(((a / 4) + 1), ((b / 4) + 1));
                        getCEdges().add(edge);
                    }
//
                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (disConnectedFlag) {
                        Edge edge = new Edge(((a / 4) + 1), ((b / 4) + 1));
                        getDEdges().add(edge);
                    }

//                    //in case that both molecule pairs are not connected a d-edge is generated
//                    if (connectedFlag && !matchBondFlag) {
//                        Edge edge = new Edge(((a / 4) + 1), ((b / 4) + 1));
//                        getDEdges().add(edge);
//                    }
                }
            }
        }

        if (DEBUG) {
            //print R and Q edges of the compatibility graph
            int c_edges_size = getCEdges().size();
            int d_edges_size = getDEdges().size();

            System.out.println("C_edges_size " + c_edges_size);
            System.out.println("D_edges_size " + d_edges_size);
        }
        return 0;
    }

    /**
     * Clear data
     */
    public void clear() {
        this.getCompGraphNodes().clear();
        this.comp_graph_nodes_C_zero.clear();
        this.c_tab1.clear();
        this.c_tab2.clear();
        this.getCEdges().clear();
        this.getDEdges().clear();
    }

}
