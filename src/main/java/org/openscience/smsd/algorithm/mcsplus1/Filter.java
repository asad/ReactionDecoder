/* Copyright (C) 2009-2018  Syed Asad Rahman <asad@ebi.ac.uk>
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
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Filter extends McGregor {

    /**
     *
     * @param f1
     * @param f2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public Filter(IAtomContainer f1,IAtomContainer f2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        super(f1, f2,  shouldMatchBonds,  shouldMatchRings,  matchAtomType);
    }

    int postfilter() {
        if ((best_MAPPING_size == 0) && (best_clique_size != 0)) {
            java.util.Iterator<List<Integer>> iter = getFinalMappings().iterator();
            List<Integer> vec = iter.next();
            best_MAPPING_size = vec.size() / 2;
        }
        if ((best_MAPPING_size == 0) && (best_clique_size == 0)) {
            return 0;
        }

        //1. Searching for carboxyl groups 
        //find carboxyl groups of molecule 1
        List<Integer> C_index_A = new ArrayList<>();
        List<List<Integer>> carb_vec_A = new ArrayList<>();
        for (int a = 0; a < atom_num_H_1; a++) {
            if ("C".equals(atomstr1.get(a).getSymbol())) {
                int O_num = 0;
                List<Integer> carboxy_C_vec = new ArrayList<>();
                boolean c_group = true;
                for (int b = 0; b < bond_number1; b++) {
                    if ((a + 1 == i_tab1.get(b * 3 + 0)) && ("O".equals(c_tab1.get(b * 2 + 1)))) {
                        if (check(i_tab1.get(b * 3 + 1), 1)) {
                            carboxy_C_vec.add(i_tab1.get(b * 3 + 1));
                            carboxy_C_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        } else {
                            c_group = false;
                        }
                    }
                    if ((a + 1 == i_tab1.get(b * 3 + 1)) && ("O".equals(c_tab1.get(b * 2 + 0)))) {
                        if (check(i_tab1.get(b * 3 + 0), 1)) {
                            carboxy_C_vec.add(i_tab1.get(b * 3 + 0));
                            carboxy_C_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        } else {
                            c_group = false;
                        }
                    }
                }
                if ((O_num == 2) && (c_group)) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (carboxy_C_vec.get(1) != 2) {
                        temp_vec.add(carboxy_C_vec.get(2));
                        temp_vec.add(carboxy_C_vec.get(3));
                        temp_vec.add(carboxy_C_vec.get(0));
                        temp_vec.add(carboxy_C_vec.get(1));
                    } else {
                        temp_vec = carboxy_C_vec;
                    }
                    C_index_A.add(a + 1);
                    carb_vec_A.add(temp_vec);
                }
            }
        }
        int C_index_A_size = C_index_A.size();

        //find carboxyl groups of molecule 2
        List<Integer> C_index_B = new ArrayList<>();
        List<List<Integer>> carb_vec_B = new ArrayList<>();
        for (int a = 0; a < atom_num_H_2; a++) {
            if ("C".equals(atomstr2.get(a).getSymbol())) {
                int O_num = 0;
                List<Integer> carboxy_C_vec = new ArrayList<>();
                boolean c_group = true;
                for (int b = 0; b < bond_number2; b++) {
                    if ((a + 1 == i_tab2.get(b * 3 + 0)) && ("O".equals(c_tab2.get(b * 2 + 1)))) {
                        if (check(i_tab2.get(b * 3 + 1), 2)) {
                            carboxy_C_vec.add(i_tab2.get(b * 3 + 1));
                            carboxy_C_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        } else {
                            c_group = false;
                        }
                    }
                    if ((a + 1 == i_tab2.get(b * 3 + 1)) && ("O".equals(c_tab2.get(b * 2 + 0)))) {
                        if (check(i_tab2.get(b * 3 + 0), 2)) {
                            carboxy_C_vec.add(i_tab2.get(b * 3 + 0));
                            carboxy_C_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        } else {
                            c_group = false;
                        }
                    }
                }
                if ((O_num == 2) && (c_group)) {
                    //Umsortieren, falls Doppelbindung nicht vorne
                    List<Integer> temp_vec = new ArrayList<>();
                    if (carboxy_C_vec.get(1) != 2) {
                        temp_vec.add(carboxy_C_vec.get(2));
                        temp_vec.add(carboxy_C_vec.get(3));
                        temp_vec.add(carboxy_C_vec.get(0));
                        temp_vec.add(carboxy_C_vec.get(1));
                    } else {
                        temp_vec = carboxy_C_vec;
                    }
                    C_index_B.add(a + 1);
                    carb_vec_B.add(temp_vec);
                }
            }
        }
        int C_index_B_size = C_index_B.size();

        List<List<Integer>> carboxy_final_MAPPINGS = new ArrayList<>();
        boolean carboxy_groups_in_both = true;
        if (C_index_A_size == 0 || C_index_B_size == 0) {
            carboxy_groups_in_both = false;
            carboxy_final_MAPPINGS.addAll(getFinalMappings());
        }
        if (carboxy_groups_in_both) {
            getFinalMappings().stream().forEach((final_solution) -> {
                boolean map_correct = true;
                int a = 0;
                while (a < best_MAPPING_size && map_correct) {
                    //gehe ein Mapping durch
                    int b = 0;
                    boolean not_found = true;
                    while (b < C_index_A_size && not_found) {
                        if (final_solution.get(a * 2 + 0).intValue() == C_index_A.get(b)) {
                            not_found = false;
                            List<Integer> vector_A = carb_vec_A.get(b);
                            int first_A = vector_A.get(0);
                            int secon_A = vector_A.get(2);
                            int first_B = 0;
                            int secon_B = 0;
                            boolean mapped_on_a_carboxy_group = false;
                            for (int c = 0; c < C_index_B_size; c++) {
                                if (final_solution.get(a * 2 + 1).intValue() == C_index_B.get(c)) {

                                    List<Integer> vector_B = carb_vec_B.get(c);
                                    first_B = vector_B.get(0);
                                    secon_B = vector_B.get(2);
                                    mapped_on_a_carboxy_group = true;
                                }
                            }

                            if (mapped_on_a_carboxy_group) {
                                boolean miss_map = true;
                                for (int c = 0; c < best_MAPPING_size; c++) {
                                    if ((final_solution.get(c * 2 + 0) == first_A) && (final_solution.get(c * 2 + 1) == first_B)) {
                                        miss_map = false;
                                    }
                                }
                                if (miss_map) {
                                    map_correct = false;
                                }
                            }
                        }
                        b++;
                    }
                    a++;
                }
                if (map_correct) {
                    carboxy_final_MAPPINGS.add(final_solution);
                }
            });
        }
        //2. Searching for phosphate groups  
        //find phosphate groups of molecule 1
        List<Integer> P_index_A = new ArrayList<>();
        List<List<Integer>> phos_vec_A = new ArrayList<>();
        for (int a = 0; a < atom_num_H_1; a++) {
            if ("P".equals(atomstr1.get(a).getSymbol())) {
                int P_num = 0;
                List<Integer> phos_P_vec = new ArrayList<>();
                for (int b = 0; b < bond_number1; b++) {
                    if ((a + 1 == i_tab1.get(b * 3 + 0)) && ("O".equals(c_tab1.get(b * 2 + 1)))) {
                        if (check(i_tab1.get(b * 3 + 1), 1)) {
                            phos_P_vec.add(i_tab1.get(b * 3 + 1));
                            phos_P_vec.add(i_tab1.get(b * 3 + 2));
                            P_num++;
                        }
                    }
                    if ((a + 1 == i_tab1.get(b * 3 + 1)) && ("O".equals(c_tab1.get(b * 2 + 0)))) {
                        if (check(i_tab1.get(b * 3 + 0), 1)) {
                            phos_P_vec.add(i_tab1.get(b * 3 + 0));
                            phos_P_vec.add(i_tab1.get(b * 3 + 2));
                            P_num++;
                        }
                    }
                }
                if (P_num == 1) {
                    List<Integer> temp_vec = new ArrayList<>();
                    temp_vec.add(phos_P_vec.get(0));
                    temp_vec.add(phos_P_vec.get(1));
                    P_index_A.add(a + 1);
                    phos_vec_A.add(temp_vec);
                }
                if (P_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) != 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                    } else {
                        temp_vec = phos_P_vec;
                    }
                    P_index_A.add(a + 1);
                    phos_vec_A.add(temp_vec);
                }
                if (P_num == 3) {
                    boolean no_double = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) == 2) {
                        temp_vec = phos_P_vec;
                        no_double = false;
                    }
                    if (phos_P_vec.get(3) == 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        no_double = false;
                    }
                    if (phos_P_vec.get(5) == 2) {
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        no_double = false;
                    }
                    if (no_double) {
                        temp_vec = phos_P_vec;
                    }
                    P_index_A.add(a + 1);
                    phos_vec_A.add(temp_vec);
                }
                if (P_num == 4) {
                    boolean no_double = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) == 2) {
                        temp_vec = phos_P_vec;
                        no_double = false;
                    }
                    if (phos_P_vec.get(3) == 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        no_double = false;
                    }
                    if (phos_P_vec.get(5) == 2) {
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        no_double = false;
                    }
                    if (phos_P_vec.get(7) == 2) {
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        no_double = false;
                    }
                    if (no_double) {
                        temp_vec = phos_P_vec;
                    }
                    P_index_A.add(a + 1);
                    phos_vec_A.add(temp_vec);
                }
            }
        }
        int P_index_A_size = P_index_A.size();
        //find phosphate groups of molecule 2
        List<Integer> P_index_B = new ArrayList<>();
        List<List<Integer>> phos_vec_B = new ArrayList<>();
        for (int a = 0; a < atom_num_H_2; a++) {
            if ("P".equals(atomstr2.get(a).getSymbol())) {
                int P_num = 0;
                List<Integer> phos_P_vec = new ArrayList<>();
                for (int b = 0; b < bond_number2; b++) {
                    if ((a + 1 == i_tab2.get(b * 3 + 0)) && ("O".equals(c_tab2.get(b * 2 + 1)))) {
                        if (check(i_tab2.get(b * 3 + 1), 2)) {
                            phos_P_vec.add(i_tab2.get(b * 3 + 1));
                            phos_P_vec.add(i_tab2.get(b * 3 + 2));
                            P_num++;
                        }
                    }
                    if ((a + 1 == i_tab2.get(b * 3 + 1)) && ("O".equals(c_tab2.get(b * 2 + 0)))) {
                        if (check(i_tab2.get(b * 3 + 0), 2)) {
                            phos_P_vec.add(i_tab2.get(b * 3 + 0));
                            phos_P_vec.add(i_tab2.get(b * 3 + 2));
                            P_num++;
                        }
                    }
                }
                if (P_num == 1) {
                    List<Integer> temp_vec = new ArrayList<>();
                    temp_vec.add(phos_P_vec.get(0));
                    temp_vec.add(phos_P_vec.get(1));
                    P_index_B.add(a + 1);
                    phos_vec_B.add(temp_vec);
                }
                if (P_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) != 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                    } else {
                        temp_vec = phos_P_vec;
                    }
                    P_index_B.add(a + 1);
                    phos_vec_B.add(temp_vec);
                }
                if (P_num == 3) {
                    boolean no_double = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) == 2) {
                        temp_vec = phos_P_vec;
                        no_double = false;
                    }
                    if (phos_P_vec.get(3) == 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        no_double = false;
                    }
                    if (phos_P_vec.get(5) == 2) {
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        no_double = false;
                    }
                    if (no_double) {
                        temp_vec = phos_P_vec;
                    }
                    P_index_B.add(a + 1);
                    phos_vec_B.add(temp_vec);
                }
                if (P_num == 4) {
                    boolean no_double = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (phos_P_vec.get(1) == 2) {
                        temp_vec = phos_P_vec;
                    }
                    if (phos_P_vec.get(3) == 2) {
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        no_double = false;
                    }
                    if (phos_P_vec.get(5) == 2) {
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        no_double = false;
                    }
                    if (phos_P_vec.get(7) == 2) {
                        temp_vec.add(phos_P_vec.get(6));
                        temp_vec.add(phos_P_vec.get(7));
                        temp_vec.add(phos_P_vec.get(0));
                        temp_vec.add(phos_P_vec.get(1));
                        temp_vec.add(phos_P_vec.get(2));
                        temp_vec.add(phos_P_vec.get(3));
                        temp_vec.add(phos_P_vec.get(4));
                        temp_vec.add(phos_P_vec.get(5));
                        no_double = false;
                    }
                    if (no_double) {
                        temp_vec = phos_P_vec;
                    }
                    P_index_B.add(a + 1);
                    phos_vec_B.add(temp_vec);
                }
            }
        }
        int P_index_B_size = P_index_B.size();

        List<List<Integer>> phosphate_final_MAPPINGS = new ArrayList<>();
        boolean phosphate_groups_in_both = true;
        if ((P_index_A_size == 0) || (P_index_B_size == 0)) {
            phosphate_groups_in_both = false;
            phosphate_final_MAPPINGS = carboxy_final_MAPPINGS;
        }
        if (phosphate_groups_in_both) {
            for (List<Integer> carb_solution : carboxy_final_MAPPINGS) {
                boolean map_correct = true;
                int a = 0;
                while ((a < best_MAPPING_size) && (map_correct)) {
                    int b = 0;
                    boolean not_found = true;
                    while ((b < P_index_A_size) && (not_found)) {
                        if (carb_solution.get(a * 2 + 0).equals(P_index_A.get(b))) {
                            not_found = false;
                            List<Integer> vector_A = phos_vec_A.get(b);
                            int first_A = vector_A.get(0);
                            //int secon_A = vector_A.get(2);
                            int first_B = 0;
                            //int secon_B = 0;
                            boolean mapped_on_a_phos_group = false;
                            for (int c = 0; c < P_index_B_size; c++) {
                                if (carb_solution.get(a * 2 + 1).equals(P_index_B.get(c))) {
                                    List<Integer> vector_B = phos_vec_B.get(c);
                                    first_B = vector_B.get(0);
                                    //secon_B = vector_B.get(2);
                                    mapped_on_a_phos_group = true;
                                }
                            }
                            if (mapped_on_a_phos_group) {
                                boolean miss_map = true;
                                for (int c = 0; c < best_MAPPING_size; c++) {
                                    if ((carb_solution.get(c * 2 + 0) == first_A) && (carb_solution.get(c * 2 + 1) == first_B)) {
                                        miss_map = false;
                                    }
                                }
                                if (miss_map) {
                                    map_correct = false;
                                }
                            }
                        }
                        b++;
                    }
                    a++;
                }
                if (map_correct) {
                    phosphate_final_MAPPINGS.add(carb_solution);
                }
            }
        }
        //3. Searching for Amino groups 2H N - C - N H2   
//find Amino-Carbon groups of molecule 1
        List<Integer> N_index_A = new ArrayList<>();
        List<List<Integer>> amino_vec_A = new ArrayList<>();
        for (int a = 0; a < atom_num_H_1; a++) {
            if (atomstr1.get(a).getSymbol().equals("C")) {
                int N_num = 0;
                List<Integer> amino_N_vec = new ArrayList<>();
                for (int b = 0; b < bond_number1; b++) {
                    if ((a + 1 == i_tab1.get(b * 3 + 0)) && ("N".equals(c_tab1.get(b * 2 + 1)))) {
                        if (check(i_tab1.get(b * 3 + 1), 1)) {
                            amino_N_vec.add(i_tab1.get(b * 3 + 1));
                            amino_N_vec.add(i_tab1.get(b * 3 + 2));
                            N_num++;
                        }
                    }
                    if ((a + 1 == i_tab1.get(b * 3 + 1)) && ("N".equals(c_tab1.get(b * 2 + 0)))) {
                        if (check(i_tab1.get(b * 3 + 0), 1)) {
                            amino_N_vec.add(i_tab1.get(b * 3 + 0));
                            amino_N_vec.add(i_tab1.get(b * 3 + 2));
                            N_num++;
                        }
                    }
                }
                if (N_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (amino_N_vec.get(1) != 2) {
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                    } else {
                        temp_vec = amino_N_vec;
                    }
                    N_index_A.add(a + 1);
                    amino_vec_A.add(temp_vec);
                }
                if (N_num == 3) {
                    boolean no_double_bond = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (amino_N_vec.get(1) == 2) {
                        temp_vec = amino_N_vec;
                        no_double_bond = false;
                    }
                    if (amino_N_vec.get(3) == 2) {
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                        temp_vec.add(amino_N_vec.get(4));
                        temp_vec.add(amino_N_vec.get(5));
                        no_double_bond = false;
                    }
                    if (amino_N_vec.get(5) == 2) {
                        temp_vec.add(amino_N_vec.get(4));
                        temp_vec.add(amino_N_vec.get(5));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        no_double_bond = false;
                    }
                    if (no_double_bond) {
                        temp_vec = amino_N_vec;
                    }
                    N_index_A.add(a + 1);
                    amino_vec_A.add(temp_vec);
                }
            }
        }
        int N_index_A_size = N_index_A.size();

        //find Amino-Carbon groups of molecule 2
        List<Integer> N_index_B = new ArrayList<>();
        List<List<Integer>> amino_vec_B = new ArrayList<>();
        for (int a = 0; a < atom_num_H_2; a++) {
            if ("C".equals(atomstr2.get(a).getSymbol())) {
                int N_num = 0;
                List<Integer> amino_N_vec = new ArrayList<>();
                for (int b = 0; b < bond_number2; b++) {
                    if ((a + 1 == i_tab2.get(b * 3 + 0)) && ("N".equals(c_tab2.get(b * 2 + 1)))) {
                        if (check(i_tab2.get(b * 3 + 1), 2)) {
                            amino_N_vec.add(i_tab2.get(b * 3 + 1));
                            amino_N_vec.add(i_tab2.get(b * 3 + 2));
                            N_num++;
                        }
                    }
                    if ((a + 1 == i_tab2.get(b * 3 + 1)) && ("N".equals(c_tab2.get(b * 2 + 0)))) {
                        if (check(i_tab2.get(b * 3 + 0), 2)) {
                            amino_N_vec.add(i_tab2.get(b * 3 + 0));
                            amino_N_vec.add(i_tab2.get(b * 3 + 2));
                            N_num++;
                        }
                    }
                }
                if (N_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (amino_N_vec.get(1) != 2) {
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                    } else {
                        temp_vec = amino_N_vec;
                    }
                    N_index_B.add(a + 1);
                    amino_vec_B.add(temp_vec);
                }
                if (N_num == 3) {
                    boolean no_double_bond = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (amino_N_vec.get(1) == 2) {
                        temp_vec = amino_N_vec;
                        no_double_bond = false;
                    }
                    if (amino_N_vec.get(3) == 2) {
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                        temp_vec.add(amino_N_vec.get(4));
                        temp_vec.add(amino_N_vec.get(5));
                        no_double_bond = false;
                    }
                    if (amino_N_vec.get(5) == 2) {
                        temp_vec.add(amino_N_vec.get(4));
                        temp_vec.add(amino_N_vec.get(5));
                        temp_vec.add(amino_N_vec.get(0));
                        temp_vec.add(amino_N_vec.get(1));
                        temp_vec.add(amino_N_vec.get(2));
                        temp_vec.add(amino_N_vec.get(3));
                        no_double_bond = false;
                    }
                    if (no_double_bond) {
                        temp_vec = amino_N_vec;
                    }
                    N_index_B.add(a + 1);
                    amino_vec_B.add(temp_vec);
                }
            }
        }
        int N_index_B_size = N_index_B.size();

        List<List<Integer>> amino_final_MAPPINGS = new ArrayList<>();
        boolean amino_groups_in_both = true;
        if ((N_index_A_size == 0) || (N_index_B_size == 0)) {
            amino_groups_in_both = false;
            amino_final_MAPPINGS = phosphate_final_MAPPINGS;
        }
        if (amino_groups_in_both) {
            for (List<Integer> phosphate_solution : phosphate_final_MAPPINGS) {
                boolean map_correct = true;
                int a = 0;
                while ((a < best_MAPPING_size) && (map_correct)) {
                    int b = 0;
                    boolean not_found = true;
                    while ((b < N_index_A_size) && (not_found)) {
                        if (phosphate_solution.get(a * 2 + 0).equals(N_index_A.get(b))) {
                            not_found = false;
                            List<Integer> vector_A = amino_vec_A.get(b);
                            int first_A = vector_A.get(0);
                            int secon_A = vector_A.get(2);
                            int first_B = 0;
                            int secon_B = 0;
                            boolean mapped_on_a_amino_group = false;
                            for (int c = 0; c < N_index_B_size; c++) {
                                if (phosphate_solution.get(a * 2 + 1).equals(N_index_B.get(c))) {
                                    List<Integer> vector_B = amino_vec_B.get(c);
                                    first_B = vector_B.get(0);
                                    secon_B = vector_B.get(2);
                                    mapped_on_a_amino_group = true;
                                }
                            }
                            if (mapped_on_a_amino_group) {
                                boolean miss_map = true;
                                for (int c = 0; c < best_MAPPING_size; c++) {
                                    if ((phosphate_solution.get(c * 2 + 0) == first_A) && (phosphate_solution.get(c * 2 + 1) == first_B)) {
                                        miss_map = false;
                                    }
                                }
                                if (miss_map) {
                                    map_correct = false;
                                }
                            }
                        }
                        b++;
                    }
                    a++;
                }
                if (map_correct) {
                    amino_final_MAPPINGS.add(phosphate_solution);
                }
            }
        }

        //4. Searching for Sulfo groups S-O3  
//find sulfo groups of molecule 1
        List<Integer> SO_index_A = new ArrayList<>();
        List<List<Integer>> sulfo_vec_A = new ArrayList<>();
        for (int a = 0; a < atom_num_H_1; a++) {
            if ("S".equals(atomstr1.get(a).getSymbol())) {
                int O_num = 0;
                List<Integer> sulfo_S_vec = new ArrayList<>();
                for (int b = 0; b < bond_number1; b++) {
                    if ((a + 1 == i_tab1.get(b * 3 + 0)) && ("O".equals(c_tab1.get(b * 2 + 1)))) {
                        if (check(i_tab1.get(b * 3 + 1), 1)) {
                            sulfo_S_vec.add(i_tab1.get(b * 3 + 1));
                            sulfo_S_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                    if ((a + 1 == i_tab1.get(b * 3 + 1)) && ("O".equals(c_tab1.get(b * 2 + 0)))) {
                        if (check(i_tab1.get(b * 3 + 0), 1)) {
                            sulfo_S_vec.add(i_tab1.get(b * 3 + 0));
                            sulfo_S_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                }
                if (O_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (sulfo_S_vec.get(1) != 2) {
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                    } else {
                        temp_vec = sulfo_S_vec;
                    }
                    SO_index_A.add(a + 1);
                    sulfo_vec_A.add(temp_vec);
                }
                if (O_num == 3) {
                    boolean no_single_bond = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (sulfo_S_vec.get(5) == 1) {
                        temp_vec = sulfo_S_vec;
                        no_single_bond = false;
                    }
                    if (sulfo_S_vec.get(1) == 1) {
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        temp_vec.add(sulfo_S_vec.get(4));
                        temp_vec.add(sulfo_S_vec.get(5));
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                        no_single_bond = false;
                    }
                    if (sulfo_S_vec.get(3) == 1) {
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                        temp_vec.add(sulfo_S_vec.get(4));
                        temp_vec.add(sulfo_S_vec.get(5));
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        no_single_bond = false;
                    }
                    if (no_single_bond) {
                        temp_vec = sulfo_S_vec;
                    }
                    SO_index_A.add(a + 1);
                    sulfo_vec_A.add(temp_vec);
                }
                if (O_num == 4) {
                    List<Integer> temp_vec = new ArrayList<>();
                    List<Integer> single_bond_posi = new ArrayList<>();
                    List<Integer> double_bond_posi = new ArrayList<>();
                    for (int c = 1; c < 8; c = c + 2) {
                        if (sulfo_S_vec.get(c) == 1) {
                            single_bond_posi.add(sulfo_S_vec.get(c - 1));
                            single_bond_posi.add(sulfo_S_vec.get(c));
                        }
                        if (sulfo_S_vec.get(c) == 2) {
                            double_bond_posi.add(sulfo_S_vec.get(c - 1));
                            double_bond_posi.add(sulfo_S_vec.get(c));
                        }
                    }
                    int d_b_posi_size = double_bond_posi.size();
                    for (int c = 0; c < d_b_posi_size; c = c + 2) {
                        temp_vec.add(double_bond_posi.get(c));
                        temp_vec.add(double_bond_posi.get(c + 1));
                    }
                    int s_b_posi_size = single_bond_posi.size();
                    for (int c = 0; c < s_b_posi_size; c = c + 2) {
                        temp_vec.add(single_bond_posi.get(c));
                        temp_vec.add(single_bond_posi.get(c + 1));
                    }
                    SO_index_A.add(a + 1);
                    sulfo_vec_A.add(temp_vec);
                }
            }
        }
        int SO_index_A_size = SO_index_A.size();

        //find sulfo groups of molecule 2
        List<Integer> SO_index_B = new ArrayList<>();
        List<List<Integer>> sulfo_vec_B = new ArrayList<>();
        for (int a = 0; a < atom_num_H_2; a++) {
            if (atomstr2.get(a).getSymbol().equals("S")) {
                int O_num = 0;
                List<Integer> sulfo_S_vec = new ArrayList<>();
                for (int b = 0; b < bond_number2; b++) {
                    if ((a + 1 == i_tab2.get(b * 3 + 0)) && ("O".equals(c_tab2.get(b * 2 + 1)))) {
                        if (check(i_tab2.get(b * 3 + 1), 2)) {
                            sulfo_S_vec.add(i_tab2.get(b * 3 + 1));
                            sulfo_S_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                    if ((a + 1 == i_tab2.get(b * 3 + 1)) && ("O".equals(c_tab2.get(b * 2 + 0)))) {
                        if (check(i_tab2.get(b * 3 + 0), 2)) {
                            sulfo_S_vec.add(i_tab2.get(b * 3 + 0));
                            sulfo_S_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                }
                if (O_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (sulfo_S_vec.get(1) != 2) {
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                    } else {
                        temp_vec = sulfo_S_vec;
                    }
                    SO_index_B.add(a + 1);
                    sulfo_vec_B.add(temp_vec);
                }
                if (O_num == 3) {
                    boolean no_single_bond = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    if (sulfo_S_vec.get(5) == 1) {
                        temp_vec = sulfo_S_vec;
                        no_single_bond = false;
                    }
                    if (sulfo_S_vec.get(1) == 1) {
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        temp_vec.add(sulfo_S_vec.get(4));
                        temp_vec.add(sulfo_S_vec.get(5));
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                        no_single_bond = false;
                    }
                    if (sulfo_S_vec.get(3) == 1) {
                        temp_vec.add(sulfo_S_vec.get(0));
                        temp_vec.add(sulfo_S_vec.get(1));
                        temp_vec.add(sulfo_S_vec.get(4));
                        temp_vec.add(sulfo_S_vec.get(5));
                        temp_vec.add(sulfo_S_vec.get(2));
                        temp_vec.add(sulfo_S_vec.get(3));
                        no_single_bond = false;
                    }
                    if (no_single_bond) {
                        temp_vec = sulfo_S_vec;
                    }
                    SO_index_B.add(a + 1);
                    sulfo_vec_B.add(temp_vec);
                }
                if (O_num == 4) {
                    List<Integer> temp_vec = new ArrayList<>();
                    List<Integer> single_bond_posi = new ArrayList<>();
                    List<Integer> double_bond_posi = new ArrayList<>();
                    for (int c = 1; c < 8; c = c + 2) {
                        if (sulfo_S_vec.get(c) == 1) {
                            single_bond_posi.add(sulfo_S_vec.get(c - 1));
                            single_bond_posi.add(sulfo_S_vec.get(c));
                        }
                        if (sulfo_S_vec.get(c) == 2) {
                            double_bond_posi.add(sulfo_S_vec.get(c - 1));
                            double_bond_posi.add(sulfo_S_vec.get(c));
                        }
                    }
                    int d_b_posi_size = double_bond_posi.size();
                    for (int c = 0; c < d_b_posi_size; c = c + 2) {
                        temp_vec.add(double_bond_posi.get(c));
                        temp_vec.add(double_bond_posi.get(c + 1));
                    }
                    int s_b_posi_size = single_bond_posi.size();
                    for (int c = 0; c < s_b_posi_size; c = c + 2) {
                        temp_vec.add(single_bond_posi.get(c));
                        temp_vec.add(single_bond_posi.get(c + 1));
                    }
                    SO_index_B.add(a + 1);
                    sulfo_vec_B.add(temp_vec);
                }
            }
        }
        int SO_index_B_size = SO_index_B.size();

        List<List<Integer>> sulfo_final_MAPPINGS = new ArrayList<>();
        boolean sulfo_groups_in_both = true;
        if ((SO_index_A_size == 0) || (SO_index_B_size == 0)) {
            sulfo_groups_in_both = false;
            sulfo_final_MAPPINGS = amino_final_MAPPINGS;
        }
        if (sulfo_groups_in_both) {
            boolean no_correct_mapping = true;
            for (List<Integer> amino_solution : amino_final_MAPPINGS) {
                boolean map_correct = true;
                int a = 0;
                while ((a < best_MAPPING_size) && (map_correct)) {
                    int b = 0;
                    boolean not_found = true;
                    while ((b < SO_index_A_size) && (not_found)) {
                        if (amino_solution.get(a * 2 + 0).equals(SO_index_A.get(b))) {
                            not_found = false;
                            List<Integer> vector_A = sulfo_vec_A.get(b);
                            int first_A = vector_A.get(0);
                            int secon_A = vector_A.get(2);
                            int first_B = 0;
                            int secon_B = 0;
                            boolean mapped_on_a_sulfo_group = false;
                            for (int c = 0; c < SO_index_B_size; c++) {
                                if (amino_solution.get(a * 2 + 1).equals(SO_index_B.get(c))) {
                                    List<Integer> vector_B = sulfo_vec_B.get(c);
                                    first_B = vector_B.get(0);
                                    secon_B = vector_B.get(2);
                                    mapped_on_a_sulfo_group = true;
                                }
                            }
                            if (mapped_on_a_sulfo_group) {
                                boolean miss_map1 = true;
                                boolean miss_map2 = true;
                                for (int c = 0; c < best_MAPPING_size; c++) {
                                    if ((amino_solution.get(c * 2 + 0) == first_A) && (amino_solution.get(c * 2 + 1) == first_B)) {
                                        miss_map1 = false;
                                    }
                                    if ((amino_solution.get(c * 2 + 0) == secon_A) && (amino_solution.get(c * 2 + 1) == secon_B)) {
                                        miss_map2 = false;
                                    }
                                }
                                if ((miss_map1) || (miss_map2)) {
                                    map_correct = false;
                                }
                                if (miss_map1) {
                                    map_correct = false;
                                }
                            }
                        }
                        b++;
                    }
                    a++;
                }
                if (map_correct) {
                    sulfo_final_MAPPINGS.add(amino_solution);
                    no_correct_mapping = false;
                }
            }
            if (no_correct_mapping) {
                List<Integer> Mol_A_Os = new ArrayList<>();
                for (int a = 0; a < SO_index_A_size; a++) {
                    List<Integer> s_v_A = sulfo_vec_A.get(a);
                    int s_v_A_size = s_v_A.size();
                    for (int b = 0; b < s_v_A_size; b = b + 2) {
                        Mol_A_Os.add(s_v_A.get(b));
                    }
                }
                int Mol_A_Os_size = Mol_A_Os.size();
                List<Integer> Mol_B_Os = new ArrayList<>();
                for (int a = 0; a < SO_index_B_size; a++) {
                    List<Integer> s_v_B = sulfo_vec_B.get(a);
                    int s_v_B_size = s_v_B.size();
                    for (int b = 0; b < s_v_B_size; b = b + 2) {
                        Mol_B_Os.add(s_v_B.get(b));
                    }
                }
                int Mol_B_Os_size = Mol_B_Os.size();
                List<List<Integer>> temp_s_f_M = new ArrayList<>();
                for (List<Integer> amino_solution : phosphate_final_MAPPINGS) {
                    List<Integer> t_map = new ArrayList<>();
                    for (int a = 0; a < best_MAPPING_size; a = a + 2) {
                        boolean store = true;
                        int b = 0;
                        while ((b < Mol_A_Os_size) && (store)) {
                            int c = 0;
                            while ((c < Mol_B_Os_size) && (store)) {
                                if ((Mol_A_Os.get(b).equals(amino_solution.get(a))) && (Mol_B_Os.get(c).equals(amino_solution.get(a + 1)))) {
                                    store = false;
                                }
                                c++;
                            }
                            b++;
                        }
                        if (store) {
                            t_map.add(amino_solution.get(a));
                            t_map.add(amino_solution.get(a + 1));
                        }
                    }
                    temp_s_f_M.add(t_map);
                }
                int temp_s_f_M_size = temp_s_f_M.size();

                List<List<Integer>> temp_s_f_M2 = new ArrayList<>();
                for (int a = 0; a < temp_s_f_M_size; a++) {
                    List<Integer> map_A = temp_s_f_M.get(a);
                    int map_A_size = map_A.size();
                    boolean unique_map = true;
                    int b = a + 1;
                    while ((b < temp_s_f_M_size) && (unique_map)) {
                        if (a == temp_s_f_M_size - 1) {
                            break;
                        }
                        boolean map_contained = true;
                        List<Integer> map_B = temp_s_f_M.get(b);
                        int map_B_size = map_B.size();
                        int c = 0;
                        while ((c < map_A_size) && (map_contained)) {
                            boolean map_not_contained = true;
                            int d = 0;
                            while ((d < map_B_size) && (map_not_contained)) {
                                if ((map_A.get(c).equals(map_B.get(d))) && (map_A.get(c + 1).equals(map_B.get(d + 1)))) {
                                    map_not_contained = false;
                                }
                                d = d + 2;
                            }
                            if (map_not_contained) {
                                map_contained = false;
                            }
                            c = c + 2;
                        }
                        if (map_contained) {
                            unique_map = false;
                        }
                        b++;
                    }
                    if (unique_map) {
                        temp_s_f_M2.add(map_A);
                    }
                }
                int temp_s_f_M2_size = temp_s_f_M2.size();
                for (int a = 0; a < temp_s_f_M2_size; a++) {
                    List<Integer> map_element = temp_s_f_M2.get(a);
                    List<Integer> new_element = temp_s_f_M2.get(a);
                    int map_element_size = map_element.size();
                    for (int b = 0; b < map_element_size; b = a + 2) {
                        for (int c = 0; c < SO_index_A_size; c++) {
                            if (map_element.get(b).equals(SO_index_A.get(c))) {
                                List<Integer> Os_A = sulfo_vec_A.get(c);
                                List<Integer> Os_B = new ArrayList<>();
                                for (int d = 0; d < SO_index_B_size; d++) {
                                    if (map_element.get(b + 1).equals(SO_index_B.get(d))) {
                                        Os_B = sulfo_vec_B.get(d);
                                    }
                                }
                                int Os_A_size = Os_A.size();
                                int Os_B_size = Os_B.size();
                                int e = 0;
                                while ((e < Os_A_size) && (e < Os_B_size)) {
                                    new_element.add(Os_A.get(e));
                                    new_element.add(Os_B.get(e));
                                    e = e + 2;
                                }
                            }
                        }
                    }
                    sulfo_final_MAPPINGS.add(new_element);
                }
            }
        }
        //5. Searching for Nitro groups N-O3   
        //find nitro groups of molecule 1
        List<Integer> NO_index_A = new ArrayList<>();
        List<List<Integer>> nitro_vec_A = new ArrayList<>();
        for (int a = 0; a < atom_num_H_1; a++) {
            if (atomstr1.get(a).getSymbol().equals("N")) {
                int O_num = 0;
                List<Integer> nitro_N_vec = new ArrayList<>();
                for (int b = 0; b < bond_number1; b++) {
                    if ((a + 1 == i_tab1.get(b * 3 + 0)) && ("O".equals(c_tab1.get(b * 2 + 1)))) {
                        if (check(i_tab1.get(b * 3 + 1), 1)) {
                            nitro_N_vec.add(i_tab1.get(b * 3 + 1));
                            nitro_N_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                    if ((a + 1 == i_tab1.get(b * 3 + 1)) && ("O".equals(c_tab1.get(b * 2 + 0)))) {
                        if (check(i_tab1.get(b * 3 + 0), 1)) {
                            nitro_N_vec.add(i_tab1.get(b * 3 + 0));
                            nitro_N_vec.add(i_tab1.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                }
                if (O_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (nitro_N_vec.get(1) != 2) {
                        temp_vec.add(nitro_N_vec.get(2));
                        temp_vec.add(nitro_N_vec.get(3));
                        temp_vec.add(nitro_N_vec.get(0));
                        temp_vec.add(nitro_N_vec.get(1));
                    } else {
                        temp_vec = nitro_N_vec;
                    }
                    NO_index_A.add(a + 1);
                    nitro_vec_A.add(temp_vec);
                }
                if (O_num == 3) {
                    boolean no_single_bond = true; //fr Fall, dass es keine Einfachbindung gibt
                    List<Integer> temp_vec = new ArrayList<>();
                    for (int c = 0; c < 6; c = c + 2) {
                        if (nitro_N_vec.get(c + 1) == 2) {
                            temp_vec.add(nitro_N_vec.get(c));
                            temp_vec.add(nitro_N_vec.get(c + 1));
                        }
                    }
                    for (int c = 0; c < 6; c = c + 2) {
                        if (nitro_N_vec.get(c + 1) == 1) {
                            temp_vec.add(nitro_N_vec.get(c));
                            temp_vec.add(nitro_N_vec.get(c + 1));
                            no_single_bond = false;
                        }
                    }
                    if (no_single_bond) {
                        temp_vec = nitro_N_vec;
                    }
                    NO_index_A.add(a + 1);
                    nitro_vec_A.add(temp_vec);
                }
            }
        }
        int NO_index_A_size = NO_index_A.size();

        //find nitro groups of molecule 2
        List<Integer> NO_index_B = new ArrayList<>();
        List<List<Integer>> nitro_vec_B = new ArrayList<>();
        for (int a = 0; a < atom_num_H_2; a++) {
            if ("N".equals(atomstr2.get(a).getSymbol())) {
                int O_num = 0;
                List<Integer> nitro_N_vec = new ArrayList<>();
                for (int b = 0; b < bond_number2; b++) {
                    if ((a + 1 == i_tab2.get(b * 3 + 0)) && ("O".equals(c_tab2.get(b * 2 + 1)))) {
                        if (check(i_tab2.get(b * 3 + 1), 2)) {
                            nitro_N_vec.add(i_tab2.get(b * 3 + 1));
                            nitro_N_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                    if ((a + 1 == i_tab2.get(b * 3 + 1)) && ("O".equals(c_tab2.get(b * 2 + 0)))) {
                        if (check(i_tab2.get(b * 3 + 0), 2)) {
                            nitro_N_vec.add(i_tab2.get(b * 3 + 0));
                            nitro_N_vec.add(i_tab2.get(b * 3 + 2));
                            O_num++;
                        }
                    }
                }
                if (O_num == 2) {
                    List<Integer> temp_vec = new ArrayList<>();
                    if (nitro_N_vec.get(1) != 2) {
                        temp_vec.add(nitro_N_vec.get(2));
                        temp_vec.add(nitro_N_vec.get(3));
                        temp_vec.add(nitro_N_vec.get(0));
                        temp_vec.add(nitro_N_vec.get(1));
                    } else {
                        temp_vec = nitro_N_vec;
                    }
                    NO_index_B.add(a + 1);
                    nitro_vec_B.add(temp_vec);
                }
                if (O_num == 3) {
                    boolean no_single_bond = true;
                    List<Integer> temp_vec = new ArrayList<>();
                    for (int c = 0; c < 6; c = c + 2) {
                        if (nitro_N_vec.get(c + 1) == 2) {
                            temp_vec.add(nitro_N_vec.get(c));
                            temp_vec.add(nitro_N_vec.get(c + 1));
                        }
                    }
                    for (int c = 0; c < 6; c = c + 2) {
                        if (nitro_N_vec.get(c + 1) == 1) {
                            temp_vec.add(nitro_N_vec.get(c));
                            temp_vec.add(nitro_N_vec.get(c + 1));
                            no_single_bond = false;
                        }
                    }
                    if (no_single_bond) {
                        temp_vec = nitro_N_vec;
                    }
                    NO_index_B.add(a + 1);
                    nitro_vec_B.add(temp_vec);
                }
            }
        }
        int NO_index_B_size = NO_index_B.size();

        List<List<Integer>> nitro_final_MAPPINGS = new ArrayList<>();
        boolean nitro_groups_in_both = true;
        if ((NO_index_A_size == 0) || (NO_index_B_size == 0)) {
            nitro_groups_in_both = false;
            nitro_final_MAPPINGS = sulfo_final_MAPPINGS;
        }
        if (nitro_groups_in_both) {
            boolean no_correct_mapping = true;
            java.util.Iterator<List<Integer>> sulfo_MAP_iter = sulfo_final_MAPPINGS.iterator();
            while (sulfo_MAP_iter.hasNext()) {
                List<Integer> sulfo_solution = sulfo_MAP_iter.next();
                boolean map_correct = true;
                int a = 0;
                while ((a < best_MAPPING_size) && (map_correct)) {
                    int b = 0;
                    boolean not_found = true;
                    while ((b < NO_index_A_size) && (not_found)) {
                        if (sulfo_solution.get(a * 2 + 0).equals(NO_index_A.get(b))) {
                            not_found = false;
                            List<Integer> vector_A = nitro_vec_A.get(b);
                            int first_A = vector_A.get(0);
                            int secon_A = vector_A.get(2);
                            int first_B = 0;
                            int secon_B = 0;
                            boolean mapped_on_a_nitro_group = false;
                            for (int c = 0; c < NO_index_B_size; c++) {
                                if (sulfo_solution.get(a * 2 + 1).equals(NO_index_B.get(c))) {
                                    List<Integer> vector_B = nitro_vec_B.get(c);
                                    first_B = vector_B.get(0);
                                    secon_B = vector_B.get(2);
                                    mapped_on_a_nitro_group = true;
                                }
                            }
                            if (mapped_on_a_nitro_group) {
                                boolean miss_map1 = true;
                                boolean miss_map2 = true;
                                for (int c = 0; c < best_MAPPING_size; c++) {
                                    if ((sulfo_solution.get(c * 2 + 0) == first_A) && (sulfo_solution.get(c * 2 + 1) == first_B)) {
                                        miss_map1 = false;
                                    }
                                    if ((sulfo_solution.get(c * 2 + 0) == secon_A) && (sulfo_solution.get(c * 2 + 1) == secon_B)) {
                                        miss_map2 = false;
                                    }
                                }
                                if ((miss_map1) || (miss_map2)) {
                                    map_correct = false;
                                }
                                if (miss_map1) {
                                    map_correct = false;
                                }
                            }
                        }
                        b++;
                    }
                    a++;
                }
                if (map_correct) {
                    nitro_final_MAPPINGS.add(sulfo_solution);
                    no_correct_mapping = false;
                }
            }
            if (no_correct_mapping) {
                List<Integer> Mol_A_Os = new ArrayList<>();
                for (int a = 0; a < NO_index_A_size; a++) {
                    List<Integer> n_v_A = nitro_vec_A.get(a);
                    int n_v_A_size = n_v_A.size();
                    for (int b = 0; b < n_v_A_size; b = b + 2) {
                        Mol_A_Os.add(n_v_A.get(b));
                    }
                }
                int Mol_A_Os_size = Mol_A_Os.size();
                List<Integer> Mol_B_Os = new ArrayList<>();
                for (int a = 0; a < NO_index_B_size; a++) {
                    List<Integer> n_v_B = nitro_vec_B.get(a);
                    int n_v_B_size = n_v_B.size();
                    for (int b = 0; b < n_v_B_size; b = b + 2) {
                        Mol_B_Os.add(n_v_B.get(b));
                    }
                }
                int Mol_B_Os_size = Mol_B_Os.size();
                List<List<Integer>> temp_n_f_M = new ArrayList<>();
                java.util.Iterator<List<Integer>> sulfo_iter = sulfo_final_MAPPINGS.iterator();
                while (sulfo_iter.hasNext()) {
                    List<Integer> sulfo_solution = sulfo_iter.next();
                    List<Integer> t_map = new ArrayList<>();
                    for (int a = 0; a < best_MAPPING_size; a = a + 2) {
                        boolean store = true;
                        int b = 0;
                        while ((b < Mol_A_Os_size) && (store)) {
                            int c = 0;
                            while ((c < Mol_B_Os_size) && (store)) {
                                if ((Mol_A_Os.get(b).equals(sulfo_solution.get(a))) && (Mol_B_Os.get(c).equals(sulfo_solution.get(a + 1)))) {
                                    store = false;
                                }
                                c++;
                            }
                            b++;
                        }
                        if (store) {
                            t_map.add(sulfo_solution.get(a));
                            t_map.add(sulfo_solution.get(a + 1));
                        }
                    }
                    temp_n_f_M.add(t_map);
                }
                int temp_n_f_M_size = temp_n_f_M.size();

                List<List<Integer>> temp_n_f_M2 = new ArrayList<>();
                for (int a = 0; a < temp_n_f_M_size; a++) {
                    List<Integer> map_A = temp_n_f_M.get(a);
                    int map_A_size = map_A.size();
                    boolean unique_map = true;
                    int b = a + 1;
                    while ((b < temp_n_f_M_size) && (unique_map)) {
                        if (a == temp_n_f_M_size - 1) {
                            break;
                        }
                        boolean map_contained = true;
                        List<Integer> map_B = temp_n_f_M.get(b);
                        int map_B_size = map_B.size();
                        int c = 0;
                        while ((c < map_A_size) && (map_contained)) {
                            boolean map_not_contained = true;
                            int d = 0;
                            while ((d < map_B_size) && (map_not_contained)) {
                                if ((map_A.get(c).equals(map_B.get(d))) && (map_A.get(c + 1).equals(map_B.get(d + 1)))) {
                                    map_not_contained = false;
                                }
                                d = d + 2;
                            }
                            if (map_not_contained) {
                                map_contained = false;
                            }
                            c = c + 2;
                        }
                        if (map_contained) {
                            unique_map = false;
                        }
                        b++;
                    }
                    if (unique_map) {
                        temp_n_f_M2.add(map_A);
                    }
                }
                boolean no_correct_solution = true;
                int temp_n_f_M2_size = temp_n_f_M2.size();
                for (int a = 0; a < temp_n_f_M2_size; a++) {
                    //gehe durch temp_n_f_M2
                    List<Integer> map_element = temp_n_f_M2.get(a);
                    List<Integer> new_element = temp_n_f_M2.get(a);
                    int map_element_size = map_element.size();
                    for (int b = 0; b < map_element_size; b = a + 2) {
                        for (int c = 0; c < NO_index_A_size; c++) {
                            if (map_element.get(b).equals(NO_index_A.get(c))) {
                                List<Integer> Os_A = nitro_vec_A.get(c);
                                List<Integer> Os_B = new ArrayList<>();
                                for (int d = 0; d < NO_index_B_size; d++) {
                                    if (map_element.get(b + 1).equals(NO_index_B.get(d))) {
                                        Os_B = nitro_vec_B.get(d);
                                    }
                                }
                                int Os_A_size = Os_A.size();
                                int Os_B_size = Os_B.size();
                                int e = 0;
                                while ((e < Os_A_size) && (e < Os_B_size)) {
                                    new_element.add(Os_A.get(e));
                                    new_element.add(Os_B.get(e));
                                    e = e + 2;
                                }
                            }
                        }
                    }
                    if (new_element.size() == best_MAPPING_size) {
                        no_correct_solution = false;
                    } else {
                        nitro_final_MAPPINGS.add(new_element);
                    }
                }
                if (no_correct_solution) {
                    nitro_final_MAPPINGS.clear();
                    nitro_final_MAPPINGS = sulfo_final_MAPPINGS;
                }
            }
        }

//6. Searching for redundant Methyl-group mappings   
        getFinalMappings().clear();
        getFinalMappings().addAll(nitro_final_MAPPINGS);

        return 0;
    }
    //Third part: postfilter system

    boolean check(int atom, int molecule) {

        int count_neighb = 0;

        if (molecule == 1) {
            for (int a = 0; a < bond_number1; a++) {
                if ((atom == i_tab1.get(a * 3 + 0)) || (atom == i_tab1.get(a * 3 + 1))) {
                    count_neighb++;
                }
            }
        }

        if (molecule == 2) {
            for (int a = 0; a < bond_number2; a++) {
                if ((atom == i_tab2.get(a * 3 + 0)) || (atom == i_tab2.get(a * 3 + 1))) {
                    count_neighb++;
                }
            }
        }

        return count_neighb == 1;
    }

}
