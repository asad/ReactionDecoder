/**
 *
 * Copyright (C) 2009-2020 Syed Asad Rahman <asad at ebi.ac.uk>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.graph.Graph;
import org.openscience.smsd.graph.Vertex;

/**
 * This class handles MCS between two identical molecules. Hence they generate
 * am MCS where all atoms are mapped.
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MappingHandler {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MappingHandler.class);

    private static final boolean DEBUG = false;

    /**
     *
     * Extract atom getMapping from the cliques and stores it in a List
     *
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    private static List<Integer> extractCliqueMapping(List<Integer> comp_graph_nodes, List<Integer> clique_List_org) {

        List<Integer> clique_mapping = Collections.synchronizedList(new ArrayList<>());
        List<Integer> clique_List = new ArrayList<>(clique_List_org);
        int clique_siz = clique_List.size();
        int vec_size = comp_graph_nodes.size();
//        System.out.println("VEC  SIZE " + vec_size);
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b += 3) {
                if (Objects.equals(clique_List.get(a), comp_graph_nodes.get(b + 2))) {
                    clique_mapping.add(comp_graph_nodes.get(b));
                    clique_mapping.add(comp_graph_nodes.get(b + 1));
                }
            }
        }

        return clique_mapping;
    }

    //extract atom getMapping from the clique List and print it on the screen
    /**
     *
     * @param _mappings
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    public static List<List<Integer>> extractMapping(List<List<Integer>> _mappings, List<Integer> comp_graph_nodes,
            List<Integer> clique_List_org) {
        try {
            List<Integer> clique_List = extractCliqueMapping(comp_graph_nodes, clique_List_org);
            _mappings.add(clique_List);
        } catch (Exception e) {
            LOGGER.debug("Error in FinalMapping List: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return _mappings;
    }

    //extract atom getMapping from the clique List and print it on the screen
    /**
     *
     * @param comp_graph_nodes
     * @param clique_List_org
     * @return
     */
    public static Map<Integer, Integer> getMapping(List<Integer> comp_graph_nodes,
            Collection<Integer> clique_List_org) {
        Map<Integer, Integer> clique_mapping = Collections.synchronizedSortedMap(new TreeMap<>());

        try {
            List<Integer> clique_List = new ArrayList<>(clique_List_org);

//        System.out.println("VEC  SIZE " + vec_size);
            for (int a = 0; a < clique_List.size(); a++) {
                for (int b = 0; b < comp_graph_nodes.size(); b += 3) {
                    if (Objects.equals(clique_List.get(a), comp_graph_nodes.get(b + 2))) {
                        clique_mapping.put(comp_graph_nodes.get(b), comp_graph_nodes.get(b + 1));
                    }
                }
            }

//            System.out.println("atomatommapping  SIZE " + atomatommapping.size());
        } catch (Exception e) {
            LOGGER.debug("Error in FinalMapping List: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return clique_mapping;
    }

    //extract atom getMapping from the clique List and print it on the screen
    /**
     *
     * @param comp_graph_nodes
     * @param s
     * @param t
     * @param cliques
     * @return
     */
    public static Map<Integer, Integer> getMapping(
            Graph comp_graph_nodes,
            IAtomContainer s,
            IAtomContainer t,
            Set<Vertex> cliques,
            AtomMatcher am,
            BondMatcher bm) {
        TreeMap<Integer, Integer> bondCliques = new TreeMap<>();

        /*
         * Retrive Bond index for mapped vertices in the compatibility graph
         */
        cliques.forEach((v) -> {
            bondCliques.put(v.getQueryBondIndex(), v.getTargetBondIndex());
        });

        if (DEBUG) {
            System.out.println("Bond clique_mapping " + bondCliques);
        }

        Map<Integer, Integer> atomatommapping = makeAtomsMapOfBondsMap(bondCliques, s, t, am, bm);
        if (DEBUG) {
            System.out.println("bondCliques " + bondCliques.size());
            System.out.println("clique_mapping " + atomatommapping);
        }
        if (DEBUG) {
            try {
                System.out.println("mcs " + new SmilesGenerator(SmiFlavor.Generic)
                        .create(toSubstructures(bondCliques.keySet(), s)));
            } catch (CDKException ex) {
                LOGGER.error(Level.SEVERE, "Unable to extract mcs ", ex.getMessage());
            }
        }
        return atomatommapping;
    }

    /**
     * Returns matched sub graph
     *
     * @param bondMap
     * @param ac
     * @return
     */
    public static IAtomContainer toSubstructures(
            Set<Integer> bondMap,
            IAtomContainer ac) {

        final IAtomContainer submol = ac.getBuilder()
                .newInstance(IAtomContainer.class);
        bondMap.stream().map((b) -> {
            submol.addAtom(ac.getBond(b).getAtom(0));
            return b;
        }).map((b) -> {
            submol.addAtom(ac.getBond(b).getAtom(1));
            return b;
        }).forEachOrdered((b) -> {
            submol.addBond(ac.getBond(b));
        });

        return submol;
    }

    /**
     *
     * Creates a new instance of Post Filter and removes redundant mapping(s).
     *
     * @param mappings
     * @return Filtered non-redundant mappings
     */
    public synchronized static Set<Map<Integer, Integer>> filter(List<List<Integer>> mappings) {
        Set<Map<Integer, Integer>> final_MAPPINGS = new TreeSet<>();

        mappings.stream().map((map) -> {
            Map<Integer, Integer> mapping = new TreeMap<>();
            for (int i = 0; i < map.size(); i = i + 2) {
                mapping.put(map.get(i), map.get(i + 1));
            }
            return mapping;
        }).forEachOrdered((mapping) -> {
            final_MAPPINGS.add(mapping);
        });
        return final_MAPPINGS;
    }

    /**
     * This makes a map of matching atoms out of a map of matching bonds as
     * produced by the get(Subgraph|Ismorphism)Map methods.
     *
     * @param l The list produced by the getMap method.
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchRings
     * @param matchAtomTypes
     * @return The mapping found projected on g1. This is a {@link List} of
     * {@link RMap} objects containing Ids of matching atoms.
     */
    public static Map<Integer, Integer> makeAtomsMapOfBondsMap(
            Map<Integer, Integer> l, IAtomContainer g1, IAtomContainer g2,
            AtomMatcher am, BondMatcher bm) {
        if (l == null) {
            return (new TreeMap<>());
        }
        Map<Integer, Integer> result = new TreeMap<>();
        for (Map.Entry<Integer, Integer> map : l.entrySet()) {
            IBond bond1 = g1.getBond(map.getKey());
            IBond bond2 = g2.getBond(map.getValue());
            IAtom[] atom1 = BondManipulator.getAtomArray(bond1);
            IAtom[] atom2 = BondManipulator.getAtomArray(bond2);
            for (int j = 0; j < 2; j++) {
                List<IBond> bondsConnectedToAtom1j = g1.getConnectedBondsList(atom1[j]);
                for (int k = 0; k < bondsConnectedToAtom1j.size(); k++) {
                    if (!bondsConnectedToAtom1j.get(k).equals(bond1)) {
                        IBond testBond = (IBond) bondsConnectedToAtom1j.get(k);
                        for (Map.Entry<Integer, Integer> m : l.entrySet()) {
                            IBond testBond2;
                            if (m.getKey() == g1.indexOf(testBond)) {
                                testBond2 = g2.getBond(m.getValue());
                                for (int n = 0; n < 2; n++) {
                                    List<IBond> bondsToTest = g2.getConnectedBondsList(atom2[n]);
                                    if (bondsToTest.contains(testBond2)) {

                                        if (j == n) {
                                            if (!result.containsKey(g1.indexOf(atom1[0]))
                                                    && !result.containsValue(g2.indexOf(atom2[0]))) {
                                                result.put(g1.indexOf(atom1[0]), g2.indexOf(atom2[0]));
                                            }
                                        } else {
                                            if (!result.containsKey(g1.indexOf(atom1[1]))
                                                    && !result.containsValue(g2.indexOf(atom2[0]))) {
                                                result.put(g1.indexOf(atom1[1]), g2.indexOf(atom2[0]));
                                            }
                                        }

                                        if (j == n) {
                                            if (!result.containsKey(g1.indexOf(atom1[1]))
                                                    && !result.containsValue(g2.indexOf(atom2[1]))) {
                                                result.put(g1.indexOf(atom1[1]), g2.indexOf(atom2[1]));
                                            }
                                        } else {
                                            if (!result.containsKey(g1.indexOf(atom1[0]))
                                                    && !result.containsValue(g2.indexOf(atom2[1]))) {
                                                result.put(g1.indexOf(atom1[0]), g2.indexOf(atom2[1]));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (result.isEmpty() && l.size() == 1) {
            result = SingleMappingCase(l, g1, g2, am, bm);
        }
        return result;
    }

    /**
     *
     * @param ac1
     * @param mapping
     * @return
     * @throws CloneNotSupportedException
     */
    public static IAtomContainer getSubgraphProjectAtoms(IAtomContainer ac1, Set<Integer> mapping) throws CloneNotSupportedException {
        IAtomContainer ac = ac1.clone();
        Set<IAtom> remove = new HashSet<>();
        for (IAtom a : ac1.atoms()) {
            if (!mapping.contains(ac1.indexOf(a))) {
                remove.add(a);
            }
        }

        remove.forEach((a) -> {
            ac.removeAtom(a);
        });

        System.out.println("");
        return ac;
    }

    /**
     *
     * @param comp_graph_nodes
     * @param ac
     * @param mapping
     * @return
     * @throws CloneNotSupportedException
     */
    public static IAtomContainer getSubgraphProjectBonds(
            Graph comp_graph_nodes,
            IAtomContainer ac, Collection<Integer> mapping) throws CloneNotSupportedException {
        IAtomContainer result = ac.clone();

        Set<Integer> commonAtoms = new HashSet<>();
        mapping.forEach((b1) -> {
            mapping.stream().map((b2) -> commonVertices(ac, ac.getBond(b1), ac.getBond(b2))).filter((commonVertices) -> (!commonVertices.isEmpty())).forEachOrdered((commonVertices) -> {
                commonAtoms.addAll(commonVertices);
            });
        });
//        System.out.println("Common Index " + commonAtoms);
        Set<IAtom> removeAtoms = new HashSet<>();
        for (IAtom a : result.atoms()) {
            if (!commonAtoms.contains(result.indexOf(a))) {
                removeAtoms.add(a);
            }
        }

        removeAtoms.forEach((a) -> {
            result.removeAtom(a);
        });

        if (DEBUG) {
            System.out.println("Number of atoms mapped " + commonAtoms.size());
            System.out.println("New AC " + result.getAtomCount());
        }
        return result;
    }

    /**
     * Returns a set with the common vertices of edge E1 and E2 in Graph g The
     * result will be a Set of size 0, 1 or 2
     *
     * @param ac
     * @param e1
     * @param e2
     * @return
     */
    public static Set<Integer> commonVertices(IAtomContainer ac, IBond e1, IBond e2) {
        Set<Integer> commonVertices = new LinkedHashSet<>();

        for (IAtom a : e1.atoms()) {
            for (IAtom b : e2.atoms()) {
                if (a == b) {
                    commonVertices.add(ac.indexOf(a));
                }
            }
        }

        return commonVertices;
    }

    private static Map<Integer, Integer> SingleMappingCase(
            Map<Integer, Integer> l, IAtomContainer g1, IAtomContainer g2,
            AtomMatcher am, BondMatcher bm) {
        Map<Integer, Integer> result = new TreeMap<>();
        if (l.size() == 1) {
            IBond bond1 = g1.getBond(l.keySet().iterator().next());
            IBond bond2 = g2.getBond(l.values().iterator().next());

            if (bond1 instanceof IQueryBond) {
                if (((IQueryBond) bond1).matches(bond2)) {
                    IQueryAtom atom1 = (IQueryAtom) (bond1.getAtom(0));
                    IQueryAtom atom2 = (IQueryAtom) (bond1.getAtom(1));
                    if (atom1.matches(bond2.getAtom(0)) && atom2.matches(bond2.getAtom(1))) {
                        result.put(g1.indexOf(bond1.getAtom(0)), g2.indexOf(bond2.getAtom(0)));
                        result.put(g1.indexOf(bond1.getAtom(1)), g2.indexOf(bond2.getAtom(1)));
                    }
                    if (atom1.matches(bond2.getAtom(1)) && atom2.matches(bond2.getAtom(0))) {
                        result.put(g1.indexOf(bond1.getAtom(0)), g2.indexOf(bond2.getAtom(1)));
                        result.put(g1.indexOf(bond1.getAtom(1)), g2.indexOf(bond2.getAtom(0)));
                    }

                }
            } else {

                IAtom a1 = bond1.getBegin();
                IAtom a2 = bond1.getEnd();
                IAtom b1 = bond2.getBegin();
                IAtom b2 = bond2.getEnd();

                if (AtomBondMatcher.matches(a1, b1, am)
                        && AtomBondMatcher.matches(a2, b2, am)) {
                    result.put(g1.indexOf(bond1.getAtom(0)), g2.indexOf(bond2.getAtom(0)));
                    result.put(g1.indexOf(bond1.getAtom(1)), g2.indexOf(bond2.getAtom(1)));
                }
                if (AtomBondMatcher.matches(a1, b2, am)
                        && AtomBondMatcher.matches(a2, b1, am)) {
                    result.put(g1.indexOf(bond1.getAtom(0)), g2.indexOf(bond2.getAtom(1)));
                    result.put(g1.indexOf(bond1.getAtom(1)), g2.indexOf(bond2.getAtom(0)));
                }
            }
        }

        return result;
    }
}
