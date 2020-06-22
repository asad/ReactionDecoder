/* Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.algorithm.mcgregor.McGregor;
import org.openscience.smsd.graph.EdgeProductGraph;
import org.openscience.smsd.graph.EdgeType;
import org.openscience.smsd.graph.Graph;
import org.openscience.smsd.graph.IClique;
import org.openscience.smsd.graph.Vertex;
import org.openscience.smsd.graph.algorithm.GraphKoch;
import org.openscience.smsd.tools.IterationManager;

/**
 * This class handles MCS plus algorithm which is a combination of c-clique
 * algorithm and McGregor algorithm.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public final class MCSPlus {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MCSPlus.class);
    private final static boolean DEBUG = false;
    private final IAtomContainer ac1;
    private final IAtomContainer ac2;
    private final List<List<Integer>> overlaps;

    private boolean timeout = false;

    private IterationManager iterationManager = null;
    private final AtomMatcher atomMatcher;
    private final BondMatcher bondMatcher;

    /**
     * @return the timeout
     */
    public synchronized boolean isTimeout() {
        return timeout;
    }

    /**
     * @return the iterationManager
     */
    private IterationManager getIterationManager() {
        return iterationManager;
    }

    /**
     * @param iterationManager the iterationManager to set
     */
    private void setIterationManager(IterationManager iterationManager) {
        this.iterationManager = iterationManager;
    }

    /**
     *
     * @param shouldMatchRings
     * @param shouldMatchBonds
     * @param ac1
     * @param ac2
     * @param matchAtomType
     */
    public MCSPlus(IAtomContainer ac1,
            IAtomContainer ac2,
            AtomMatcher am,
            BondMatcher bm) {
        this.atomMatcher = am;
        this.bondMatcher = bm;

        this.ac1 = ac1;
        this.ac2 = ac2;
        this.overlaps = calculateMCS();
    }

    /**
     *
     * @param ac1
     * @param ac2
     */
    public MCSPlus(IQueryAtomContainer ac1,
            IAtomContainer ac2,
            AtomMatcher am,
            BondMatcher bm) {
        this.atomMatcher = am;
        this.bondMatcher = bm;

        this.ac1 = ac1;
        this.ac2 = ac2;
        this.overlaps = calculateMCS();
    }

    /**
     *
     * @param ac1
     * @param ac2
     * @return
     * @throws CDKException
     */
    private List<List<Integer>> calculateMCS() {

        List<List<Integer>> finalMapping = new ArrayList<>();

        if (DEBUG) {
            System.out.println("ac1 : " + ac1.getAtomCount());
            System.out.println("ac2 : " + ac2.getAtomCount());
        }
        setIterationManager(new IterationManager((ac1.getAtomCount() * ac2.getAtomCount())));
        try {
            EdgeProductGraph gcg
                    = EdgeProductGraph.create(ac1, ac2, atomMatcher, bondMatcher);
            int search_cliques = gcg.searchCliques();
            Graph comp_graph_nodes = gcg.getCompatibilityGraph();
            if (DEBUG) {
                System.out.println("**************************************************");
                System.out.println("--Compatibility Graph--");
                System.out.println("C_edges: " + comp_graph_nodes.getCEdges().size());
                System.out.println("D_edges: " + comp_graph_nodes.getDEdges().size());
                System.out.println("Vertices: " + comp_graph_nodes.V());
                System.out.println("Edges: " + comp_graph_nodes.E());
                System.out.println("**************************************************");
            }

            IClique init = null;
            boolean connected = ConnectivityChecker.isConnected(ac1)
                    && ConnectivityChecker.isConnected(ac2);
            init = new GraphKoch(comp_graph_nodes);
            init.findMaximalCliques();

            Stack<Set<Vertex>> maxCliqueSet = init.getMaxCliquesSet();
            if (DEBUG) {
                System.out.println("Max_Cliques_Set: " + maxCliqueSet.size());
                System.out.println("**************************************************");
            }
            List<Map<Integer, Integer>> mappings = new ArrayList<>();

            while (!maxCliqueSet.empty()) {
                Map<Integer, Integer> indexindexMapping;
                indexindexMapping = MappingHandler.getMapping(comp_graph_nodes, ac1, ac2, maxCliqueSet.peek(),
                        atomMatcher, bondMatcher);
                if (indexindexMapping != null) {
                    mappings.add(indexindexMapping);
                    if (DEBUG) {
                        System.out.println("mappings " + mappings);
                    }
                }
                maxCliqueSet.pop();
            }

            //clear all the compatibility graph content
            gcg.clear();

            for (Map<Integer, Integer> m : mappings) {
                //find mapped atoms of both molecules and store these in mappedAtoms
                List<Integer> exact_mapped_atoms = new ArrayList<>();
                if (DEBUG) {
                    System.out.println("\nClique Mapped Atoms");
                }
                m.entrySet().stream().map((map) -> {
                    if (DEBUG) {
                        System.out.println("i:" + map.getKey() + " j:" + map.getValue());
                    }
                    exact_mapped_atoms.add(map.getKey());
                    return map;
                }).forEach((map) -> {
                    exact_mapped_atoms.add(map.getValue());
                });
                finalMapping.add(exact_mapped_atoms);
            }

            if (DEBUG) {
                System.out.println("mappings: " + mappings);
            }
            List<List<Integer>> extendMappings = null;
            if (ac1 instanceof IQueryAtomContainer) {
                extendMappings = searchMcGregorMapping((IQueryAtomContainer) ac1, ac2, mappings);
            } else {
                extendMappings = searchMcGregorMapping(ac1, ac2, mappings);
            }

            if (extendMappings != null && !extendMappings.isEmpty()) {
                finalMapping.addAll(extendMappings);
            }
            if (DEBUG) {
                //int size = !extendMappings.isEmpty() ? (extendMappings.size() / 2) : 0;
                System.out.println("extendMappings: " + extendMappings);
            }
        } catch (IOException ex) {
            LOGGER.error(Level.SEVERE, null, ex);
        }
        return finalMapping;
    }

    private List<List<Integer>> searchMcGregorMapping(
            IAtomContainer ac1,
            IAtomContainer ac2,
            List<Map<Integer, Integer>> allMCSCopy) throws IOException {

        List<List<Integer>> cliques = new ArrayList<>();

        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : allMCSCopy) {
            Map<Integer, Integer> extendMapping = new TreeMap<>(firstPassMappings);
            McGregor mgit;
            if (ac1.getAtomCount() >= ac2.getAtomCount()
                    && extendMapping.size() < ac2.getAtomCount()) {
                if (DEBUG) {
                    System.out.println("McGregor 1");
                }
                mgit = new McGregor(ac1, ac2, cliques, atomMatcher, bondMatcher);
                mgit.startMcGregorIteration(ac1, mgit.getMCSSize(), extendMapping);
                cliques = mgit.getMappings();
            } else if (ac1.getAtomCount() < ac2.getAtomCount()
                    && extendMapping.size() < ac1.getAtomCount()) {
                if (DEBUG) {
                    System.out.println("McGregor 2");
                }
                extendMapping.clear();
                ROPFlag = false;
                firstPassMappings.entrySet().stream().forEach((map) -> {
                    extendMapping.put(map.getValue(), map.getKey());
                });
                mgit = new McGregor(ac2, ac1, cliques, atomMatcher, bondMatcher);
                mgit.startMcGregorIteration(ac2, mgit.getMCSSize(), extendMapping);
                cliques = mgit.getMappings();
            } else {
                ROPFlag = true;
                //find mapped atoms of both molecules and store these in mappedAtoms
                List<Integer> exact_mapped_atoms = new ArrayList<>();
                if (DEBUG) {
                    System.out.println("\nExact Mapped Atoms");
                }
                extendMapping.entrySet().stream().map((map) -> {
                    if (DEBUG) {
                        System.out.println("i:" + map.getKey() + " j:" + map.getValue());
                    }
                    exact_mapped_atoms.add(map.getKey());
                    return map;
                }).forEach((map) -> {
                    exact_mapped_atoms.add(map.getValue());
                });
                cliques.add(exact_mapped_atoms);
            }
            if (checkTimeout()) {
                break;
            }
        }
        List<List<Integer>> finalMappings = setMappings(ROPFlag, cliques);
        if (DEBUG) {
            System.out.println("After MG --First Mapping-- " + finalMappings.get(0).size());
            System.out.println("After set Sol count MG " + finalMappings.size());
        }
        return finalMappings;
    }

    private List<List<Integer>> searchMcGregorMapping(
            IQueryAtomContainer ac1,
            IAtomContainer ac2,
            List<Map<Integer, Integer>> allMCSCopy) throws IOException {

        List<List<Integer>> cliques = new ArrayList<>();

        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : allMCSCopy) {
            Map<Integer, Integer> extendMapping = new TreeMap<>(firstPassMappings);
            McGregor mgit;
            mgit = new McGregor((IQueryAtomContainer) ac1, ac2, cliques, atomMatcher, bondMatcher);
            mgit.startMcGregorIteration((IQueryAtomContainer) ac1, mgit.getMCSSize(), extendMapping);
//            System.out.println("\nStart McGregor search");
            //Start McGregor search
            cliques = mgit.getMappings();
//            System.out.println("\nSol count after MG " + cliquesBondMap.size());
            if (checkTimeout()) {
                break;
            }
        }
        List<List<Integer>> finalMappings = setMappings(ROPFlag, cliques);
//        System.out.println("After set Sol count MG " + finalMappings.size());
        return finalMappings;
    }

    private List<List<Integer>> setMappings(
            boolean RONP,
            List<List<Integer>> mappings) {
        int counter = 0;
        int mcsSize = 0;
        List<List<Integer>> finalMappings = new ArrayList<>();
        for (List<Integer> mapping : mappings) {
            List<Integer> indexindexMapping = new ArrayList<>();
            for (int index = 0; index < mapping.size(); index += 2) {
                Integer qIndex;
                Integer tIndex;

                if (RONP) {
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != null && tIndex != null) {
                    indexindexMapping.add(qIndex);
                    indexindexMapping.add(tIndex);
                }
            }
            if (!indexindexMapping.isEmpty() && indexindexMapping.size() > mcsSize) {
                mcsSize = indexindexMapping.size();
                finalMappings.clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty() && !finalMappings.contains(indexindexMapping)
                    && (indexindexMapping.size()) == mcsSize) {
                finalMappings.add(counter, indexindexMapping);
                counter++;
            }
        }
        return finalMappings;
    }

    private boolean checkTimeout() {
        if (getIterationManager().isMaxIteration()) {
            this.timeout = true;
//            System.out.println("MCS+ iterations " + getIterationManager().getCounter());
            return true;
        }
        getIterationManager().increment();
        return false;
    }

    /**
     * @return the overlaps
     */
    public synchronized List<List<Integer>> getOverlaps() {
        return Collections.unmodifiableList(overlaps);
    }
}
