/*
 * Copyright (C) 2014 Syed Asad Rahman <asad at ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package org.openscience.smsd.algorithm.vflib;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.mcsplus.BKKCKCF;
import org.openscience.smsd.algorithm.mcsplus.GenerateCompatibilityGraph;
import org.openscience.smsd.algorithm.rgraph.CDKRMapHandler;
import org.openscience.smsd.interfaces.Algorithm;

/**
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib
 * {@link org.openscience.smsd.algorithm.vflib.VF2MCS} and reports
 * MCS between run source and target graphs. Then these solutions are extended
 * using McGregor {@link org.openscience.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MCSSeedGenerator implements Callable<List<AtomAtomMapping>> {

    private final IAtomContainer source;
    private final IAtomContainer target;
    private final List<AtomAtomMapping> allCliqueAtomMCS;
    private final boolean ringMatch;
    private final Algorithm algorithm;
    private final static ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(MCSSeedGenerator.class);
    private final boolean bondMatch;
    private final boolean matchAtomType;

    /**
     *
     * @param source
     * @param target
     * @param bondMatch
     * @param ringMatch
     * @param matchAtomType
     * @param algorithm
     */
    public MCSSeedGenerator(IAtomContainer source, IAtomContainer target, boolean bondMatch, boolean ringMatch, boolean matchAtomType, Algorithm algorithm) {
        this.source = source;
        this.target = target;
        this.allCliqueAtomMCS = new ArrayList<>();
        this.ringMatch = ringMatch;
        this.algorithm = algorithm;
        this.matchAtomType = matchAtomType;
        this.bondMatch = bondMatch;
    }

    public MCSSeedGenerator(IQueryAtomContainer source, IAtomContainer target, Algorithm algorithm) {
        this.source = source;
        this.target = target;
        this.allCliqueAtomMCS = new ArrayList<>();
        this.ringMatch = true;
        this.algorithm = algorithm;
        this.matchAtomType = true;
        this.bondMatch = true;
    }

    @Override
    public List<AtomAtomMapping> call() throws Exception {
//        System.out.println("ac1: " + this.source.getAtomCount());
//        System.out.println("ac2: " + this.target.getAtomCount());
        switch (algorithm) {
            case CDKMCS:
                //            System.out.println("Calling CDKMCS " + bondMatch + " " + ringMatch);
                List<AtomAtomMapping> addUIT = addUIT();
//            System.out.println("addUIT " + addUIT.iterator().next().getCount());
                return addUIT;
            case MCSPlus:
                //            System.out.println("Calling MCSPLUS " + bondMatch + " " + ringMatch + " " + matchAtomType);
                List<AtomAtomMapping> addKochCliques = addKochCliques();
//            System.out.println("MCSPLUS " + addKochCliques.iterator().next().getCount());
                return addKochCliques;
            default:
                return Collections.unmodifiableList(allCliqueAtomMCS);
        }
    }

    protected synchronized List<AtomAtomMapping> addKochCliques() throws IOException {
        IAtomContainer ac1;
        IAtomContainer ac2;
        boolean flagExchange = false;

        if (source instanceof IQueryAtomContainer) {
            ac1 = (IQueryAtomContainer) source;
            ac2 = target;
        } else if (source.getAtomCount() < target.getAtomCount()) {
            ac1 = source;
            ac2 = target;
        } else {
            flagExchange = true;
            ac1 = target;
            ac2 = source;
        }

        GenerateCompatibilityGraph gcg
                = new GenerateCompatibilityGraph(ac1, ac2, bondMatch, ringMatch, matchAtomType);
        List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();
        List<Integer> cEdges = gcg.getCEgdes();
        List<Integer> dEdges = gcg.getDEgdes();
        BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
        Stack<List<Integer>> maxCliqueSet = new Stack<>();
        maxCliqueSet.addAll(init.getMaxCliqueSet());
        Collections.sort(maxCliqueSet, (List<Integer> a1, List<Integer> a2) -> a2.size() - a1.size() // assumes you want biggest to smallest
        );
        while (!maxCliqueSet.empty()) {
            List<Integer> peek = maxCliqueSet.peek();
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

            for (Integer value : peek) {
                int[] index = getIndex(value, comp_graph_nodes);
                Integer qIndex = index[0];
                Integer tIndex = index[1];
                if (qIndex != -1 && tIndex != -1) {
                    IAtom qAtom;
                    IAtom tAtom;
                    if (flagExchange) {
                        qAtom = source.getAtom(tIndex);
                        tAtom = target.getAtom(qIndex);
                    } else {
                        qAtom = source.getAtom(qIndex);
                        tAtom = target.getAtom(tIndex);
                    }
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }

            if (!atomatomMapping.isEmpty()) {
                allCliqueAtomMCS.add(atomatomMapping);
            }
            maxCliqueSet.pop();
        }
        gcg.clear();
        return Collections.unmodifiableList(allCliqueAtomMCS);
    }

    /**
     *
     * @return
     */
    private List<AtomAtomMapping> addUIT() throws CDKException {
        CDKRMapHandler rmap = new CDKRMapHandler();
        List<Map<Integer, Integer>> solutions;

        boolean rOnPFlag;
        if (source instanceof IQueryAtomContainer) {
            rOnPFlag = false;
            solutions = rmap.calculateOverlapsAndReduce(target, (IQueryAtomContainer) source);
        } else if (source.getAtomCount() > target.getAtomCount()) {
            rOnPFlag = true;
            solutions = rmap.calculateOverlapsAndReduce(source, target, bondMatch, ringMatch, matchAtomType);
        } else {
            rOnPFlag = false;
            solutions = rmap.calculateOverlapsAndReduce(target, source, bondMatch, ringMatch, matchAtomType);
        }
        return setUITMappings(rOnPFlag, solutions);
    }

    private List<AtomAtomMapping> setUITMappings(boolean RONP, List<Map<Integer, Integer>> sol) {
        /*
         * Sort biggest clique to smallest
         */
        Collections.sort(sol, new Map1ValueComparator(SortOrder.DESCENDING));
        sol.stream().map((Map<Integer, Integer> solution) -> {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            solution.keySet().stream().forEach((qAtomIndex) -> {
                IAtom qAtom;
                IAtom tAtom;
                int qIndex;
                int tIndex;
                
                if (RONP) {
                    qAtom = source.getAtom(qAtomIndex);
                    tAtom = target.getAtom(solution.get(qAtomIndex));
                } else {
                    tAtom = target.getAtom(qAtomIndex);
                    qAtom = source.getAtom(solution.get(qAtomIndex));
                }
                
                qIndex = source.indexOf(qAtom);
                tIndex = target.indexOf(tAtom);
                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            });
            return atomatomMapping;
        }).filter((atomatomMapping) -> (!atomatomMapping.isEmpty())).forEach((atomatomMapping) -> {
            allCliqueAtomMCS.add(atomatomMapping);
        });
        return Collections.unmodifiableList(allCliqueAtomMCS);
    }

    private int[] getIndex(int cliqueIndex, List<Integer> comp_graph_nodes) {
        int[] v = new int[2];
        v[0] = -1;
        v[1] = -1;
        for (int i = 0; i < comp_graph_nodes.size(); i += 3) {
            if (cliqueIndex == comp_graph_nodes.get(i + 2)) {
                v[0] = comp_graph_nodes.get(i);
                v[1] = comp_graph_nodes.get(i + 1);
            }
        }
        return v;
    }
}
