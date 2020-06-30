/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.Stack;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class EdgeProductGraph implements Serializable {

    private final boolean DEBUG = false;

    /**
     * Generates a compatibility graph between two molecules
     *
     * @param source
     * @param target
     * @param am
     * @param bm
     * @return
     * @throws IOException
     */
    public static EdgeProductGraph
            create(IAtomContainer source,
                    IAtomContainer target,
                    AtomMatcher am,
                    BondMatcher bm)
            throws IOException {
        return new EdgeProductGraph(source, target, am, bm);
    }
    private final AtomMatcher atomMatcher;
    private final BondMatcher bondMatcher;

    /**
     * @return the Compatibility Graph
     */
    public Graph getCompatibilityGraph() {
        return g;
    }

    private static final long serialVersionUID = 96986606860861L;
    private final Graph g;
    private final IAtomContainer source;
    private final IAtomContainer target;

    /**
     * Generates a compatibility graph between two molecules
     *
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @throws java.io.IOException
     */
    private EdgeProductGraph(
            IAtomContainer source,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm) throws IOException {
        this.atomMatcher = am;
        this.bondMatcher = bm;
        this.source = source;
        this.target = target;
        this.g = new Graph();

    }

    public int searchCliques() {
        compatibilityGraphNodes();
        int edges = compatibilityGraphDirected();
        if (DEBUG) {
            System.out.println("**************************************************");
            System.out.println("--Compatibility Graph--");
            System.out.println("C_edges: " + g.getCEdges().size());
            System.out.println("D_edges: " + g.getDEdges().size());
            System.out.println("Vertices: " + g.V());
            System.out.println("Edges: " + g.E());
        }
        return g.V();
    }

    private void compatibilityGraphNodes() {
        int compatibilityNodeCounter = 1;
        Iterable<IBond> qbonds = source.bonds();
        Iterable<IBond> tbonds = target.bonds();
        for (IBond a : qbonds) {
            for (IBond b : tbonds) {
                //Asad-Imp for large graphs
                //Only add the edge product vertex if the edge labels and vertex labels are the same
                //IMP: directed manner i.e. if {a-b = a-b} then true else false 
                //Only add the edge product vertex if the edge labels and end vertex labels are the same
                if (AtomBondMatcher.matchAtomAndBond(a, b, atomMatcher, bondMatcher, true)) {
                    Vertex node = new Vertex(compatibilityNodeCounter);
                    if (DEBUG) {
                        System.out.print("Q: " + source.indexOf(a) + ", " + a.getBegin().getSymbol() + "- 1 -" + a.getEnd().getSymbol());
                        System.out.println(", T: " + target.indexOf(b) + ", " + b.getBegin().getSymbol() + "- 2 -" + b.getEnd().getSymbol());
                    }
                    node.setCompatibilityBondPair(source.indexOf(a), target.indexOf(b));
                    g.addNode(node);
                    compatibilityNodeCounter++;

                }
            }
        }

        if (DEBUG) {
            System.out.println("Vertices " + g.V());
        }
    }

    private void addEdge(Vertex n1, Vertex n2) {
        EdgeType edgetype = edgePairsCompatible(n1, n2);
        if (edgetype != null) {
            if (DEBUG) {
                System.out.println("n1: " + n1.getID()
                        + ", " + "n2: " + n2.getID() + ", Edge " + edgetype);
            }
            if (edgetype == EdgeType.C_EDGE
                    || edgetype == EdgeType.D_EDGE) {
                //Assume it to be a undirected graph
                g.addEdge(n1, n2, edgetype);
            }
        }
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphDirected() {
        int counter = 1;
        Stack<Vertex> nodesToCompare = new Stack<>();
        Iterator<Vertex> iterator = g.iterator();
        while (iterator.hasNext()) {
            nodesToCompare.add(iterator.next());
        }
        while (!nodesToCompare.empty()) {
            Vertex n1 = nodesToCompare.pop();
            for (Vertex n2 : nodesToCompare) {
                addEdge(n1, n2);
                if (DEBUG && counter % 100000 == 0) {
                    System.out.println(nodesToCompare.size() + "    Found clique #" + counter);
                    System.out.println("    C-Edges #" + g.getCEdges().size() + ".\n");
                    System.out.println("    D-Edges #" + g.getDEdges().size() + ".\n");
                    System.out.println("    Edges #" + g.E() + ".\n");
                }
                counter++;
            }
        }

        if (DEBUG) {
            System.out.println("Edges " + g.E());
        }
        return g.E();
    }

    /**
     * Returns true when two edge pairs (e1,e2) and (f1,f2) are compatible
     *
     * There is an edge between two vertices eH,fH in VH with eH =(e1,e2) and fH
     * =(f1,f2), if 1) e1 != f1 and e2 != f2, and 2) if either e1,f1 in G1 are
     * connected via a vertex of the same label as the vertex shared by e2,f2 in
     * G2, 3) or e1,f1 and e2,f2 are not adjacent in G1 and in G2, respectively
     *
     */
    private EdgeType edgePairsCompatible(Vertex p1, Vertex p2) {
        //either e1,f1 in G1 are connected via a vertex of the same label as the vertex shared by e2,f2 in G2
        //or e1,f1 and e2,f2 are not adjacent in G1 and in G2, respectively
        IBond e1, e2, f1, f2;
        e1 = source.getBond(p1.getQueryBondIndex()); //Edge in G1
        e2 = target.getBond(p1.getTargetBondIndex()); //Edge in G2
        f1 = source.getBond(p2.getQueryBondIndex()); //Edge in G1
        f2 = target.getBond(p2.getTargetBondIndex()); //Edge in G2

        //check condition 1)
        if (e1 == f1 || e2 == f2) {
            //condition 1 not satisfied, edges are not compatible
            return null;
        }

        Set<IAtom> possibleVerticesG1 = commonVertices(source, e1, f1);
        Set<IAtom> possibleVerticesG2 = commonVertices(target, e2, f2);
        if (DEBUG) {
            System.out.println("possibleVerticesG1 " + possibleVerticesG1.size());
            System.out.println("possibleVerticesG2 " + possibleVerticesG2.size());
        }
        if (possibleVerticesG1.isEmpty() && possibleVerticesG2.isEmpty()) {
            //e1,f1 and e2,f2 are not adjacent in G1 and in G2, respectively
            //Create a D_Edge
            return EdgeType.D_EDGE;
        }
        if (!possibleVerticesG1.isEmpty() && !possibleVerticesG2.isEmpty()) {
            for (IAtom v1 : possibleVerticesG1) {
                for (IAtom v2 : possibleVerticesG2) {
                    if (AtomBondMatcher.matches(v1, v2, atomMatcher)) {
                        // e1,f1 in G1 are connected via a vertex of
                        // the same label as the vertex shared by e2,f2 in G2.
                        //A C_edge should be created
                        return EdgeType.C_EDGE;
                    }
                }
            }
        }

        //The edge pairs are not compatible
        return null;
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
    public Set<IAtom> commonVertices(IAtomContainer ac, IBond e1, IBond e2) {
        Set<IAtom> commonVertices = new LinkedHashSet<>();
        if (e1.getBegin().equals(e2.getBegin())) {
            commonVertices.add(e1.getBegin());
        }
        if (e1.getBegin().equals(e2.getEnd())) {
            commonVertices.add(e1.getBegin());
        }

        if (e1.getEnd().equals(e2.getBegin())) {
            commonVertices.add(e1.getEnd());
        }
        if (e1.getEnd().equals(e2.getEnd())) {
            commonVertices.add(e1.getEnd());
        }

        return commonVertices;
    }

    /**
     * Creates the subgraph of g1 containing all the edges from the edge product
     * in the vertices of this EdgeProductGraph
     *
     * @param edgeProductVertices if (and only if) these vertices induce a
     * complete subgraph in this EdgeProductGraph, then the result will be the a
     * common subgraph of g1 and g2.
     * @return a subgraph of g1
     * @throws java.lang.CloneNotSupportedException
     */
    public IAtomContainer toQuerySubgraph(Set<Vertex> edgeProductVertices) throws CloneNotSupportedException {

        IAtomContainer ac = ExtAtomContainerManipulator.cloneWithIDs(source);

        //Add the left Edge (including vertices) from all the EdgeProducts in vertices
        Set<IAtom> atomsMapped = new HashSet<>();
        edgeProductVertices.stream().map((ep) -> ep.getQueryBondIndex()).map((bondIndex) -> ac.getBond(bondIndex)).map((bond) -> {
            atomsMapped.add(bond.getBegin());
            return bond;
        }).forEachOrdered((bond) -> {
            atomsMapped.add(bond.getEnd());
        });
        Set<IAtom> atomsToBeRemoved = new HashSet<>();
        for (IAtom a : ac.atoms()) {
            atomsToBeRemoved.add(a);
        }

        atomsToBeRemoved.removeAll(atomsMapped);
        atomsToBeRemoved.forEach((a) -> {
            ac.removeAtom(a);
        });

        return ac;
    }

    /**
     * Creates the subgraph of g2 containing all the edges from the edge product
     * in the vertices of this EdgeProductGraph
     *
     * @param edgeProductVertices if (and only if) these vertices induce a
     * complete subgraph in this EdgeProductGraph, then the result will be the a
     * common subgraph of g2 and g1.
     * @return a subgraph of g2
     * @throws java.lang.CloneNotSupportedException
     */
    public IAtomContainer toTargetSubgraph(Set<Vertex> edgeProductVertices) throws CloneNotSupportedException {

        IAtomContainer ac = ExtAtomContainerManipulator.cloneWithIDs(target);

        //Add the left Edge (including vertices) from all the EdgeProducts in vertices
        Set<IAtom> atomsMapped = new HashSet<>();
        edgeProductVertices.stream().map((ep)
                -> ep.getTargetBondIndex()).map((bondIndex)
                -> ac.getBond(bondIndex)).map((bond) -> {
            atomsMapped.add(bond.getBegin());
            return bond;
        }).forEachOrdered((bond) -> {
            atomsMapped.add(bond.getEnd());
        });
        Set<IAtom> atomsToBeRemoved = new HashSet<>();
        for (IAtom a : ac.atoms()) {
            atomsToBeRemoved.add(a);
        }

        atomsToBeRemoved.removeAll(atomsMapped);
        atomsToBeRemoved.forEach((a) -> {
            ac.removeAtom(a);
        });

        return ac;
    }

    /**
     * Clear maps
     */
    public void clear() {
        g.clear();
    }
}
