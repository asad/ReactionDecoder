/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.io.Serializable;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.smsd.tools.Utility;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class EdgeProductCompatibilityGraph implements Serializable {

    /**
     * @return the c_edges
     */
    public Set<Edge> getCEdges() {
        return cEdges;
    }

    /**
     * @return the d_edges
     */
    public Set<Edge> getDEdges() {
        return dEdges;
    }

    /**
     * @return the Compatibility Graph
     */
    public Graph getCompatibilityGraph() {
        return g;
    }

    private static final long serialVersionUID = 96986606860861L;
    private final Set<Edge> cEdges;
    private final Set<Edge> dEdges;
    private final Graph g;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchBonds;
    private final boolean shouldMatchRings;
    private final boolean matchAtomType;
    private final boolean DEBUG = false;

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
    public EdgeProductCompatibilityGraph(
            IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws IOException {
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.matchAtomType = matchAtomType;
        this.source = source;
        this.target = target;
        this.cEdges = Collections.synchronizedSet(new HashSet<>());
        this.dEdges = Collections.synchronizedSet(new HashSet<>());
        this.g = new Graph();

    }

    public int searchCliques() {
        compatibilityGraphNodes();
        int edges = compatibilityGraph();
        if (DEBUG) {
            System.out.println("**************************************************");
            System.out.println("--Compatibility Graph--");
            System.out.println("C_edges: " + getCEdges().size());
            System.out.println("D_edges: " + getDEdges().size());
            System.out.println("Vertices: " + getCompatibilityGraph().V());
            System.out.println("Edges: " + edges);
        }
        return getCompatibilityGraph().V();
    }

    private void compatibilityGraphNodes() {
        int compatibilityNodeCounter = 0;
        for (IBond a : source.bonds()) {
            for (IBond b : target.bonds()) {
                if (Utility.isMatchFeasible(a, b, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                    Vertex node = new Vertex(compatibilityNodeCounter);
                    if (DEBUG) {
                        System.out.println(a.getBegin().getSymbol() + "- 1 -" + a.getEnd().getSymbol());
                        System.out.println(b.getBegin().getSymbol() + "- 2 -" + b.getEnd().getSymbol());
                    }
                    node.setCompatibilityBondPair(source.indexOf(a), target.indexOf(b));
                    g.addNode(node);
                    compatibilityNodeCounter++;
                }
            }
        }
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraph() {

        for (Vertex n1 : g.nodes()) {
            for (Vertex n2 : g.nodes()) {
                //avoid self loop
                if (n1 == n2) {
                    continue;
                }
                EdgeType edgetype = edgePairsCompatible(n1, n2);
                if (edgetype != null) {
                    Edge edge = new Edge(n1, n2);
                    edge.setEdgeType(edgetype);
                    if (edgetype == EdgeType.C_EDGE) {
                        cEdges.add(edge);
                    }
                    if (edgetype == EdgeType.D_EDGE) {
                        dEdges.add(edge);
                    }
                    g.addEdge(edge, false);
                }
            }
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
        e1 = source.getBond(p1.getQueryBond()); //Edge in G1
        e2 = target.getBond(p1.getTargetBond()); //Edge in G2
        f1 = source.getBond(p2.getQueryBond()); //Edge in G1
        f2 = target.getBond(p2.getTargetBond()); //Edge in G2

        //check condition 1)
        if (e1 == f1 || e2 == f2) {
            //condition 1 not satisfied, edges are not compatible
            return null;
        }

        Set<IAtom> possibleVerticesG1 = commonVertices(source, e1, f1);
        Set<IAtom> possibleVerticesG2 = commonVertices(target, e2, f2);

        if (possibleVerticesG1.isEmpty() && possibleVerticesG2.isEmpty()) {
            //e1,f1 and e2,f2 are not adjacent in G1 and in G2, respectively
            //Create a D_Edge
            return EdgeType.D_EDGE;
        }

        for (IAtom v1 : possibleVerticesG1) {
            for (IAtom v2 : possibleVerticesG2) {
//                if (DefaulAtomBondMatcher.isAtomMatch(
//                        DefaulAtomBondMatcher.AtomMatcher(v1, shouldMatchRings, matchAtomType), v2)) {
                if (v1.getSymbol().equals(v2.getSymbol())) {
                    // e1,f1 in G1 are connected via a vertex of
                    // the same label as the vertex shared by e2,f2 in G2.
                    //A C_edge shuold be created
                    return EdgeType.C_EDGE;
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

        for (IAtom a : e1.atoms()) {
            for (IAtom b : e2.atoms()) {
                if (a == b) {
                    commonVertices.add(a);
                }
            }
        }

        return commonVertices;
    }

    /**
     * Clear maps
     */
    public void clear() {
        g.clear();
        cEdges.clear();
        dEdges.clear();
    }
}
