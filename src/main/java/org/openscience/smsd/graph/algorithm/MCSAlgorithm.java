/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph.algorithm;

import org.openscience.smsd.graph.EdgeProductGraph;
import org.openscience.smsd.graph.Vertex;
import java.io.IOException;
import java.util.Set;
import java.util.Stack;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;

/**
 * This class implements calls MCS graph algorithms
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman @ bioinceptionlabs.com>
 */
public class MCSAlgorithm {

    private final static boolean DEBUG = false;

    /**
     *
     * @param source
     * @param target
     * @param am
     * @param bm
     * @return
     * @throws IOException
     * @throws java.lang.CloneNotSupportedException
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static IAtomContainer koch(IAtomContainer source, IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm) throws IOException, CloneNotSupportedException, CDKException {

        EdgeProductGraph compatibilityGraph = EdgeProductGraph.create(source, target, am, bm);
        compatibilityGraph.searchCliques();
        boolean disconnected = ConnectivityChecker.isConnected(source) && ConnectivityChecker.isConnected(target);

        GraphKoch graphKoch = new GraphKoch(compatibilityGraph.getCompatibilityGraph());
        graphKoch.findMaximalCliques();
        if (DEBUG) {
            System.out.println("graphKoch.getMaxCliquesSet() " + graphKoch.getMaxCliquesSet());
        }
        Stack<Set<Vertex>> maxCliquesSet = graphKoch.getMaxCliquesSet();
        if (maxCliquesSet == null || maxCliquesSet.isEmpty()) { //interrupted
            return null;
        } else {
            if (DEBUG) {
                System.out.println("\nKoch3 Q:" + SmilesGenerator.generic().create(compatibilityGraph.toQuerySubgraph(maxCliquesSet.peek())));
                System.out.println("\nKoch3 P:" + SmilesGenerator.generic().create(compatibilityGraph.toTargetSubgraph(maxCliquesSet.peek())));
            }
            return compatibilityGraph.toQuerySubgraph(maxCliquesSet.peek());
        }
    }
}
