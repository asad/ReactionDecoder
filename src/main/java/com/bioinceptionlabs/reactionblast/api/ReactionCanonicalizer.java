/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.api;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Canonical reaction signature generator using Weisfeiler-Lehman (WL) graph
 * hashing on the Imaginary Transition State (ITS) graph.
 *
 * The ITS graph merges reactant and product molecular graphs, with edge labels
 * encoding bond changes (formed, cleaved, order changed). The WL hash produces
 * a canonical, invariant fingerprint that is:
 * - Deterministic (same reaction always gives same hash)
 * - Permutation-invariant (independent of atom ordering)
 * - Hierarchical (deeper iterations capture wider neighborhood)
 *
 * Based on the Weisfeiler-Lehman graph isomorphism test (1968) and its
 * application to molecular graphs. Implementation is IP-free (public domain
 * algorithm, no dependency on external tools like Nauty).
 *
 * References:
 * - Weisfeiler, Lehman (1968): "A reduction of a graph to a canonical form"
 * - Shervashidze et al. (2011): "Weisfeiler-Lehman Graph Kernels" (JMLR)
 * - Leber (2008): R-matrix canonicalization for enzymatic reactions
 * - Phan et al. (2025): SynKit graph-based reaction canonicalization
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class ReactionCanonicalizer {

    private ReactionCanonicalizer() {}

    private static final int WL_ITERATIONS = 3;

    /**
     * Compute a canonical hash for a reaction based on its bond changes.
     * The hash is invariant to atom ordering and deterministic.
     *
     * @param formedCleavedBonds bond formation/cleavage patterns (e.g., "C-O:1")
     * @param orderChangedBonds bond order change patterns (e.g., "C=C:1")
     * @param stereoChangedBonds stereo change patterns
     * @param reactionCentreFP reaction centre fingerprint patterns
     * @return canonical hex hash string (SHA-256 based)
     */
    public static String computeCanonicalHash(
            List<String> formedCleavedBonds,
            List<String> orderChangedBonds,
            List<String> stereoChangedBonds,
            List<String> reactionCentreFP) {

        // Build ITS graph as adjacency representation
        // Nodes = unique atom types at reaction centre
        // Edges = bond changes with labels
        ITSGraph its = buildITSGraph(formedCleavedBonds, orderChangedBonds,
                stereoChangedBonds, reactionCentreFP);

        // Apply WL hash iterations
        String wlHash = wlGraphHash(its, WL_ITERATIONS);

        return wlHash;
    }

    /**
     * Build an Imaginary Transition State graph from bond change fingerprints.
     * The ITS graph encodes the reaction centre as a labeled graph where:
     * - Nodes are atom types involved in changes
     * - Edges are bond changes with labels (FORMED, CLEAVED, ORDER_CHANGE)
     */
    static ITSGraph buildITSGraph(
            List<String> formedCleaved,
            List<String> orderChanges,
            List<String> stereoChanges,
            List<String> reactionCentre) {

        ITSGraph graph = new ITSGraph();

        // Parse bond change patterns: "X-Y:weight" or "X=Y:weight"
        for (String pattern : formedCleaved) {
            addBondChange(graph, pattern, "FC");
        }
        for (String pattern : orderChanges) {
            addBondChange(graph, pattern, "OC");
        }
        for (String pattern : stereoChanges) {
            addBondChange(graph, pattern, "SC");
        }
        for (String pattern : reactionCentre) {
            addBondChange(graph, pattern, "RC");
        }

        return graph;
    }

    /**
     * Parse a bond change pattern like "C-O:1" or "C=C:2" and add to graph.
     */
    private static void addBondChange(ITSGraph graph, String pattern, String changeType) {
        // Strip weight suffix
        int colon = pattern.lastIndexOf(':');
        String bondPattern = colon > 0 ? pattern.substring(0, colon) : pattern;
        String weight = colon > 0 ? pattern.substring(colon + 1) : "1";

        // Parse atom pair from bond pattern: "X-Y", "X=Y", "X#Y", "X%Y", "X@Y"
        String[] atoms = bondPattern.split("[-=#%@]");
        if (atoms.length == 2) {
            // Extract bond type symbol
            char bondType = '-';
            for (char c : bondPattern.toCharArray()) {
                if (c == '-' || c == '=' || c == '#' || c == '%' || c == '@') {
                    bondType = c;
                    break;
                }
            }

            String nodeA = atoms[0].trim();
            String nodeB = atoms[1].trim();
            String edgeLabel = changeType + ":" + bondType + ":" + weight;

            graph.addNode(nodeA);
            graph.addNode(nodeB);
            graph.addEdge(nodeA, nodeB, edgeLabel);
        }
    }

    /**
     * Weisfeiler-Lehman graph hash.
     * Iteratively refines node labels by aggregating sorted neighbor labels.
     * The final hash is a sorted concatenation of all refined labels.
     *
     * @param graph the ITS graph
     * @param iterations number of WL refinement iterations
     * @return canonical hash string
     */
    static String wlGraphHash(ITSGraph graph, int iterations) {
        if (graph.nodes.isEmpty()) return "EMPTY";

        // Initial labels = node type (atom symbol)
        Map<String, String> labels = new HashMap<>();
        for (String node : graph.nodes.keySet()) {
            labels.put(node, graph.nodes.get(node));
        }

        // Collect multiset labels at each iteration
        List<String> allLabels = new ArrayList<>();

        for (int iter = 0; iter < iterations; iter++) {
            Map<String, String> newLabels = new HashMap<>();

            for (String node : graph.nodes.keySet()) {
                // Get sorted neighbor labels with edge labels
                List<String> neighborInfo = new ArrayList<>();
                for (ITSGraph.Edge edge : graph.getEdges(node)) {
                    String neighborLabel = labels.get(edge.target);
                    neighborInfo.add(edge.label + "|" + neighborLabel);
                }
                Collections.sort(neighborInfo);

                // New label = old label + sorted neighbor info
                String newLabel = labels.get(node) + "(" + String.join(",", neighborInfo) + ")";
                newLabels.put(node, newLabel);
            }

            labels = newLabels;

            // Collect all labels at this iteration
            List<String> iterLabels = new ArrayList<>(labels.values());
            Collections.sort(iterLabels);
            allLabels.addAll(iterLabels);
        }

        // Final canonical string = sorted concatenation of all iteration labels
        Collections.sort(allLabels);
        String canonical = String.join(";", allLabels);

        // Hash to fixed-length string
        return sha256Hex(canonical);
    }

    /**
     * SHA-256 hash of a string, returned as hex.
     */
    private static String sha256Hex(String input) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            byte[] hash = md.digest(input.getBytes());
            StringBuilder hex = new StringBuilder();
            for (byte b : hash) {
                hex.append(String.format("%02x", b));
            }
            return hex.toString();
        } catch (NoSuchAlgorithmException e) {
            // SHA-256 is always available in Java
            throw new RuntimeException(e);
        }
    }

    /**
     * Internal ITS graph representation.
     * Nodes are atom types, edges are labeled bond changes.
     */
    static class ITSGraph {
        final Map<String, String> nodes = new TreeMap<>(); // id → label
        final List<Edge> edges = new ArrayList<>();

        void addNode(String id) {
            nodes.putIfAbsent(id, id);
        }

        void addEdge(String source, String target, String label) {
            edges.add(new Edge(source, target, label));
            edges.add(new Edge(target, source, label)); // undirected
        }

        List<Edge> getEdges(String node) {
            List<Edge> result = new ArrayList<>();
            for (Edge e : edges) {
                if (e.source.equals(node)) result.add(e);
            }
            return result;
        }

        static class Edge {
            final String source, target, label;
            Edge(String source, String target, String label) {
                this.source = source;
                this.target = target;
                this.label = label;
            }
        }
    }
}
