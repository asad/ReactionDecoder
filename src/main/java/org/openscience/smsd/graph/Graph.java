/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 */
public final class Graph implements Iterable<Vertex> {

    private static final String NEWLINE = System.getProperty("line.separator");

    private final Map<Vertex, Set<Vertex>> adj;
    private final Map<Vertex, Set<Vertex>> c_adj;
    private final Map<Vertex, Set<Vertex>> d_adj;
    private final Map<EdgeType, Set<Edge>> adj_type_Map;
    private final List<Vertex> vertices;

    /**
     * Initializes an empty graph with {@code V} vertices and 0 edges.param V
     * the number of vertices
     *
     */
    public Graph() {
        this.vertices = new ArrayList<>();
        this.adj = new TreeMap<>();
        this.c_adj = new TreeMap<>();
        this.d_adj = new TreeMap<>();
        this.adj_type_Map = new HashMap<>();
    }

    /**
     * Returns the number of vertices in this graph.
     *
     * @return the number of vertices in this graph
     */
    public int V() {
        return vertices.size();
    }

    /**
     * Returns the number of edges in this graph.
     *
     * @return the number of edges in this graph
     */
    public int E() {
        return edges().size();
    }

    /**
     * Returns Nodes in this graph.
     *
     * @return vertices in this graph
     */
    public Set<Vertex> nodes() {
        Set<Vertex> nodes = new HashSet<>();
        nodes.addAll(vertices);
        return nodes;
    }

    /**
     * Returns edges in this graph.
     *
     * @return edges in this graph
     */
    public Set<Edge> edges() {
        Set<Edge> edgesSet = new HashSet<>();
        adj_type_Map.values().forEach((edges) -> {
            edgesSet.addAll(edges);
        });
        return edgesSet;
    }

    private void validateVertex(Vertex v) {
        if (!vertices.contains(v)) {
            throw new IllegalArgumentException("vertex " + v + " not found in the graph");
        }
    }

    public void addEdge(Vertex v, Vertex u, EdgeType e) {

        validateVertex(v);
        validateVertex(u);
        Edge edge = new Edge(vertices.indexOf(v), vertices.indexOf(u));
        edge.setEdgeType(e);
        addEdge(edge);
    }

    /**
     * Adds the undirected edge v-w to this graph. Assumes that the nodes
     * assigned in the edge is already present
     *
     * @param e edge to be added
     */
    private void addEdge(Edge e) {

        /*
         * Add edges to the map
         */
        addEdge(adj, e);

        /*
         * Add C edges to the map
         */
        if (e.getEdgeType() == EdgeType.C_EDGE) {
            addEdge(c_adj, e);
        }
        /*
         * Add D edges to the map
         */
        if (e.getEdgeType() == EdgeType.D_EDGE) {
            addEdge(d_adj, e);
        }
        /*
         * Add Edge type to the map
         */
        if (!adj_type_Map.containsKey(e.getEdgeType())) {
            adj_type_Map.put(e.getEdgeType(), new HashSet<>());
        }
        adj_type_Map.get(e.getEdgeType()).add(e);

    }

    private void addEdge(Map map, Edge e) {
        addEdge(map, vertices.get(e.getSource()), vertices.get(e.getSink()));
        addEdge(map, vertices.get(e.getSink()), vertices.get(e.getSource()));

    }

    private void addEdge(Map<Vertex, Set<Vertex>> map, Vertex u, Vertex v) {
        if (!map.containsKey(u)) {
            map.put(u, new HashSet<>());
        }
        map.get(u).add(v);
    }

    /**
     * Adds Vertex to this graph.
     *
     * @param node Vertex to be added
     */
    public void addNode(Vertex node) {
        if (!adj.containsKey(node)) {
            adj.put(node, new HashSet<>());
            vertices.add(node);
        } else {
            throw new IllegalArgumentException("Node " + node + " found in the graph");
        }
    }

    /**
     * Returns the vertices adjacent to vertex {@code v}.
     *
     * @param v the vertex
     * @return the vertices adjacent to vertex {@code v}, as an iterable
     */
    public Set<Vertex> getNeighbours(Vertex v) {
        validateVertex(v);
        return new TreeSet<>(adj.get(v));
    }

    /**
     * Returns the getDegree of vertex {@code v}.
     *
     * @param v the vertex
     * @return the getDegree of vertex {@code v}
     */
    public int getDegree(Vertex v) {
        validateVertex(v);
        return adj.get(v).size();
    }

    /**
     * Returns a string representation of this graph.
     *
     * @return the number of vertices <em>V</em>, followed by the number of
     * edges <em>E</em>, followed by the <em>V</em> adjacency lists
     */
    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(vertices.size()).append(" vertices, ").append(edges().size()).append(" edges ").append(NEWLINE);
        adj.entrySet().stream().map((m) -> {
            s.append(m.getKey()).append(": ");
            return m;
        }).map((m) -> {
            m.getValue().forEach((w) -> {
                s.append(w).append(" ");
            });
            return m;
        }).forEachOrdered((_item) -> {
            s.append(NEWLINE);
        });
        return s.toString();
    }

    /**
     * Clean graph
     */
    public void clear() {
        this.vertices.clear();
        this.adj.clear();
        this.c_adj.clear();
        this.d_adj.clear();
        this.adj_type_Map.clear();
    }

    /**
     *
     * @param u
     * @param v
     * @return if an edge exists between vertex
     */
    public boolean hasEdge(Vertex u, Vertex v) {
        return adj.containsKey(u) && adj.get(u).contains(v) ? true
                : adj.containsKey(v) && adj.get(v).contains(u);
    }

    /**
     * Returns edges of the vertex
     *
     * @param currentVertex
     * @return
     */
    public Iterable<Edge> edgesOf(Vertex currentVertex) {
        validateVertex(currentVertex);
        Integer v = vertices.indexOf(currentVertex);
        Set<Edge> edgesOfVertex = new LinkedHashSet<>();
        edges().stream().map((e) -> {
            if (e.getSource().equals(v)) {
                edgesOfVertex.add(e);
            }
            return e;
        }).filter((e) -> (e.getSink().equals(v))).forEachOrdered((e) -> {
            edgesOfVertex.add(e);
        });
        return edgesOfVertex;
    }

    /**
     * Returns true if there is c edge else false
     *
     * @param u
     * @param v
     * @return true if there is c edge else false
     */
    public boolean isCEdge(Vertex u, Vertex v) {
        validateVertex(u);
        validateVertex(v);
        return c_adj.containsKey(u) && c_adj.get(u).contains(v) ? true
                : c_adj.containsKey(v) && c_adj.get(v).contains(u);
    }

    /**
     * Returns true if there is d edge else false
     *
     * @param u
     * @param v
     * @return true if there is d edge else false
     */
    public boolean isDEdge(Vertex u, Vertex v) {
        validateVertex(u);
        validateVertex(v);
        return d_adj.containsKey(u) && d_adj.get(u).contains(v) ? true
                : d_adj.containsKey(v) && d_adj.get(v).contains(u);
    }

    /**
     * Returns an edge connecting source vertex to target vertex if such
     * vertices and such edge exist in this graph. Otherwise returns null. If
     * any of the specified vertices is null returns null In undirected graphs,
     * the returned edge may have its source and target vertices in the opposite
     * order.
     *
     * @param edge
     * @return
     */
    public Vertex getEdgeSource(Edge edge) {
        return vertices.get(edge.getSource());
    }

    /**
     * Returns the target vertex of an edge. For an undirected graph, source and
     * target are distinguishable designations (but without any mathematical
     * meaning)
     *
     * @param edge
     * @return
     */
    public Vertex getEdgeTarget(Edge edge) {
        return vertices.get(edge.getSink());
    }

    /**
     *
     * @param v
     * @return
     */
    public boolean removeVertex(Vertex v) {
        adj.keySet().stream().filter((key) -> (!adj.get(key).isEmpty()
                && adj.get(key).contains(v))).forEachOrdered((key) -> {
            adj.get(key).remove(v);
        });
        c_adj.keySet().stream().filter((key) -> (!c_adj.get(key).isEmpty()
                && c_adj.get(key).contains(v))).forEachOrdered((key) -> {
            c_adj.get(key).remove(v);
        });

        d_adj.keySet().stream().filter((key) -> (!d_adj.get(key).isEmpty()
                && d_adj.get(key).contains(v))).forEachOrdered((key) -> {
            d_adj.get(key).remove(v);
        });

        adj_type_Map.entrySet().forEach((c) -> {
            c.getValue().stream().filter((e) -> (vertices.get(e.getSource()) == v
                    || vertices.get(e.getSink()) == v)).forEachOrdered((e) -> {
                adj_type_Map.get(c.getKey()).remove(e);
            });
        });

        adj.remove(v);
        c_adj.remove(v);
        d_adj.remove(v);

        return this.vertices.remove(v);

    }

    /**
     * Return Edges of Type EdgeType (C-Edges/D-Edges etc.)
     *
     * @param e
     * @return Set of edges of type C-Edges/D-Edges etc
     */
    private Set<Edge> getEdgesOfType(EdgeType e) {
        Set<Edge> edgesOfTypes = new HashSet<>();
        if (adj_type_Map.containsKey(e)) {
            edgesOfTypes.addAll(adj_type_Map.get(e));
        }
        return edgesOfTypes;
    }

    public Set<Vertex> getCEdgeNeighbours(Vertex u) {
        validateVertex(u);
        return c_adj.containsKey(u) ? new HashSet<>(c_adj.get(u)) : new HashSet<>();
    }

    @Override
    public Iterator<Vertex> iterator() {
        return vertices.iterator();
    }

    public Set<Edge> getCEdges() {
        return getEdgesOfType(EdgeType.C_EDGE);
    }

    public Set<Edge> getDEdges() {
        return getEdgesOfType(EdgeType.D_EDGE);
    }

    /**
     * Return Index of Length
     * @param index
     * @return
     */
    public Vertex resolveVertex(Integer index) {
        if (vertices.size() > index) {
            return vertices.get(index);
        }
        return null;
    }
}
