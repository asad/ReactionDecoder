package org.openscience.smsd.algorithm.vflib.vf2.mcs;

/**
 * Defines a state for matching (subgraph-)isomorphism from a query graph
 * (<i>G1</i>) to a target graph (<i>G2</i>). The mutable state allows
 * generation and adding and removal of mappings. A mapping {n, m} indicates a
 * query vertex (from <i>G1</i>), n, is paired (mapped) with the target vertex,
 * m (from <i>G2</i>). Candidate pairs are generated using {@link #hasNextCandidate(int)}
 * and {@link #nextM(int)}. Each candidate pair {n, m} is then {@link #addMapping}ed if
 * the mapping was feasible.
 *
 * @author John May
 * @author Syed Asad Rahman
 * 
 */
abstract class State {

    /**
     * Given the previous candidate generate the next query candidate. The first
     * candidate passed is always -1.
     *
     * @param n the previous candidate
     * @return next candidate
     */
    abstract int hasNextCandidate(int n);

    /**
     * Given the previous candidate generate the next target candidate. The
     * first candidate passed is always -1.
     *
     * @param n the current n vertex
     * @param m the previous candidate
     * @return next candidate
     */
    abstract int nextCandidate(int n, int m);

    /**
     * The max query candidate (number of vertices in the query).
     *
     * @return <i>|V| ∈ G1</i>
     */
    abstract int maxQueryCandidate();

    /**
     * The max target candidate (number of vertices in the target).
     *
     * @return <i>|V| ∈ G2</i>
     */
    abstract int maxTargetCandidate();

    /**
     * Add a mapping between n (a vertex G1) and m (a vertex in G2). If the
     * mapping was not feasible the mapping is not added.
     *
     * @param n a vertex in G1
     * @param m a vertex in G2
     * @return the mapping was added
     */
    abstract boolean addMapping(int n, int m);

    /**
     * Remove a mapping (backtrack) between n (a vertex G1) and m (a vertex in
     * G2).
     *
     * @param n a vertex in G1
     * @param m a vertex in G2
     */
    abstract void backTrack(int n, int m);

    /**
     * Access a copy of the current mapping.
     *
     * @return mapping of vertices from <i>G1</i> to <i>G2</i>
     */
    abstract int[] mapping();

    /**
     * Current size of the state. If <i>size</i> is the current number of mapped
     * candidates.
     *
     * @return the size of the state
     */
    abstract int size();
}
