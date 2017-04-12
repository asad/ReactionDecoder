package org.openscience.smsd.algorithm.vflib.vf2.sub;

import org.openscience.smsd.algorithm.vflib.vf2.sub.State;
import java.util.Arrays;

/**
 * A state for the Vento-Foggia (VF) algorithm. The state allows adding and
 * removing of mappings as well as generating the new candidate mappings {@link
 * #hasNextCandidate(int)} and {@link #nextCandidate(int, int)}. The feasibility check is left for
 * subclasses to implement.
 *
 * @author John May
 * @author Syed Asad Rahman
 *  isomorphism
 */
abstract class AbstractVFState extends State {

    /**
     * Value indicates a vertex is unmapped.
     */
    protected static final int UNMAPPED = -1;

    /**
     * Adjacency list representation of the containers.
     */
    protected final int[][] g1, g2;

    /**
     * Mapping - m1 is the the mapping from g1 to g1, m2 is from g2 to g1.
     */
    protected final int[] m1, m2;

    /**
     * The (terminal) vertices which are adjacent to each mapped pair.
     */
    protected final int[] t1, t2;

    /**
     * Size of current solution - the number of vertices matched.
     */
    protected int size;

    /**
     * Create a state which will be used to match g1 in g2.
     *
     * @param g1 find this graph
     * @param g2 search this graph
     */
    public AbstractVFState(final int[][] g1, final int[][] g2) {
        this.g1 = g1;
        this.g2 = g2;
        this.m1 = new int[g1.length];
        this.m2 = new int[g2.length];
        this.t1 = new int[g1.length];
        this.t2 = new int[g2.length];
        size = 0;
        Arrays.fill(m1, UNMAPPED);
        Arrays.fill(m2, UNMAPPED);
    }

    /**
     * Given the current query candidate (n), find the next candidate. The next
     * candidate is the next vertex > n (in some ordering) that is unmapped and
     * is adjacent to a mapped vertex (terminal). If there is no such vertex
     * (disconnected) the next unmapped vertex is returned. If there are no more
     * candidates m == |V| of G1.
     *
     * @param n previous candidate n
     * @return the next value of n
     */
    @Override
    final int hasNextCandidate(int n) {
        if (size == 0) {
            return 0;
        }
        for (int i = n + 1; i < g1.length; i++) {
            if (m1[i] == UNMAPPED && t1[i] > 0) {
                return i;
            }
        }
        for (int i = n + 1; i < g1.length; i++) {
            if (m1[i] == UNMAPPED) {
                return i;
            }
        }
        return maxQueryCandidate();
    }

    /**
     * Given the current target candidate (m), find the next candidate. The next
     * candidate is the next vertex > m (in some ordering) that is unmapped and
     * is adjacent to a mapped vertex (terminal). If there is no such vertex
     * (disconnected) the next unmapped vertex is returned. If there are no more
     * candidates m == |V| of G2.
     *
     * @param m previous candidate m
     * @return the next value of m
     */
    @Override
    final int nextCandidate(int n, int m) {
        if (size == 0) {
            return m + 1;
        }
        // if the query vertex 'n' is in the terminal set (t1) then the
        // target vertex must be in the terminal set (t2)
        for (int i = m + 1; i < g2.length; i++) {
            if (m2[i] == UNMAPPED && (t1[n] == 0 || t2[i] > 0)) {
                return i;
            }
        }
        return maxTargetCandidate();
    }

    /**
     * @inheritDoc
     */
    @Override
    final int maxQueryCandidate() {
        return g1.length;
    }

    /**
     * @inheritDoc
     */
    @Override
    final int maxTargetCandidate() {
        return g2.length;
    }

    /**
     * @inheritDoc
     */
    @Override
    final boolean addMapping(int n, int m) {
        if (!isMatchFeasible(n, m)) {
            return false;
        }
        m1[n] = m;
        m2[m] = n;
        size = size + 1;
        for (int w : g1[n]) {
            if (t1[w] == 0) {
                t1[w] = size;
            }
        }
        for (int w : g2[m]) {
            if (t2[w] == 0) {
                t2[w] = size;
            }
        }
        return true;
    }

    /**
     * @inheritDoc
     */
    @Override
    final void backTrack(int n, int m) {
        m1[n] = m2[m] = UNMAPPED;
        size = size - 1;
        for (int w : g1[n]) {
            if (t1[w] > size) {
                t1[w] = 0;
            }
        }
        for (int w : g2[m]) {
            if (t2[w] > size) {
                t2[w] = 0;
            }
        }
    }

    /**
     * Is the candidate pair {n, m} isMatchFeasible. Verifies if the adding candidate
     * pair {n, m} to the state would lead to an invalid mapping.
     *
     * @param n query vertex
     * @param m target vertex
     * @return the mapping is isMatchFeasible
     */
    abstract boolean isMatchFeasible(int n, int m);

    /**
     * @inheritDoc
     */
    @Override
    int[] mapping() {
        return Arrays.copyOf(m1, m1.length);
    }

    /**
     * @inheritDoc
     */
    @Override
    int size() {
        return size;
    }
}
