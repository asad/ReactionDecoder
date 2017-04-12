package org.openscience.smsd.algorithm.vflib.vf2.sub;

/**
 * Vento-Foggia (VF) state for matching subgraph-monomorphisms. The original
 * algorithm {
 *
 * @cdk.cite Cordella04} is for matching vertex-induced subgraph-isomorphisms. A
 * matching is vertex-induced if adjacency relation of mapped vertices is
 * isomorphic. Under such a procedure propane ({@code CCC}) is not considered a
 * substructure of cyclopropane ({@code C1CC1}). The term subgraph-isomorphism
 * is often conflated and it is really subgraph-monomorphism that is required
 * for substructure matching.
 * <p/>
 *
 * Note: no heuristics or sorting are performed at all and should be checked
 * externally
 *
 * @author John May
 * @author Syed Asad Rahman
 * 
 */
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.smsd.algorithm.vflib.vf2.AtomMatcher;
import org.openscience.smsd.algorithm.vflib.vf2.BondMatcher;

import static org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;

final class VFSubState extends AbstractVFState {

    /**
     * The query (container1) and target (container2) of the subgraph matching.
     */
    private final IAtomContainer queryContainer, targetContainer;

    /**
     * Lookup for the query bonds (bonds1) and target bonds (bonds2) of the
     * subgraph matching.
     */
    private final EdgeToBondMap queryBond, targetBond;

    /**
     * Defines how atoms are matched.
     */
    private final AtomMatcher atomMatcher;

    /**
     * Defines how bonds are matched.
     */
    private final BondMatcher bondMatcher;

    /**
     * Create a VF state for matching subgraph-monomorphism. The query is passed
     * first and should read as, find queryContainer in targetContainer.
     *
     * @param container1 the molecule to search for (query)
     * @param container2 the molecule to search in (target)
     * @param g1 adjacency list of the query
     * @param g2 adjacency list of the target
     * @param bonds1 bond lookup of the query
     * @param bonds2 bond lookup of the target
     * @param atomMatcher what semantic attributes (symbol, charge, query)
     * determines atoms to be compatible
     * @param bondMatcher what semantic attributes (order/aromatic, query)
     * determines bonds to be compatible
     */
    VFSubState(IAtomContainer container1, IAtomContainer container2, int[][] g1, int[][] g2, EdgeToBondMap bonds1,
            EdgeToBondMap bonds2, AtomMatcher atomMatcher, BondMatcher bondMatcher) {
        super(g1, g2);
        this.queryContainer = container1;
        this.targetContainer = container2;
        this.queryBond = bonds1;
        this.targetBond = bonds2;
        this.atomMatcher = atomMatcher;
        this.bondMatcher = bondMatcher;
    }

    /**
     * Check the feasibility of the candidate pair {n, m}. A candidate pair is
     * syntactically isMatchFeasible iff all k-look-ahead rules hold. These look
     * ahead rules check adjacency relation of the mapping. If an edge is mapped
     * in g1 it should also be mapped in g2 and vise-versa (0-look-ahead). If an
     * edge in g1 is unmapped but the edge is adjacent to an another mapped
     * vertex (terminal) then the number of such edges should be less or equal
     * in g1 compared to g2 (1-look-ahead). If the edge is unmapped and
     * non-terminal then the number of such edges should be less or equal in g1
     * compared to g2 (2-look-ahead).
     * <p/>
     * The above feasibility rules are for subgraph-isomorphism and have been
     * adapted for subgraph-monomorphism. For a monomorphism a mapped edge in g2
     * does not have to be present in g1. The 2-look-ahead also requires summing
     * the terminal and remaining counts (or sorting the vertices).
     * <p/>
     * The semantic feasibility verifies that the labels the label n, m are
     * compatible and that the label on each matched edge is compatible.
     *
     * @param n a candidate vertex from g1
     * @param m a candidate vertex from g2
     * @return the mapping is isMatchFeasible
     */
    @Override
    boolean isMatchFeasible(int n, int m) {

        // verify atom semantic feasibility
        if (!atomMatcher.matches(queryContainer.getAtom(n), targetContainer.getAtom(m))) {
            return false;
        }

        // unmapped terminal vertices n and m are adjacent to
        int nTerminal1 = 0, nTerminal2 = 0;
        // unmapped non-terminal (remaining) vertices n and m are adjacent to
        int nRemain1 = 0, nRemain2 = 0;

        // 0-look-ahead: check each adjacent edge for being mapped, and count
        // terminal or remaining
        for (int n_prime : g1[n]) {
            int m_prime = m1[n_prime];

            // v is already mapped, there should be an edge {m, w} in g2.
            if (m_prime != UNMAPPED) {
                IBond bond2 = targetBond.get(m, m_prime);
                if (bond2 == null) // the bond is not present in the target
                {
                    return false;
                }
                // verify bond semantic feasibility
                if (!bondMatcher.matches(queryBond.get(n, n_prime), bond2)) {
                    return false;
                }
            } else if (t1[n_prime] > 0) {
                nTerminal1++;
            } else {
                nRemain1++;
            }
        }

        // monomorphism: each mapped edge in g2 doesn't need to be in g1 so
        // only the terminal and remaining edges are counted
        for (int m_prime : g2[m]) {
            if (m2[m_prime] == UNMAPPED) {
                if (t2[m_prime] > 0) {
                    nTerminal2++;
                } else {
                    nRemain2++;
                }
            }
        }

        // 1-look-ahead : the mapping {n, m} is isMatchFeasible iff the number of
        // terminal vertices (t1) adjacent to n is less than or equal to the
        // number of terminal vertices (t2) adjacent to m.
        //
        // 2-look-ahead: the mapping {n, m} is isMatchFeasible iff the number of
        // vertices adjacent to n that are neither in m1 or t1 is less than or
        // equal to the number of the number of vertices adjacent to m that
        // are neither in m2 or t2. To allow mapping of monomorphisms we addMapping the
        // number of adjacent terminal vertices.
        return nTerminal1 <= nTerminal2 && (nRemain1 + nTerminal1) <= (nRemain2 + nTerminal2);
    }

}
