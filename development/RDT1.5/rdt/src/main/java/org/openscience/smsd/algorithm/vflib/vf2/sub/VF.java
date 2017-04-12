package org.openscience.smsd.algorithm.vflib.vf2.sub;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.algorithm.vflib.vf2.AtomMatcher;
import org.openscience.smsd.algorithm.vflib.vf2.BondMatcher;
import org.openscience.smsd.algorithm.vflib.vf2.DefaultAtomMatcher;
import org.openscience.smsd.algorithm.vflib.vf2.DefaultBondMatcher;

/**
 * A structure pattern which utilises the Vento-Foggia (VF) algorithm {
 *
 * @cdk.cite Cordella04}.
 *
 * <p/>
 *
 * Find and count the number molecules which contain the query substructure.
 *
 * <blockquote><pre>
 * IAtomContainer query   = ...;
 * Pattern        pattern = VF.findSubstructure(query);
 *
 * int hits = 0;
 * for (IAtomContainer m : ms)
 *     if (pattern.matches(m))
 *         hits++;
 * </pre></blockquote>
 * <p/>
 *
 * Finding the matching to molecules which contain the query substructure. It is
 * more efficient to obtain the {@link #match} and check it's size rather than
 * test if it {@link #matches}. These methods automatically verify
 * stereochemistry.
 *
 * <blockquote><pre>
 * IAtomContainer query   = ...;
 * Pattern        pattern = VF.findSubstructure(query);
 *
 * int hits = 0;
 * for (IAtomContainer m : ms) {
 *     int[] match = pattern.match(m);
 *     if (match.length > 0)
 *         hits++;
 * }
 * </pre></blockquote>
 *
 * @author John May
 * @author Syed Asad Rahman
 * 
 */
public final class VF extends Pattern {

    /**
     * The query structure.
     */
    private final IAtomContainer query;

    /**
     * The query structure adjacency list.
     */
    private final int[][] g1;

    /**
     * The bonds of the query structure.
     */
    private final EdgeToBondMap bonds1;

    /**
     * The atom matcher to determine atom feasibility.
     */
    private final AtomMatcher atomMatcher;

    /**
     * The bond matcher to determine atom feasibility.
     */
    private final BondMatcher bondMatcher;

    /**
     * Search for a subgraph.
     */
    private final Patterns searchType;

    /**
     * Non-public constructor for-now the atom/bond semantics are fixed.
     *
     * @param query the query structure
     * @param atomMatcher how atoms should be matched
     * @param bondMatcher how bonds should be matched
     * @param substructure substructure search
     */
    private VF(IAtomContainer query, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType, Patterns searchType) {
        this.query = query;
        this.atomMatcher = new DefaultAtomMatcher(shouldMatchRings, matchAtomType);
        this.bondMatcher = new DefaultBondMatcher(shouldMatchBonds, shouldMatchRings, matchAtomType);
        this.bonds1 = EdgeToBondMap.withSpaceFor(query);
        this.g1 = GraphUtil.toAdjList(query, bonds1);
        this.searchType = searchType;
    }

    /**
     * @inheritDoc
     */
    @Override
    public List<Map<IAtom, IAtom>> matchAll(final IAtomContainer target) {
        EdgeToBondMap bonds2 = EdgeToBondMap.withSpaceFor(target);
        int[][] g2 = GraphUtil.toAdjList(target, bonds2);
        Iterable<int[]> iterable = new VFIterable(query, target, g1, g2, bonds1, bonds2, atomMatcher, bondMatcher, searchType);
        List<Map<IAtom, IAtom>> mappings = new ArrayList<>();
        for (int[] map : iterable) {
            Map<IAtom, IAtom> atomAtomMapping = new HashMap<>();
            for (int i = 0; i < map.length; i++) {
                if (map[i] < 0) {
                    continue;
                }
                atomAtomMapping.put(query.getAtom(i), target.getAtom(map[i]));
            }
            mappings.add(atomAtomMapping);
        }
        return mappings;
    }

    /**
     * Create a pattern which can be used to find molecules which contain the
     * {@code query} structure.
     *
     * @param query the substructure to find
     * @return a pattern for finding the {@code query}
     */
    public static Pattern findSubstructure(IQueryAtomContainer query) {
        boolean isQuery = query instanceof IQueryAtomContainer;
        return findSubstructure(query, false, false, false);
    }

    /**
     * Create a pattern which can be used to find molecules which are the same
     * as the {@code query} structure.
     *
     * @param query the substructure to find
     * @return a pattern for finding the {@code query}
     */
    public static Pattern findIdentical(IQueryAtomContainer query) {
        boolean isQuery = query instanceof IQueryAtomContainer;
        return findIdentical(query, false, false, false);
    }

    /**
     * Create a pattern which can be used to find molecules which contain the
     * {@code query} structure.
     *
     * @param query the substructure to find
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return a pattern for finding the {@code query}
     */
    public static Pattern findSubstructure(IAtomContainer query, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        return new VF(query, shouldMatchBonds, shouldMatchRings, matchAtomType, Patterns.SUBGRAPH);
    }

    /**
     * Create a pattern which can be used to find molecules which are the same
     * as the {@code query} structure.
     *
     * @param query the substructure to find
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return a pattern for finding the {@code query}
     */
    public static Pattern findIdentical(IAtomContainer query, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        return new VF(query, shouldMatchBonds, shouldMatchRings, matchAtomType, Patterns.IDENTICAL);
    }

    private static final class VFIterable implements Iterable<int[]> {

        /**
         * Query and target containers.
         */
        private final IAtomContainer container1, container2;

        /**
         * Query and target adjacency lists.
         */
        private final int[][] g1, g2;

        /**
         * Query and target bond lookup.
         */
        private final EdgeToBondMap bonds1, bonds2;

        /**
         * How are atoms are matched.
         */
        private final AtomMatcher atomMatcher;

        /**
         * How are bonds are match.
         */
        private final BondMatcher bondMatcher;

        /**
         * The query is a subgraph.
         */
        private final Patterns searchType;

        /**
         * Create a match for the following parameters.
         *
         * @param container1 query structure
         * @param container2 target structure
         * @param g1 query adjacency list
         * @param g2 target adjacency list
         * @param bonds1 query bond map
         * @param bonds2 target bond map
         * @param atomMatcher how atoms are matched
         * @param bondMatcher how bonds are matched
         * @param subgraph perform subgraph search
         */
        private VFIterable(IAtomContainer container1, IAtomContainer container2, int[][] g1, int[][] g2,
                EdgeToBondMap bonds1, EdgeToBondMap bonds2, AtomMatcher atomMatcher, BondMatcher bondMatcher,
                Patterns searchType) {
            this.container1 = container1;
            this.container2 = container2;
            this.g1 = g1;
            this.g2 = g2;
            this.bonds1 = bonds1;
            this.bonds2 = bonds2;
            this.atomMatcher = atomMatcher;
            this.bondMatcher = bondMatcher;
            this.searchType = searchType;
        }

        /**
         * @inheritDoc
         */
        @Override
        public Iterator<int[]> iterator() {
            if (null != searchType) {
                switch (searchType) {
                    case IDENTICAL:
                        return new StateStream(
                                new VFState(container1, container2, g1, g2, bonds1, bonds2, atomMatcher, bondMatcher), searchType);
                    default:
                        return new StateStream(new VFSubState(container1, container2, g1, g2, bonds1, bonds2, atomMatcher, bondMatcher), searchType);
                }
            }
            return null;
        }
    }
}
