package org.openscience.smsd.algorithm.vflib.vf2.sub;

import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;

/**
 * A structural pattern for finding an exact matching in a target compound.
 *
 * @author John May
 * @author Syed Asad Rahman
 * 
 */
public abstract class Pattern {

    /**
     * Patterns Algorithm Patterns.
     */
    public static enum Patterns {
        IDENTICAL, SUBGRAPH
    };

    /**
     * Determine if there is a mapping of this pattern in the {@code target}.
     * Depending on the implementation stereochemistry may be checked
     * (recommended).
     *
     * <blockquote><pre>
     * Pattern        pattern = ...; // create pattern
     * for (IAtomContainer m : ms) {
     *     if (pattern.matches(m)) {
     *         // found mapping!
     *     }
     * }
     * </pre></blockquote>
     *
     * @param target the container to search for the pattern in
     * @return the mapping from the pattern to the target
     */
    public final boolean matches(IAtomContainer target) {
        return !matchAll(target).isEmpty();
    }

    /**
     * Find all mappings of this pattern in the {@code target}.
     *
     * <blockquote><pre>
     * Pattern pattern = Pattern.findSubstructure(query);
     * for (IAtomContainer m : ms) {
     *     for (int[] mapping : pattern.matchAll(m)) {
     *         // found mapping
     *     }
     * }
     * </pre></blockquote>
     *
     * Using the fluent interface (see {@link AtomAtomMapping}) we can search
     * and manipulate the mappings. Here's an example of finding the first 5
     * mappings and creating an array. If the mapper is lazy other states are
     * simply not explored.
     *
     * <blockquote><pre>
     * // find only the first 5 mappings and store them in an array
     * Pattern pattern  = Pattern.findSubstructure(query);
     * List <Map<IAtom, IAtom>> mappings = pattern.matchAll(target)
     *                           .limit(5)
     *                           .toArray();
     * </pre></blockquote>
     *
     * @param target the container to search for the pattern in
     * @return the mapping from the pattern to the target
     * @see AtomAtomMapping
     */
    public abstract List<Map<IAtom, IAtom>> matchAll(IAtomContainer target);

    /**
     * Create a pattern which can be used to find molecules which contain the
     * {@code query} structure. The default structure search implementation is
     * {@link VF}.
     *
     * @param query the substructure to find
     * @return a pattern for finding the {@code query}
     * @see VF
     */
    public static Pattern findSubstructure(IQueryAtomContainer query) {
        return VF.findSubstructure(query);
    }

    /**
     * Create a pattern which can be used to find molecules which are the same
     * as the {@code query} structure. The default structure search
     * implementation is {@link VF}.
     *
     * @param query the substructure to find
     * @return a pattern for finding the {@code query}
     * @see VF
     */
    public static Pattern findIdentical(IQueryAtomContainer query) {
        return VF.findSubstructure(query);
    }

    /**
     * Create a pattern which can be used to find molecules which contain the
     * {@code query} structure. The default structure search implementation is
     * {@link VF}.
     *
     * @param query the substructure to find
     * @param shouldMatchBonds match bonds
     * @param shouldMatchRings match rings
     * @param matchAtomType match Atom types
     * @return a pattern for finding the {@code query}
     * @see VF
     */
    public static Pattern findSubstructure(IAtomContainer query, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        return VF.findSubstructure(query, shouldMatchBonds, shouldMatchRings, matchAtomType);
    }

    /**
     * Create a pattern which can be used to find molecules which are the same
     * as the {@code query} structure. The default structure search
     * implementation is {@link VF}.
     *
     * @param query the substructure to find
     * @param shouldMatchBonds match bonds
     * @param shouldMatchRings match rings
     * @param matchAtomType match Atom types
     * @return a pattern for finding the {@code query}
     * @see VF
     */
    public static Pattern findIdentical(IAtomContainer query, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        return VF.findSubstructure(query, shouldMatchBonds, shouldMatchRings, matchAtomType);
    }
}
