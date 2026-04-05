/* Copyright (C) 2009-2020  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SearchEngine;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import static java.util.Collections.sort;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.ChemicalFilters;
/**
 *
 * java1.8+
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 */
public class BaseMapping extends ChemicalFilters implements ChemicalFilters.IAtomMapping {

    private boolean subgraph;
    private List<Double> stereoScoreList;
    private List<Integer> fragmentSizeList;
    private List<Double> bondEnergiesList;
    private final static ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(BaseMapping.class);
    final AtomMatcher atomMatcher;
    final BondMatcher bondMatcher;

    /**
     * Translate the legacy AtomMatcher/BondMatcher selection into SMSD 6.9.0
     * chemistry options so old call sites keep their historical semantics.
     */
    protected ChemOptions buildChemOptions() {
        ChemOptions options = new ChemOptions();
        configureAtomMatcher(options, atomMatcher);
        configureBondMatcher(options, bondMatcher);
        return options;
    }

    /**
     * Legacy algorithm enums are kept for source compatibility; the SMSD 6.9.0
     * engine is configured through McsOptions. Old constructors still get a
     * stable default and new callers can override specific flags.
     */
    protected SearchEngine.McsOptions buildMcsOptions(Algorithm algorithmType,
            SearchEngine.McsOptions overrides) {
        SearchEngine.McsOptions options = new SearchEngine.McsOptions();
        if (overrides != null) {
            options.induced = overrides.induced;
            options.connectedOnly = overrides.connectedOnly;
            options.timeoutMs = overrides.timeoutMs;
            options.extraSeeds = overrides.extraSeeds;
            options.seedNeighborhoodRadius = overrides.seedNeighborhoodRadius;
            options.seedMaxAnchors = overrides.seedMaxAnchors;
            options.useTwoHopNLFInExtension = overrides.useTwoHopNLFInExtension;
            options.useThreeHopNLFInExtension = overrides.useThreeHopNLFInExtension;
            options.disconnectedMCS = overrides.disconnectedMCS;
            options.maximizeBonds = overrides.maximizeBonds;
            options.minFragmentSize = overrides.minFragmentSize;
            options.maxFragments = overrides.maxFragments;
            options.atomWeights = overrides.atomWeights != null
                    ? Arrays.copyOf(overrides.atomWeights, overrides.atomWeights.length) : null;
            options.templateFuzzyAtoms = overrides.templateFuzzyAtoms;
            options.reactionAware = overrides.reactionAware;
            options.nearMcsDelta = overrides.nearMcsDelta;
            options.nearMcsCandidates = overrides.nearMcsCandidates;
            options.postFilter = overrides.postFilter;
            options.bondChangeAware = overrides.bondChangeAware;
            options.excludedTargetAtoms = overrides.excludedTargetAtoms != null
                    ? new LinkedHashSet<>(overrides.excludedTargetAtoms) : null;
        }
        if (options.timeoutMs <= 0) {
            options.timeoutMs = 10_000L;
        }
        return options;
    }

    private static void configureAtomMatcher(ChemOptions options, AtomMatcher matcher) {
        if (matcher == null) {
            return;
        }

        switch (matcher.toString()) {
            case "AnyMatcher":
                options.matchAtomType = false;
                options.matchFormalCharge = false;
                options.matchIsotope = false;
                options.ringMatchesRingOnly = false;
                break;
            case "RingElementMatcher":
            case "RingAtomTypeMatcher":
                options.matchAtomType = true;
                options.matchFormalCharge = true;
                options.matchIsotope = false;
                options.ringMatchesRingOnly = true;
                break;
            case "AtomTypeElementMatcher":
                options.matchAtomType = true;
                options.matchFormalCharge = true;
                options.matchIsotope = false;
                options.ringMatchesRingOnly = false;
                break;
            case "QueryMatcher":
                options.matchAtomType = true;
                options.matchFormalCharge = true;
                options.matchIsotope = false;
                options.ringMatchesRingOnly = false;
                break;
            case "ElementMatcher":
            default:
                options.matchAtomType = true;
                options.matchFormalCharge = true;
                options.matchIsotope = true;
                options.ringMatchesRingOnly = false;
                break;
        }
    }

    private static void configureBondMatcher(ChemOptions options, BondMatcher matcher) {
        if (matcher == null) {
            return;
        }

        switch (matcher.toString()) {
            case "AnyMatcher":
                options.matchBondOrder = ChemOptions.BondOrderMode.ANY;
                options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
                break;
            case "RingMatcher":
                options.matchBondOrder = ChemOptions.BondOrderMode.ANY;
                options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
                break;
            case "StrictOrderMatcher":
                options.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
                options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
                break;
            case "QueryMatcher":
                options.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
                options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
                break;
            case "OrderMatcher":
            default:
                options.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
                options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
                break;
        }
    }

    /**
     * Normalize a defensive clone before search so stricter SMSD versions do not
     * abort on legacy aromatic flag inconsistencies.
     */
    protected IAtomContainer normalizeForSearch(IAtomContainer container) throws CDKException {
        if (container == null || container instanceof IQueryAtomContainer) {
            return container;
        }
        try {
            IAtomContainer copy = container.clone();
            copy.setID(container.getID());
            copy.setProperties(container.getProperties());

            for (int i = 0; i < copy.getAtomCount(); i++) {
                IAtom sourceAtom = container.getAtom(i);
                if (sourceAtom != null) {
                    copy.getAtom(i).setID(sourceAtom.getID());
                }
            }

            for (IBond bond : copy.bonds()) {
                if (bond.getOrder() == null || bond.getOrder() == IBond.Order.UNSET) {
                    bond.setOrder(IBond.Order.SINGLE);
                }
                if (bond.isAromatic()) {
                    if (bond.getBegin() != null) {
                        bond.getBegin().setIsAromatic(true);
                    }
                    if (bond.getEnd() != null) {
                        bond.getEnd().setIsAromatic(true);
                    }
                }
            }
            return copy;
        } catch (CloneNotSupportedException ex) {
            throw new CDKException("Failed to normalize search molecule", ex);
        }
    }

    /**
     *
     * @param bm bond matcher
     * @param am atom matcher
     * @param mol1
     * @param mol2
     */
    public BaseMapping(IAtomContainer mol1, IAtomContainer mol2, AtomMatcher am, BondMatcher bm) {
        super(mol1, mol2);
        this.atomMatcher = am;
        this.bondMatcher = bm;
    }

    /**
     * @param mol1
     * @param mol2
     * @param am atom matcher
     * @param bm bond matcher
     */
    public BaseMapping(IQueryAtomContainer mol1, IAtomContainer mol2,
            AtomMatcher am, BondMatcher bm) {
        super(mol1, mol2);
        this.atomMatcher = am;
        this.bondMatcher = bm;

    }

    @Override
    public void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (getMappingCount() > 0) {
            this.fragmentSizeList = null;
            this.stereoScoreList = null;
            this.bondEnergiesList = null;

            if (fragmentFilter) {
                try {
                    sortResultsByFragments();
                    this.fragmentSizeList = getSortedFragment();
                } catch (RuntimeException ex) {
                    LOGGER.error(Level.SEVERE, "Fragment filter failed", ex);
                }
            }

            if (stereoFilter) {
                try {
                    sortResultsByStereoAndBondMatch();
                    this.stereoScoreList = getStereoMatches();
                } catch (CDKException | RuntimeException ex) {
                    LOGGER.error(Level.SEVERE, "Stereo filter failed", ex);
                }
            }

            if (energyFilter) {
                try {
                    sortResultsByEnergies();
                    this.bondEnergiesList = getSortedEnergy();
                } catch (CDKException | RuntimeException ex) {
                    LOGGER.error(Level.SEVERE, "Energy filter failed", ex);
                }
            }

            applyDeterministicMappingOrder();
        }
    }

    /**
     * Equivalent SMSD mappings often differ only in how symmetric atoms are
     * paired. Reorder those ties deterministically so callers get a stable
     * "first" mapping and score index 0 stays aligned with that choice.
     */
    private void applyDeterministicMappingOrder() {
        List<AtomAtomMapping> mappings = getMCSList();
        if (mappings.size() < 2) {
            return;
        }

        Map<AtomAtomMapping, MappingSortKey> sortKeys = new IdentityHashMap<>();
        for (AtomAtomMapping mapping : mappings) {
            sortKeys.put(mapping, buildMappingSortKey(mapping));
        }

        List<Integer> order = new ArrayList<>(mappings.size());
        for (int i = 0; i < mappings.size(); i++) {
            order.add(i);
        }
        order.sort((leftIndex, rightIndex)
                -> compareMappings(
                        sortKeys.get(mappings.get(leftIndex)),
                        sortKeys.get(mappings.get(rightIndex))));

        boolean changed = false;
        for (int i = 0; i < order.size(); i++) {
            if (order.get(i) != i) {
                changed = true;
                break;
            }
        }
        if (!changed) {
            return;
        }

        List<AtomAtomMapping> reorderedMappings = new ArrayList<>(mappings.size());
        for (Integer index : order) {
            reorderedMappings.add(mappings.get(index));
        }
        mappings.clear();
        mappings.addAll(reorderedMappings);

        fragmentSizeList = reorderScores(fragmentSizeList, order);
        stereoScoreList = reorderScores(stereoScoreList, order);
        bondEnergiesList = reorderScores(bondEnergiesList, order);
    }

    private int compareMappings(MappingSortKey left, MappingSortKey right) {
        int byMappedAtoms = Integer.compare(right.mappedAtoms, left.mappedAtoms);
        if (byMappedAtoms != 0) {
            return byMappedAtoms;
        }

        int byMappedBonds = Integer.compare(right.mappedBonds, left.mappedBonds);
        if (byMappedBonds != 0) {
            return byMappedBonds;
        }

        int byPairDistance = Integer.compare(left.pairDistanceScore, right.pairDistanceScore);
        if (byPairDistance != 0) {
            return byPairDistance;
        }

        Iterator<Map.Entry<IAtom, IAtom>> leftIterator = left.sortedPairs.iterator();
        Iterator<Map.Entry<IAtom, IAtom>> rightIterator = right.sortedPairs.iterator();

        while (leftIterator.hasNext() && rightIterator.hasNext()) {
            Map.Entry<IAtom, IAtom> leftPair = leftIterator.next();
            Map.Entry<IAtom, IAtom> rightPair = rightIterator.next();

            int byQueryKey = compareAtomOrder(leftPair.getKey(), rightPair.getKey(), left.query);
            if (byQueryKey != 0) {
                return byQueryKey;
            }

            int byTargetKey = compareAtomOrder(leftPair.getValue(), rightPair.getValue(), left.target);
            if (byTargetKey != 0) {
                return byTargetKey;
            }
        }
        return Integer.compare(left.sortedPairs.size(), right.sortedPairs.size());
    }

    private MappingSortKey buildMappingSortKey(AtomAtomMapping mapping) {
        List<Map.Entry<IAtom, IAtom>> sortedPairs = getSortedPairs(mapping);
        return new MappingSortKey(
                mapping.getQuery(),
                mapping.getTarget(),
                mapping.getCount(),
                makeBondMapOfAtomMap(mapping.getQuery(), mapping.getTarget(), mapping).size(),
                calculatePairDistanceScore(sortedPairs, mapping.getQuery(), mapping.getTarget()),
                sortedPairs);
    }

    private List<Map.Entry<IAtom, IAtom>> getSortedPairs(AtomAtomMapping mapping) {
        List<Map.Entry<IAtom, IAtom>> pairs = new ArrayList<>(mapping.getMappingsByAtoms().entrySet());
        pairs.sort((left, right) -> {
            int byQueryKey = compareAtomOrder(left.getKey(), right.getKey(), mapping.getQuery());
            if (byQueryKey != 0) {
                return byQueryKey;
            }
            return compareAtomOrder(left.getValue(), right.getValue(), mapping.getTarget());
        });
        return pairs;
    }

    private int compareAtomOrder(IAtom left, IAtom right, IAtomContainer container) {
        return buildAtomOrderKey(left, container).compareTo(buildAtomOrderKey(right, container));
    }

    private AtomOrderKey buildAtomOrderKey(IAtom atom, IAtomContainer container) {
        if (atom == null) {
            return new AtomOrderKey(9, Integer.MAX_VALUE, Integer.MAX_VALUE,
                    Integer.MAX_VALUE, Integer.MAX_VALUE, "");
        }

        AtomOrderKey benchmarkOrder = buildStableAtomOrderKey(atom.getProperty("benchmarkAtomId"), 0);
        if (benchmarkOrder != null) {
            return benchmarkOrder;
        }

        AtomOrderKey sourceOrder = buildStableAtomOrderKey(atom.getProperty("sourceAtomId"), 1);
        if (sourceOrder != null) {
            return sourceOrder;
        }

        Integer oldRank = integerProperty(atom.getProperty("OLD_RANK"));
        if (oldRank != null) {
            return new AtomOrderKey(2, oldRank, 0, 0, 0, "");
        }

        Integer originalIndex = integerProperty(atom.getProperty("index"));
        if (originalIndex != null) {
            return new AtomOrderKey(3, originalIndex, 0, 0, 0, "");
        }

        Integer label = integerProperty(atom.getProperty("label"));
        if (label != null) {
            return new AtomOrderKey(4, label, 0, 0, 0, "");
        }

        if (atom.getMapIdx() > 0) {
            return new AtomOrderKey(5, atom.getMapIdx(), 0, 0, 0, "");
        }

        Integer numericId = integerProperty(atom.getID());
        if (numericId != null) {
            return new AtomOrderKey(6, numericId, 0, 0, 0, "");
        }
        if (atom.getID() != null) {
            return new AtomOrderKey(7, 0, 0, 0, 0, atom.getID());
        }

        int fallbackIndex = container != null ? container.indexOf(atom) : Integer.MAX_VALUE;
        return new AtomOrderKey(8, fallbackIndex, 0, 0, 0, "");
    }

    private AtomOrderKey buildStableAtomOrderKey(Object propertyValue, int priority) {
        if (propertyValue == null) {
            return null;
        }

        String text = propertyValue.toString();
        String[] parts = text.split(":");
        if (parts.length != 3) {
            return null;
        }

        Integer moleculeIndex = tryParseInt(parts[1]);
        Integer atomIndex = tryParseInt(parts[2]);
        if (moleculeIndex == null || atomIndex == null) {
            return null;
        }

        int side = "P".equalsIgnoreCase(parts[0]) ? 1 : 0;
        return new AtomOrderKey(priority, side, moleculeIndex, atomIndex, 0, text);
    }

    private int calculatePairDistanceScore(List<Map.Entry<IAtom, IAtom>> pairs,
            IAtomContainer query, IAtomContainer target) {
        int total = 0;
        for (Map.Entry<IAtom, IAtom> pair : pairs) {
            total += Math.abs(
                    stableSequenceIndex(pair.getKey(), query)
                    - stableSequenceIndex(pair.getValue(), target));
        }
        return total;
    }

    private int stableSequenceIndex(IAtom atom, IAtomContainer container) {
        if (atom == null) {
            return Integer.MAX_VALUE / 4;
        }

        Integer stableIndex = stableSequenceIndex(atom.getProperty("benchmarkAtomId"));
        if (stableIndex != null) {
            return stableIndex;
        }

        stableIndex = stableSequenceIndex(atom.getProperty("sourceAtomId"));
        if (stableIndex != null) {
            return stableIndex;
        }

        Integer oldRank = integerProperty(atom.getProperty("OLD_RANK"));
        if (oldRank != null) {
            return oldRank;
        }

        Integer originalIndex = integerProperty(atom.getProperty("index"));
        if (originalIndex != null) {
            return originalIndex;
        }

        Integer label = integerProperty(atom.getProperty("label"));
        if (label != null) {
            return label;
        }

        if (atom.getMapIdx() > 0) {
            return atom.getMapIdx();
        }

        Integer numericId = integerProperty(atom.getID());
        if (numericId != null) {
            return numericId;
        }

        return container != null ? container.indexOf(atom) : Integer.MAX_VALUE / 4;
    }

    private Integer stableSequenceIndex(Object propertyValue) {
        if (propertyValue == null) {
            return null;
        }

        String text = propertyValue.toString();
        String[] parts = text.split(":");
        if (parts.length != 3) {
            return null;
        }

        Integer moleculeIndex = tryParseInt(parts[1]);
        Integer atomIndex = tryParseInt(parts[2]);
        if (moleculeIndex == null || atomIndex == null) {
            return null;
        }
        return moleculeIndex * 10_000 + atomIndex;
    }

    private Integer integerProperty(Object value) {
        if (value instanceof Integer intVal) {
            return intVal;
        }
        if (value != null) {
            return tryParseInt(value.toString());
        }
        return null;
    }

    private Integer tryParseInt(String value) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException ex) {
            return null;
        }
    }

    private static final class AtomOrderKey implements Comparable<AtomOrderKey> {

        private final int priority;
        private final int primary;
        private final int secondary;
        private final int tertiary;
        private final int quaternary;
        private final String text;

        private AtomOrderKey(int priority, int primary, int secondary,
                int tertiary, int quaternary, String text) {
            this.priority = priority;
            this.primary = primary;
            this.secondary = secondary;
            this.tertiary = tertiary;
            this.quaternary = quaternary;
            this.text = text == null ? "" : text;
        }

        @Override
        public int compareTo(AtomOrderKey other) {
            int byPriority = Integer.compare(priority, other.priority);
            if (byPriority != 0) {
                return byPriority;
            }
            int byPrimary = Integer.compare(primary, other.primary);
            if (byPrimary != 0) {
                return byPrimary;
            }
            int bySecondary = Integer.compare(secondary, other.secondary);
            if (bySecondary != 0) {
                return bySecondary;
            }
            int byTertiary = Integer.compare(tertiary, other.tertiary);
            if (byTertiary != 0) {
                return byTertiary;
            }
            int byQuaternary = Integer.compare(quaternary, other.quaternary);
            if (byQuaternary != 0) {
                return byQuaternary;
            }
            return text.compareTo(other.text);
        }
    }

    private static final class MappingSortKey {

        private final IAtomContainer query;
        private final IAtomContainer target;
        private final int mappedAtoms;
        private final int mappedBonds;
        private final int pairDistanceScore;
        private final List<Map.Entry<IAtom, IAtom>> sortedPairs;

        private MappingSortKey(IAtomContainer query, IAtomContainer target,
                int mappedAtoms, int mappedBonds, int pairDistanceScore,
                List<Map.Entry<IAtom, IAtom>> sortedPairs) {
            this.query = query;
            this.target = target;
            this.mappedAtoms = mappedAtoms;
            this.mappedBonds = mappedBonds;
            this.pairDistanceScore = pairDistanceScore;
            this.sortedPairs = sortedPairs;
        }
    }

    private static <T> List<T> reorderScores(List<T> scores, List<Integer> order) {
        if (scores == null || scores.size() != order.size()) {
            return scores;
        }
        List<T> reordered = new ArrayList<>(scores.size());
        for (Integer index : order) {
            reordered.add(scores.get(index));
        }
        return reordered;
    }

    @Override
    public Integer getFragmentSize(int Key) {
        return (fragmentSizeList != null && !fragmentSizeList.isEmpty())
                ? fragmentSizeList.get(Key) : null;
    }

    @Override
    public Integer getStereoScore(int Key) {
        return (stereoScoreList != null && !stereoScoreList.isEmpty()) ? stereoScoreList.get(Key).intValue() : null;
    }

    @Override
    public Double getEnergyScore(int Key) {
        return (bondEnergiesList != null && !bondEnergiesList.isEmpty()) ? bondEnergiesList.get(Key) : null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTanimotoSimilarity() {
        int decimalPlaces = 4;
        double rAtomCount;
        double pAtomCount;
        double tanimotoAtom = 0.0;

        if (getMappingCount() > 0) {
            applyDeterministicMappingOrder();
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                rAtomCount = (double) this.getMCSList().iterator().next().getQuery().getAtomCount();
                pAtomCount = (double) this.getMCSList().iterator().next().getTarget().getAtomCount();

                double matchCount = (double) firstAtomMCS.getCount();
                tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
                BigDecimal tan = new BigDecimal(tanimotoAtom);
                tan = tan.setScale(decimalPlaces, RoundingMode.HALF_UP);
                tanimotoAtom = tan.doubleValue();
            }
        }
        return tanimotoAtom;
    }

    /**
     * {@inheritDoc}
     *
     */
    @Override
    public boolean isStereoMisMatch() {
        boolean flag = false;
        IAtomContainer reactant = getQuery();
        IAtomContainer product = getTarget();
        int stereoMisMatchScore = 0;
        if (getMappingCount() > 0) {
            applyDeterministicMappingOrder();
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();
            for (IAtom indexI : firstAtomMCS.getMappingsByAtoms().keySet()) {
                IAtom indexJ = firstAtomMCS.getMappingsByAtoms().get(indexI);
                for (IAtom indexIPlus : firstAtomMCS.getMappingsByAtoms().keySet()) {
                    IAtom indexJPlus = firstAtomMCS.getMappingsByAtoms().get(indexIPlus);
                    if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                        IAtom sourceAtom1 = indexI;
                        IAtom sourceAtom2 = indexIPlus;
                        IBond rBond = reactant.getBond(sourceAtom1, sourceAtom2);

                        IAtom targetAtom1 = indexJ;
                        IAtom targetAtom2 = indexJPlus;
                        IBond pBond = product.getBond(targetAtom1, targetAtom2);

                        if ((rBond != null && pBond != null)
                                && (rBond.getStereo() != pBond.getStereo())) {
                            stereoMisMatchScore++;
                        }
                    }
                }
            }
        }
        if (stereoMisMatchScore > 0) {
            flag = true;
        }
        return flag;
    }

    @Override
    public int getMappingCount() {
        return this.getMCSList().isEmpty() ? 0 : this.getMCSList().size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getEuclideanDistance() {
        int decimalPlaces = 4;
        double sourceAtomCount;
        double targetAtomCount;
        double euclidean = -1.;

        if (getMappingCount() > 0) {
            applyDeterministicMappingOrder();
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                sourceAtomCount = (double) this.getMCSList().iterator()
                        .next().getQuery().getAtomCount();
                targetAtomCount = (double) this.getMCSList().iterator()
                        .next().getTarget().getAtomCount();

                double common = (double) firstAtomMCS.getCount();
                euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
                BigDecimal dist = new BigDecimal(euclidean);
                dist = dist.setScale(decimalPlaces, RoundingMode.HALF_UP);
                euclidean = dist.doubleValue();
            }
        }
        return euclidean;
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public List<AtomAtomMapping> getAllAtomMapping() {
        applyDeterministicMappingOrder();
        return Collections.unmodifiableList(new ArrayList<>(getMCSList()));
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public AtomAtomMapping getFirstAtomMapping() {
        applyDeterministicMappingOrder();
        return getMCSList().isEmpty() ? new AtomAtomMapping(getQuery(), getTarget())
                : getMCSList().iterator().next();
    }

    /**
     * Returns true if Query is a subgraph of the Target.
     *
     * @return true if Query is a subgraph of the Target
     */
    public boolean isSubgraph() {
        return this.subgraph;
    }

    public void clearMaps() {
        this.getMCSList().clear();
    }

    /**
     * @return the allBondMCS
     */
    public List<Map<IBond, IBond>> getAllBondMaps() {
        if (!getMCSList().isEmpty()) {
            applyDeterministicMappingOrder();
            return makeBondMapsOfAtomMaps(getQuery(), getTarget(), getMCSList());
        }
        return new ArrayList<>();
    }

    /**
     * @param subgraph the subgraph to set
     */
    public void setSubgraph(boolean subgraph) {
        this.subgraph = subgraph;
    }

    /**
     * Returns bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mappings mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    public List<Map<IBond, IBond>> makeBondMapsOfAtomMaps(IAtomContainer ac1,
            IAtomContainer ac2, List<AtomAtomMapping> mappings) {
        List<Map<IBond, IBond>> bondMaps = new ArrayList<>();
        mappings.stream().forEach((mapping) -> {
            bondMaps.add(makeBondMapOfAtomMap(ac1, ac2, mapping));
        });
        return bondMaps;
    }

    /**
     *
     * Returns bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mapping mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    private Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2,
            AtomAtomMapping mapping) {

        Map<IBond, IBond> bondbondMappingMap = new HashMap<>();

        mapping.getMappingsByAtoms().entrySet().stream().forEach((Map.Entry<IAtom, IAtom> map1) -> {
            mapping.getMappingsByAtoms().entrySet().stream().filter((map2) -> (map1.getKey()
                    != map2.getKey())).forEach((Map.Entry<IAtom, IAtom> map2) -> {
                IBond bond1;
                bond1 = ac1.getBond(map1.getKey(), map2.getKey());
                IBond bond2 = ac2.getBond(map1.getValue(), map2.getValue());
                if (bond1 != null && bond2 != null && !bondbondMappingMap.containsKey(bond1)) {
                    bondbondMappingMap.put(bond1, bond2);
                }
            });
        });
        return bondbondMappingMap;
    }

    int expectedMaxGraphmatch(IAtomContainer q, IAtomContainer t) {

        /*
         a={c,c,c,o,n}
         b={c,c,c,p}
       
         expectedMaxGraphmatch=3;
         */
        List<String> atomUniqueCounter1 = new ArrayList<>();
        List<String> atomUniqueCounter2 = new ArrayList<>();

        for (IAtom a : q.atoms()) {
            String hyb = a.getAtomicNumber() == null
                    ? a.getSymbol() : a.getAtomicNumber() + "";
            atomUniqueCounter1.add(hyb);
        }

        for (IAtom b : t.atoms()) {
            String hyb = b.getAtomicNumber() == null
                    ? b.getSymbol() : b.getAtomicNumber() + "";
            atomUniqueCounter2.add(hyb);
        }

        sort(atomUniqueCounter1);
        sort(atomUniqueCounter2);

        if (atomUniqueCounter1.isEmpty()) {
            return 0;
        }
        List<String> common = new LinkedList<>(atomUniqueCounter1);
        common.retainAll(atomUniqueCounter2);

        atomUniqueCounter1.clear();
        atomUniqueCounter2.clear();
        return common.size();
    }

    // ==================== Inner enum Algorithm ====================

    /**
     * This class represents various algorithm type supported by SMSD.
     *
     * @author Syed Asad Rahman
     */
    public static enum Algorithm {

        /** Default SMSD algorithm. */
        DEFAULT(0, "Default SMSD algorithm"),
        /** MCS Plus algorithm. */
        MCSPlus(2, "MCS Plus algorithm"),
        /** VF-Koch-McGregor Lib based MCS algorithm. */
        VFLibMCS(3, "VF-Koch-McGregor Lib based MCS algorithm"),
        /** CDK UIT MCS. */
        CDKMCS(4, "CDK UIT MCS");

        private final int type;
        private final String description;

        Algorithm(int aStatus, String desc) {
            this.type = aStatus;
            this.description = desc;
        }

        /** Returns type of algorithm. */
        public int type() {
            return this.type;
        }

        /** Returns short description of the algorithm. */
        public String description() {
            return this.description;
        }
    }

}
