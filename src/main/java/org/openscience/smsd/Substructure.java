/*
 * Copyright (C) 2009-2020  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;

/**
 * Substructure search adapter that delegates to SMSD 3.0.0.
 *
 * Maintains the same public API as the original Substructure class
 * but uses the new SMSD 3.0.0 engine (VF2++) for substructure detection.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class Substructure extends BaseMapping {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Substructure.class);
    private int vfMappingSize = -1;
    private final int maxMatches;
    private final long timeoutMs;
    private SearchEngine.SubstructureResult lastSearchResult;

    /**
     * Constructor for VF Substructure Algorithm.
     *
     * @param query query molecule
     * @param target target molecule
     * @param am atom matcher
     * @param bm bond matcher
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IAtomContainer query,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm,
            boolean findAllSubgraph) throws CDKException {
        this(query, target, am, bm, findAllSubgraph, findAllSubgraph ? 10 : 1, 30_000L);
    }

    /**
     * Constructor for VF Substructure Algorithm with explicit enumeration flags.
     *
     * @param query query molecule
     * @param target target molecule
     * @param am atom matcher
     * @param bm bond matcher
     * @param findAllSubgraph report all subgraphs
     * @param maxMatches maximum number of mappings to enumerate
     * @param timeoutMs timeout in milliseconds
     * @throws CDKException
     */
    public Substructure(
            IAtomContainer query,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm,
            boolean findAllSubgraph,
            int maxMatches,
            long timeoutMs) throws CDKException {
        super(query, target, am, bm);
        this.maxMatches = Math.max(1, maxMatches);
        this.timeoutMs = timeoutMs > 0 ? timeoutMs : 30_000L;
        super.setSubgraph(findSubgraphs(findAllSubgraph));
    }

    /**
     * Constructor for IQueryAtomContainer.
     *
     * @param query query container
     * @param target target molecule
     * @param am atom matcher
     * @param bm bond matcher
     * @param findAllSubgraphFlag report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm,
            boolean findAllSubgraphFlag) throws CDKException {
        this(query, target, am, bm, findAllSubgraphFlag, findAllSubgraphFlag ? 10 : 1, 30_000L);
    }

    /**
     * Constructor for IQueryAtomContainer with explicit enumeration flags.
     *
     * @param query query container
     * @param target target molecule
     * @param am atom matcher
     * @param bm bond matcher
     * @param findAllSubgraphFlag report all subgraphs
     * @param maxMatches maximum number of mappings to enumerate
     * @param timeoutMs timeout in milliseconds
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm,
            boolean findAllSubgraphFlag,
            int maxMatches,
            long timeoutMs) throws CDKException {
        super(query, target, am, bm);
        this.maxMatches = Math.max(1, maxMatches);
        this.timeoutMs = timeoutMs > 0 ? timeoutMs : 30_000L;
        super.setSubgraph(findSubgraphs(findAllSubgraphFlag));
    }

    /**
     * Constructor for IQueryAtomContainer with default matchers.
     *
     * @param query query container
     * @param target target molecule
     * @param findAllSubgraphFlag report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllSubgraphFlag) throws CDKException {
        this(query, target, findAllSubgraphFlag, findAllSubgraphFlag ? 10 : 1, 30_000L);
    }

    /**
     * Constructor for IQueryAtomContainer with default matchers and explicit
     * subgraph enumeration flags.
     *
     * @param query query container
     * @param target target molecule
     * @param findAllSubgraphFlag report all subgraphs
     * @param maxMatches maximum number of mappings to enumerate
     * @param timeoutMs timeout in milliseconds
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllSubgraphFlag,
            int maxMatches,
            long timeoutMs) throws CDKException {
        super(query, target, AtomMatcher.forQuery(), BondMatcher.forQuery());
        this.maxMatches = Math.max(1, maxMatches);
        this.timeoutMs = timeoutMs > 0 ? timeoutMs : 30_000L;
        super.setSubgraph(findSubgraphs(findAllSubgraphFlag));
    }

    /**
     * Delegates substructure search to SMSD 3.0.0.
     */
    private boolean findSubgraphs(boolean findAllMatches) throws CDKException {
        if (getTarget() == null || getQuery() == null) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (getQuery().getAtomCount() > getTarget().getAtomCount()) {
            return false;
        }

        int rAtomCount = getQuery().getAtomCount();
        int pAtomCount = getTarget().getAtomCount();
        int expectedMaxGraphmatch = expectedMaxGraphmatch(getQuery(), getTarget());

        // Single atom case
        if (expectedMaxGraphmatch == 1 && rAtomCount <= pAtomCount) {
            return singleMapping();
        }

        try {
            IAtomContainer searchQuery = normalizeForSearch(getQuery());
            IAtomContainer searchTarget = normalizeForSearch(getTarget());
            ChemOptions chemOptions = buildChemOptions();
            int limit = findAllMatches ? maxMatches : 1;
            lastSearchResult = SearchEngine.findAllSubstructuresWithStats(
                    searchQuery, searchTarget, chemOptions, limit, timeoutMs);

            if (lastSearchResult != null
                    && lastSearchResult.exists()
                    && lastSearchResult.mappings() != null
                    && !lastSearchResult.mappings().isEmpty()) {
                for (Map<Integer, Integer> mapping : lastSearchResult.mappings()) {
                    AtomAtomMapping aam = convertMapping(getQuery(), getTarget(), mapping);
                    if (!aam.isEmpty() && aam.getCount() >= vfMappingSize) {
                        if (aam.getCount() > vfMappingSize) {
                            vfMappingSize = aam.getCount();
                            getMCSList().clear();
                        }
                        if (!hasMap(aam, getMCSList())) {
                            getMCSList().add(aam);
                        }
                    }
                }
                return !getMCSList().isEmpty();
            }
        } catch (Exception e) {
            LOGGER.error(Level.SEVERE, "Error in SMSD substructure search", e);
            throw new CDKException("Substructure search failed: " + e.getMessage(), e);
        }

        return false;
    }

    private AtomAtomMapping convertMapping(IAtomContainer query, IAtomContainer target,
            Map<Integer, Integer> indexMapping) {
        AtomAtomMapping aam = new AtomAtomMapping(query, target);
        for (Map.Entry<Integer, Integer> entry : indexMapping.entrySet()) {
            int qIdx = entry.getKey();
            int tIdx = entry.getValue();
            if (qIdx >= 0 && qIdx < query.getAtomCount()
                    && tIdx >= 0 && tIdx < target.getAtomCount()) {
                IAtom qAtom = query.getAtom(qIdx);
                IAtom tAtom = target.getAtom(tIdx);
                if (qAtom != null && tAtom != null) {
                    aam.put(qAtom, tAtom);
                }
            }
        }
        return aam;
    }

    private boolean hasMap(AtomAtomMapping map, List<AtomAtomMapping> mapGlobal) {
        return mapGlobal.stream().anyMatch(test -> test.equals(map));
    }

    public SearchEngine.SubstructureResult getLastSearchResult() {
        return lastSearchResult;
    }

    private boolean singleMapping() {
        IAtomContainer query = getQuery();
        IAtomContainer target = getTarget();
        for (IAtom qAtom : query.atoms()) {
            for (IAtom tAtom : target.atoms()) {
                if (AtomMatcher.matchSymbol(qAtom, tAtom)) {
                    AtomAtomMapping aam = new AtomAtomMapping(query, target);
                    aam.put(qAtom, tAtom);
                    getMCSList().add(aam);
                    return true;
                }
            }
        }
        return false;
    }
}
