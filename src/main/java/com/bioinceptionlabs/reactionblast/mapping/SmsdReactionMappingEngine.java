package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinception.smsd.core.SearchEngine;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.BaseMapping.Algorithm;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;

/**
 * Default internal mapping engine backed by SMSD.
 */
public final class SmsdReactionMappingEngine implements ReactionMappingEngine {

    private static final ReactionMappingEngine INSTANCE = new SmsdReactionMappingEngine();

    public static ReactionMappingEngine getInstance() {
        return INSTANCE;
    }

    private SmsdReactionMappingEngine() {
    }

    @Override
    public BaseMapping findMcs(IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher) throws CDKException {
        return new Isomorphism(query, target, algorithmType, atomMatcher, bondMatcher);
    }

    @Override
    public BaseMapping findMcs(IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            SearchEngine.McsOptions mcsOptions) throws CDKException {
        return new Isomorphism(query, target, algorithmType, atomMatcher, bondMatcher, mcsOptions);
    }

    @Override
    public BaseMapping findSubstructure(IAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches) throws CDKException {
        return new Substructure(query, target, atomMatcher, bondMatcher, findAllMatches);
    }

    @Override
    public BaseMapping findSubstructure(IAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException {
        return new Substructure(query, target, atomMatcher, bondMatcher,
                findAllMatches, maxMatches, timeoutMs);
    }

    @Override
    public BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches) throws CDKException {
        return new Substructure(query, target, atomMatcher, bondMatcher, findAllMatches);
    }

    @Override
    public BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException {
        return new Substructure(query, target, atomMatcher, bondMatcher,
                findAllMatches, maxMatches, timeoutMs);
    }

    @Override
    public BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllMatches) throws CDKException {
        return new Substructure(query, target, findAllMatches);
    }

    @Override
    public BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException {
        return new Substructure(query, target, findAllMatches, maxMatches, timeoutMs);
    }

    @Override
    public void applyDefaultFilters(BaseMapping mapping) {
        if (mapping != null) {
            try {
                mapping.setChemFilters(true, true, true);
            } catch (NullPointerException ex) {
                // SMSD 6.10.1 energy filter can NPE on certain molecule pairs
                // where bond-energy lookup returns null. Fall back to stereo+fragment only.
                try {
                    mapping.setChemFilters(true, true, false);
                } catch (Exception fallback) {
                    // last resort — no filters at all
                }
            }
        }
    }
}
