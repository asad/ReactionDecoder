package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinception.smsd.core.SearchEngine;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.BaseMapping.Algorithm;

/**
 * Internal mapping/search abstraction for ReactionDecoder.
 *
 * Keeps SMSD construction policy in one place so the rest of the codebase can
 * work with a stable mapping interface instead of scattering raw constructor
 * calls for {@code Isomorphism} and {@code Substructure}.
 */
public interface ReactionMappingEngine {

    BaseMapping findMcs(IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher) throws CDKException;

    BaseMapping findMcs(IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            SearchEngine.McsOptions mcsOptions) throws CDKException;

    BaseMapping findSubstructure(IAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches) throws CDKException;

    BaseMapping findSubstructure(IAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException;

    BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches) throws CDKException;

    BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException;

    BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllMatches) throws CDKException;

    BaseMapping findSubstructure(IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllMatches,
            int maxMatches,
            long timeoutMs) throws CDKException;

    void applyDefaultFilters(BaseMapping mapping);
}
