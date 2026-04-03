/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinception.smsd.core.MolGraph;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

/**
 * Utility methods for building structure-based molecule and pair cache keys.
 *
 * The keys intentionally ignore occurrence-specific atom IDs so the same
 * molecular structure can reuse one MCS template across stoichiometric copies.
 */
public final class MappingKeyUtil {

    private static final SmilesGenerator CANONICAL_SMILES
            = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.Stereo);

    private MappingKeyUtil() {
    }

    public static String computeStructureKey(IAtomContainer molecule) {
        if (molecule == null) {
            return "null";
        }
        try {
            return new MolGraph(molecule).toCanonicalSmiles();
        } catch (RuntimeException ex) {
            try {
                return CANONICAL_SMILES.create(molecule);
            } catch (CDKException | RuntimeException ignored) {
                return molecule.getAtomCount() + ":" + molecule.getBondCount();
            }
        }
    }

    public static String buildPairKey(IAtomContainer query, IAtomContainer target,
            String modeLabel,
            boolean atomType, boolean bondMatch,
            boolean ringMatch, boolean ringSizeMatch) {
        return buildPairKey(
                computeStructureKey(query),
                computeStructureKey(target),
                modeLabel,
                atomType, bondMatch, ringMatch, ringSizeMatch);
    }

    public static String buildPairKey(String queryStructureKey, String targetStructureKey,
            String modeLabel,
            boolean atomType, boolean bondMatch,
            boolean ringMatch, boolean ringSizeMatch) {
        StringBuilder key = new StringBuilder();
        key.append(queryStructureKey)
                .append(">>")
                .append(targetStructureKey)
                .append('|')
                .append(modeLabel)
                .append('|')
                .append(atomType ? '1' : '0')
                .append(bondMatch ? '1' : '0')
                .append(ringMatch ? '1' : '0')
                .append(ringSizeMatch ? '1' : '0');
        return key.toString();
    }
}
