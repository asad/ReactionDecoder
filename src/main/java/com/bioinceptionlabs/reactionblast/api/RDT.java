/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package com.bioinceptionlabs.reactionblast.api;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mechanism.BondChangeCalculator;
import com.bioinceptionlabs.reactionblast.mechanism.MappingSolution;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;

/**
 * Simple, clean public API for Reaction Decoder Tool.
 *
 * <pre>
 * // Map a reaction from SMILES
 * ReactionResult result = RDT.map("CC(=O)O.OCC>>CC(=O)OCC.O");
 *
 * // Check results
 * System.out.println(result.getBondChanges());    // [C-O, O-H, C=O, ...]
 * System.out.println(result.getMappedSmiles());    // mapped SMILES
 * System.out.println(result.getBondChangeCount()); // number of bond changes
 * </pre>
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class RDT {

    private RDT() {}

    /**
     * Map a reaction from SMILES and extract bond changes.
     *
     * @param reactionSmiles reaction SMILES (reactants>>products)
     * @return ReactionResult with mapping, bond changes, and fingerprints
     * @throws IllegalArgumentException if SMILES is invalid
     */
    public static ReactionResult map(String reactionSmiles) {
        return map(reactionSmiles, true, true);
    }

    /**
     * Map a reaction with control over stereo perception and ring handling.
     *
     * @param reactionSmiles reaction SMILES (reactants>>products)
     * @param generate2D perceive 2D stereo centers
     * @param complexMapping handle ring system mapping
     * @return ReactionResult with mapping, bond changes, and fingerprints
     * @throws IllegalArgumentException if SMILES is invalid
     */
    public static ReactionResult map(String reactionSmiles, boolean generate2D, boolean complexMapping) {
        if (reactionSmiles == null || !reactionSmiles.contains(">>")) {
            throw new IllegalArgumentException("Invalid reaction SMILES: must contain '>>'");
        }
        try {
            SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
            IReaction reaction = sp.parseReactionSmiles(reactionSmiles);
            reaction.setID("RDT_" + Integer.toHexString(reactionSmiles.hashCode()));

            ReactionMechanismTool rmt = new ReactionMechanismTool(
                    reaction, true, generate2D, false, complexMapping, true, new StandardizeReaction());

            return extractResult(rmt, reactionSmiles);
        } catch (Exception e) {
            throw new RuntimeException("Mapping failed for: " + reactionSmiles, e);
        }
    }

    /**
     * Compare two reactions for similarity based on bond change fingerprints.
     *
     * @param smiles1 first reaction SMILES
     * @param smiles2 second reaction SMILES
     * @return Tanimoto similarity (0.0 = no overlap, 1.0 = identical changes)
     */
    public static double compare(String smiles1, String smiles2) {
        ReactionResult r1 = map(smiles1);
        ReactionResult r2 = map(smiles2);
        return r1.similarity(r2);
    }

    private static ReactionResult extractResult(ReactionMechanismTool rmt, String inputSmiles) {
        MappingSolution solution = rmt.getSelectedSolution();
        if (solution == null) {
            return new ReactionResult(inputSmiles, null, 0, 0, 0,
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), "NONE");
        }

        BondChangeCalculator bcc = solution.getBondChangeCalculator();
        List<String> formedCleaved;
        List<String> orderChanges;
        List<String> stereoChanges;
        List<String> reactionCentre;
        try {
            formedCleaved = extractFeatures(bcc.getFormedCleavedWFingerprint());
            orderChanges = extractFeatures(bcc.getOrderChangesWFingerprint());
            stereoChanges = extractFeatures(bcc.getStereoChangesWFingerprint());
            reactionCentre = extractFeatures(bcc.getReactionCenterWFingerprint());
        } catch (Exception e) {
            formedCleaved = new ArrayList<>();
            orderChanges = new ArrayList<>();
            stereoChanges = new ArrayList<>();
            reactionCentre = new ArrayList<>();
        }

        String mappedSmiles = null;
        try {
            org.openscience.cdk.smiles.SmilesGenerator sg = new org.openscience.cdk.smiles.SmilesGenerator(
                    org.openscience.cdk.smiles.SmiFlavor.Stereo | org.openscience.cdk.smiles.SmiFlavor.AtomAtomMap);
            mappedSmiles = sg.create(bcc.getReaction());
        } catch (Exception ignored) {}

        String algorithm = solution.getAlgorithmID() != null
                ? solution.getAlgorithmID().name() : "UNKNOWN";

        return new ReactionResult(
                inputSmiles,
                mappedSmiles,
                formedCleaved.size(),
                orderChanges.size(),
                stereoChanges.size(),
                formedCleaved,
                orderChanges,
                stereoChanges,
                reactionCentre,
                algorithm);
    }

    private static List<String> extractFeatures(IPatternFingerprinter fp) {
        List<String> features = new ArrayList<>();
        if (fp != null) {
            for (var feature : fp.getFeatures()) {
                features.add(feature.getPattern() + ":" + (int) feature.getWeight());
            }
        }
        return features;
    }
}
