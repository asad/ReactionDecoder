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

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Immutable result of a reaction mapping. Contains bond changes,
 * mapped SMILES, and fingerprint features as plain Java types
 * (no CDK/toolkit dependency).
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class ReactionResult {

    private final String inputSmiles;
    private final String mappedSmiles;
    private final int formedCleavedCount;
    private final int orderChangeCount;
    private final int stereoChangeCount;
    private final List<String> formedCleavedBonds;
    private final List<String> orderChangedBonds;
    private final List<String> stereoChangedBonds;
    private final List<String> reactionCentreFingerprint;
    private final String algorithmUsed;
    private final String reactionSignature;
    private final String canonicalHash;

    ReactionResult(String inputSmiles, String mappedSmiles,
                   int formedCleavedCount, int orderChangeCount, int stereoChangeCount,
                   List<String> formedCleavedBonds, List<String> orderChangedBonds,
                   List<String> stereoChangedBonds, List<String> reactionCentreFingerprint,
                   String algorithmUsed) {
        this.inputSmiles = inputSmiles;
        this.mappedSmiles = mappedSmiles;
        this.formedCleavedCount = formedCleavedCount;
        this.orderChangeCount = orderChangeCount;
        this.stereoChangeCount = stereoChangeCount;
        this.formedCleavedBonds = Collections.unmodifiableList(formedCleavedBonds);
        this.orderChangedBonds = Collections.unmodifiableList(orderChangedBonds);
        this.stereoChangedBonds = Collections.unmodifiableList(stereoChangedBonds);
        this.reactionCentreFingerprint = Collections.unmodifiableList(reactionCentreFingerprint);
        this.algorithmUsed = algorithmUsed;
        this.reactionSignature = buildReactionSignature();
        this.canonicalHash = ReactionCanonicalizer.computeCanonicalHash(
                formedCleavedBonds, orderChangedBonds, stereoChangedBonds, reactionCentreFingerprint);
    }

    /** Original input SMILES */
    public String getInputSmiles() { return inputSmiles; }

    /** Mapped reaction SMILES with atom-atom mapping numbers */
    public String getMappedSmiles() { return mappedSmiles; }

    /** Number of bonds formed or cleaved */
    public int getFormedCleavedCount() { return formedCleavedCount; }

    /** Number of bond order changes */
    public int getOrderChangeCount() { return orderChangeCount; }

    /** Number of stereochemistry changes */
    public int getStereoChangeCount() { return stereoChangeCount; }

    /** Total bond changes (formed/cleaved + order changes) */
    public int getTotalBondChanges() { return formedCleavedCount + orderChangeCount; }

    /** Whether mapping was successful */
    public boolean isMapped() { return mappedSmiles != null; }

    /** Bond formation/cleavage patterns, e.g. ["C-O:1", "O-H:-1"] */
    public List<String> getFormedCleavedBonds() { return formedCleavedBonds; }

    /** Bond order change patterns, e.g. ["C=C:1"] */
    public List<String> getOrderChangedBonds() { return orderChangedBonds; }

    /** Stereo change patterns */
    public List<String> getStereoChangedBonds() { return stereoChangedBonds; }

    /** Reaction centre fingerprint — patterns at the reaction centre */
    public List<String> getReactionCentreFingerprint() { return reactionCentreFingerprint; }

    /** Algorithm that produced this mapping (RINGS, MIN, MAX, MIXTURE) */
    public String getAlgorithm() { return algorithmUsed; }

    /**
     * Canonical, hierarchical reaction signature (R-string).
     * Deterministic, invariant, and searchable. Encodes the complete
     * electron shift pattern as a canonical string.
     *
     * Format: FC[patterns]|OC[patterns]|SC[patterns]|RC[patterns]
     * Where FC=formed/cleaved, OC=order changes, SC=stereo, RC=reaction centre.
     * Patterns are sorted alphabetically within each level.
     *
     * Two reactions with identical signatures have identical bond changes
     * (same R-matrix in the Dugundji-Ugi model / Leber canonicalization).
     *
     * @return canonical reaction signature string, or empty string if unmapped
     */
    public String getReactionSignature() { return reactionSignature; }

    /**
     * Canonical WL graph hash of the ITS (Imaginary Transition State) graph.
     * SHA-256 based, permutation-invariant, deterministic.
     * Two reactions with identical hashes have identical bond change patterns.
     *
     * Use for database indexing, deduplication, and exact-match search.
     *
     * @return 64-character hex hash string
     */
    public String getCanonicalHash() { return canonicalHash; }

    /**
     * Build the canonical reaction signature from sorted fingerprint patterns.
     * Strips weights, sorts alphabetically, joins with semicolons.
     */
    private String buildReactionSignature() {
        if (!isMapped()) return "";
        StringBuilder sb = new StringBuilder();
        sb.append("FC[").append(canonicalPatterns(formedCleavedBonds)).append("]");
        sb.append("|OC[").append(canonicalPatterns(orderChangedBonds)).append("]");
        sb.append("|SC[").append(canonicalPatterns(stereoChangedBonds)).append("]");
        sb.append("|RC[").append(canonicalPatterns(reactionCentreFingerprint)).append("]");
        return sb.toString();
    }

    /**
     * Extract pattern names (strip weights), sort, join with semicolons.
     */
    private static String canonicalPatterns(List<String> features) {
        List<String> patterns = new java.util.ArrayList<>();
        for (String f : features) {
            int colon = f.lastIndexOf(':');
            patterns.add(colon > 0 ? f.substring(0, colon) : f);
        }
        Collections.sort(patterns);
        return String.join(";", patterns);
    }

    /**
     * Compute Tanimoto similarity between this reaction and another
     * based on bond change fingerprints. Returns 0.0 (no overlap) to 1.0 (identical).
     *
     * @param other another ReactionResult to compare against
     * @return Tanimoto similarity coefficient
     */
    public double similarity(ReactionResult other) {
        if (other == null || !this.isMapped() || !other.isMapped()) return 0.0;
        return tanimoto(this.getAllFingerprints(), other.getAllFingerprints());
    }

    /**
     * Get all fingerprint features as a combined set (for similarity).
     * Keeps the full "PATTERN:WEIGHT" strings so that "C-O:2" and "C-O:1"
     * are treated as distinct, giving accurate Tanimoto scores for reactions
     * that differ only in stoichiometry.
     */
    private Set<String> getAllFingerprints() {
        Set<String> all = new HashSet<>();
        all.addAll(formedCleavedBonds);
        all.addAll(orderChangedBonds);
        all.addAll(stereoChangedBonds);
        all.addAll(reactionCentreFingerprint);
        return all;
    }

    /**
     * Tanimoto coefficient: |A ∩ B| / |A ∪ B|
     */
    private static double tanimoto(Set<String> a, Set<String> b) {
        if (a.isEmpty() && b.isEmpty()) return 1.0;
        if (a.isEmpty() || b.isEmpty()) return 0.0;
        Set<String> intersection = new HashSet<>(a);
        intersection.retainAll(b);
        Set<String> union = new HashSet<>(a);
        union.addAll(b);
        return (double) intersection.size() / union.size();
    }

    @Override
    public String toString() {
        return "ReactionResult{" +
                "mapped=" + isMapped() +
                ", algorithm=" + algorithmUsed +
                ", bondChanges=" + getTotalBondChanges() +
                ", formed/cleaved=" + formedCleavedBonds +
                ", orderChanges=" + orderChangedBonds +
                ", stereoChanges=" + stereoChangedBonds +
                ", reactionCentre=" + reactionCentreFingerprint +
                ", signature=" + reactionSignature +
                '}';
    }
}
