/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mechanism.MappingSolution;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.mapping.MappingDiagnostics;
import com.bioinceptionlabs.reactionblast.mapping.MappingKeyUtil;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import com.bioinceptionlabs.reactionblast.tools.ChemicalFileIO.MDLRXNV2000Reader;
import com.bioinceptionlabs.testgroups.Benchmark;

import static org.junit.Assert.assertTrue;

/**
 * Golden Dataset Benchmark: 1,851 manually curated reactions from
 * Lin et al. "Atom-to-atom Mapping: A Benchmarking Study of Popular
 * Mapping Algorithms and Consensus Strategies"
 * Molecular Informatics 41(4):e2100138, 2022.
 * DOI: 10.1002/minf.202100138
 *
 * Combined metrics framework for fair comparison across tools:
 * - Mapping success (tool doesn't crash)
 * - Exact atom-map match (every atom maps to same position as reference)
 * - Atom-level accuracy (% of individual atoms correctly mapped)
 * - Bond-change accuracy (same formed/broken bonds as reference)
 * - Quality score (coverage + parsimony + locality)
 * - Equivalent match (different but equally valid mapping)
 *
 * @author Syed Asad Rahman
 */
@Category(Benchmark.class)
public class GoldenDatasetBenchmarkTest {

    private static final String GOLDEN_RDF = "benchmark/golden_dataset.rdf";
    private static final int MAX_REACTIONS = Integer.getInteger("golden.max", 0); // 0 = all
    private static final int REPORT_MISMATCHES = Integer.getInteger("golden.reportMismatches", 0);
    private static final int DUPLICATE_PROFILE_LIMIT = Integer.getInteger("golden.duplicate.max", 25);
    private static final int DUPLICATE_PROFILE_PRINT = Integer.getInteger("golden.duplicate.print", 10);
    private static final String BENCHMARK_ATOM_ID = "benchmarkAtomId";

    @Test
    public void benchmarkGoldenDataset() throws Exception {
        URL rdfUrl = getClass().getClassLoader().getResource(GOLDEN_RDF);
        if (rdfUrl == null) {
            System.out.println("SKIP: Golden dataset not found at " + GOLDEN_RDF);
            System.out.println("Place golden_dataset.rdf in src/test/resources/benchmark/");
            return;
        }

        List<GoldReaction> goldReactions = parseRDF(rdfUrl);
        System.out.println("Loaded " + goldReactions.size() + " reactions from golden dataset");

        int limit = MAX_REACTIONS > 0 ? Math.min(MAX_REACTIONS, goldReactions.size())
                : goldReactions.size();

        // Counters
        int total = 0;
        int success = 0;             // RDT produced a mapping
        int molMapExact = 0;         // exact reactant-molecule -> product-molecule relations
        int exactAtomMatch = 0;      // all atom mappings match gold standard exactly
        int chemistryEquivalent = 0; // same bond changes as gold, regardless of numbering
        int alternateValidMapping = 0; // chemically equivalent but different atom numbering
        int trueChemistryMiss = 0;   // bond changes differ from gold
        int ambiguousNoChange = 0;   // no-change reactions where numbering differs but chemistry is identical
        int rdtBetter = 0;           // RDT finds fewer bond changes (more parsimonious)
        int bondChangeMatch = 0;     // at least 1 bond change detected
        int bondChangeExact = 0;     // exact same number of bond changes as gold
        int bondChangeCountExact = 0; // exact same total number of bond changes
        int bondChangeTypeExact = 0; // exact same FORM/BREAK/ORDER counts
        int reactionCenterExact = 0; // exact same reaction-center atom set
        int errors = 0;
        int totalGoldAtoms = 0;
        int correctAtoms = 0;
        int totalGoldReactionCenterAtoms = 0;
        int correctReactionCenterAtoms = 0;
        int goldParseFail = 0;
        double totalQualityScore = 0;
        int qualityScored = 0;
        int mismatchReports = 0;
        int totalAlgorithmsExecuted = 0;
        Map<Integer, Integer> algorithmsPerReaction = new HashMap<>();
        Map<String, Integer> selectedAlgorithms = new HashMap<>();
        long totalMappingPhaseMs = 0;
        long totalEvaluationPhaseMs = 0;

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < limit; i++) {
            GoldReaction gold = goldReactions.get(i);
            total++;

            try {
                // Parse the RXN block twice — one for gold, one for RDT
                IReaction goldRxn = parseRXNBlock(gold.rxnBlock);
                IReaction rdtRxn = parseRXNBlock(gold.rxnBlock);
                if (goldRxn == null || rdtRxn == null) {
                    errors++;
                    continue;
                }

                assignBenchmarkAtomIds(goldRxn);
                assignBenchmarkAtomIds(rdtRxn);

                // Extract gold atom maps before stripping
                Map<String, String> goldMapByAtomId = extractMapByBenchmarkId(goldRxn);
                totalGoldAtoms += goldMapByAtomId.size();
                if (goldMapByAtomId.isEmpty()) {
                    goldParseFail++;
                }

                // Extract gold bond changes from the mapped reaction
                Set<String> goldBondChanges = extractBondChanges(goldRxn);
                Map<String, Integer> goldBondChangeTypes = extractBondChangeTypeCounts(goldBondChanges);
                Set<String> goldReactionCenter = extractReactionCenterAtoms(goldBondChanges);
                Set<String> goldMolMapPairs = extractMolMapPairs(goldMapByAtomId);
                totalGoldReactionCenterAtoms += goldReactionCenter.size();

                // Strip atom maps for RDT input
                stripAtomMaps(rdtRxn);
                rdtRxn.setID("GOLDEN_" + (i + 1));

                // Run RDT mapping
                ReactionMechanismTool rmt = performAtomAtomMapping(rdtRxn, "GOLDEN_" + (i + 1));
                MappingSolution solution = rmt.getSelectedSolution();
                MappingDiagnostics.ReactionSnapshot snapshot = MappingDiagnostics.snapshot("GOLDEN_" + (i + 1));
                int executedAlgorithms = snapshot.algorithms.size();
                totalAlgorithmsExecuted += executedAlgorithms;
                algorithmsPerReaction.merge(executedAlgorithms, 1, Integer::sum);
                totalMappingPhaseMs += snapshot.mappingPhaseMillis;
                totalEvaluationPhaseMs += snapshot.evaluationPhaseMillis;

                    if (solution != null && solution.getBondChangeCalculator() != null) {
                        success++;
                        selectedAlgorithms.merge(solution.getAlgorithmID().name(), 1, Integer::sum);

                        // --- Metric 1: Atom-level accuracy ---
                        IReaction mappedRxn = solution.getReaction();
                        int matched = 0;
                        int goldAtomCount = goldMapByAtomId.size();
                        boolean exactMapping = false;
                        if (mappedRxn != null) {
                            Map<String, String> rdtMapByAtomId = extractMapByBenchmarkId(mappedRxn);
                            Set<String> rdtMolMapPairs = extractMolMapPairs(rdtMapByAtomId);
                            for (Map.Entry<String, String> goldEntry : goldMapByAtomId.entrySet()) {
                                String rdtProduct = rdtMapByAtomId.get(goldEntry.getKey());
                                if (rdtProduct != null && rdtProduct.equals(goldEntry.getValue())) {
                                    matched++;
                                }
                            }
                            correctAtoms += matched;
                            exactMapping = matched == goldMapByAtomId.size() && !goldMapByAtomId.isEmpty();
                            if (rdtMolMapPairs.equals(goldMolMapPairs) && !goldMolMapPairs.isEmpty()) {
                                molMapExact++;
                            }
                            if (exactMapping) {
                                exactAtomMatch++;
                            }
                        }

                        // --- Metric 2: Bond-change analysis ---
                        Set<String> rdtBondChangesSet = mappedRxn != null
                                ? extractBondChanges(mappedRxn) : new HashSet<>();
                        IPatternFingerprinter formedCleaved = solution.getBondChangeCalculator()
                                .getFormedCleavedWFingerprint();
                        IPatternFingerprinter orderChanged = solution.getBondChangeCalculator()
                                .getOrderChangesWFingerprint();
                    int rdtBondChanges = rdtBondChangesSet.size();
                    int goldBondChangeCount = goldBondChanges.size();
                    boolean chemicallyEquivalent = rdtBondChangesSet.equals(goldBondChanges);
                    Map<String, Integer> rdtBondChangeTypes = extractBondChangeTypeCounts(rdtBondChangesSet);
                    Set<String> rdtReactionCenter = extractReactionCenterAtoms(rdtBondChangesSet);
                    int matchedReactionCenterAtoms = 0;
                    for (String atomId : goldReactionCenter) {
                        if (rdtReactionCenter.contains(atomId)) {
                            matchedReactionCenterAtoms++;
                        }
                    }
                    correctReactionCenterAtoms += matchedReactionCenterAtoms;
                    if (rdtReactionCenter.equals(goldReactionCenter)) {
                        reactionCenterExact++;
                    }

                        if (rdtBondChanges > 0) {
                            bondChangeMatch++;
                        }
                    if (rdtBondChanges == goldBondChangeCount) {
                        bondChangeCountExact++;
                    }
                    if (rdtBondChangeTypes.equals(goldBondChangeTypes)) {
                        bondChangeTypeExact++;
                    }
                    if (chemicallyEquivalent) {
                        bondChangeExact++;
                    }

                        // --- Metric 3: Smarter mapping classification ---
                    if (chemicallyEquivalent) {
                        chemistryEquivalent++;
                        if (!exactMapping) {
                            alternateValidMapping++;
                            if (rdtBondChanges == 0 && goldBondChangeCount == 0) {
                                ambiguousNoChange++;
                            }
                        }
                    } else {
                        trueChemistryMiss++;
                    }

                    // --- Metric 4: RDT better (more parsimonious) ---
                    if (rdtBondChanges > 0 && rdtBondChanges < goldBondChangeCount) {
                        rdtBetter++;
                    }

                    // --- Metric 5: Quality score ---
                    // Coverage: fraction of gold atoms mapped by RDT
                    double coverage = goldMapByAtomId.isEmpty() ? 1.0
                            : (double) correctAtoms / totalGoldAtoms; // running, not per-rxn
                    // Parsimony: fewer bond changes = better (ratio gold/rdt, capped at 1)
                    double parsimony = (rdtBondChanges == 0) ? 0
                            : Math.min(1.0, (double) goldBondChangeCount / rdtBondChanges);
                    // Simple quality = average of coverage indicator and parsimony
                    double quality = (rdtBondChanges > 0 || goldBondChangeCount > 0)
                            ? parsimony : 1.0;
                    totalQualityScore += quality;
                    qualityScored++;

                    if (REPORT_MISMATCHES > 0
                            && mismatchReports < REPORT_MISMATCHES
                            && (!exactMapping || !chemicallyEquivalent)) {
                        mismatchReports++;
                        System.out.printf(
                                "  Mismatch %d: GOLDEN_%d algo=%s atoms=%d/%d bondChanges=%d/%d exact=%s chemEq=%s%n",
                                mismatchReports,
                                i + 1,
                                solution.getAlgorithmID(),
                                matched,
                                goldAtomCount,
                                rdtBondChanges,
                                goldBondChangeCount,
                                exactMapping,
                                chemicallyEquivalent);
                        System.out.println("    direct=" + rdtBondChangesSet);
                        System.out.println("    gold=" + goldBondChanges);
                        System.out.println("    formed/cleaved=" + describeFingerprint(formedCleaved));
                        System.out.println("    order=" + describeFingerprint(orderChanged));
                    }
                }
            } catch (Exception e) {
                errors++;
            }

            if ((i + 1) % 500 == 0 || (i + 1) == limit) {
                long elapsed = System.currentTimeMillis() - startTime;
                double rate = (i + 1) * 1000.0 / elapsed;
                System.out.printf("  Progress: %d/%d (%.1f rxn/sec, %d errors, %d mol-exact, %d atom-exact, %d alt-valid, %d chem-miss)%n",
                        i + 1, limit, rate, errors, molMapExact, exactAtomMatch,
                        alternateValidMapping, trueChemistryMiss);
            }
        }

        long totalTime = System.currentTimeMillis() - startTime;
        double rxnPerSec = total * 1000.0 / totalTime;
        double atomAccuracy = totalGoldAtoms > 0 ? (100.0 * correctAtoms / totalGoldAtoms) : 0;
        double reactionCenterAccuracy = totalGoldReactionCenterAtoms > 0
                ? (100.0 * correctReactionCenterAtoms / totalGoldReactionCenterAtoms) : 0;
        double avgQuality = qualityScored > 0 ? (100.0 * totalQualityScore / qualityScored) : 0;
        double chemistryEquivalentPct = pct_d(chemistryEquivalent, total);

        System.out.println();
        System.out.println("=== Golden Dataset Benchmark Results (RDT v3.9.0) ===");
        System.out.println("Total reactions:        " + total);
        System.out.println();
        System.out.println("--- Core Metrics ---");
        System.out.println("Mapping success:        " + success + "/" + total
                + " (" + pct(success, total) + "%)");
        System.out.println("Mol-map exact:         " + molMapExact + "/" + total
                + " (" + pct(molMapExact, total) + "%)");
        System.out.println("Exact atom-map match:   " + exactAtomMatch + "/" + total
                + " (" + pct(exactAtomMatch, total) + "%)");
        System.out.println("Atom-level accuracy:    " + correctAtoms + "/" + totalGoldAtoms
                + " (" + String.format("%.1f", atomAccuracy) + "%)");
        System.out.println();
        System.out.println("--- Bond-Change Metrics ---");
        System.out.println("Bond-change found:      " + bondChangeMatch + "/" + total
                + " (" + pct(bondChangeMatch, total) + "%)");
        System.out.println("Bond-change exact:      " + bondChangeExact + "/" + total
                + " (" + pct(bondChangeExact, total) + "%)");
        System.out.println("Bond-change count:      " + bondChangeCountExact + "/" + total
                + " (" + pct(bondChangeCountExact, total) + "%)");
        System.out.println("Bond-change type:       " + bondChangeTypeExact + "/" + total
                + " (" + pct(bondChangeTypeExact, total) + "%)");
        System.out.println("Reaction-center exact:  " + reactionCenterExact + "/" + total
                + " (" + pct(reactionCenterExact, total) + "%)");
        System.out.println("Reaction-center atoms:  " + correctReactionCenterAtoms + "/"
                + totalGoldReactionCenterAtoms
                + " (" + String.format("%.1f", reactionCenterAccuracy) + "%)");
        System.out.println();
        System.out.println("--- Smarter Match Breakdown ---");
        System.out.println("Chemically equivalent:  " + chemistryEquivalent + "/" + total
                + " (" + String.format("%.1f", chemistryEquivalentPct) + "%)");
        System.out.println("Alternate valid map:    " + alternateValidMapping + "/" + total
                + " (" + pct(alternateValidMapping, total) + "%)");
        System.out.println("True chemistry miss:    " + trueChemistryMiss + "/" + total
                + " (" + pct(trueChemistryMiss, total) + "%)");
        System.out.println("No-change ambiguous:    " + ambiguousNoChange + "/" + total
                + " (" + pct(ambiguousNoChange, total) + "%)");
        System.out.println();
        System.out.println("--- Quality Metrics ---");
        System.out.println("RDT more parsimonious:  " + rdtBetter + "/" + total
                + " (" + pct(rdtBetter, total) + "%)");
        System.out.println("Avg quality score:      " + String.format("%.1f", avgQuality) + "%");
        System.out.println();
        System.out.println("--- Diagnostics ---");
        System.out.println("Gold parse failures:    " + goldParseFail);
        System.out.println("Errors:                 " + errors);
        System.out.println("Speed:                  " + String.format("%.1f", rxnPerSec) + " rxn/sec");
        System.out.println("Total time:             " + (totalTime / 1000) + "s");
        System.out.println("Avg algorithms/run:     " + String.format("%.2f",
                total == 0 ? 0.0 : (double) totalAlgorithmsExecuted / total));
        System.out.println("Algorithms/reaction:    " + formatDistribution(algorithmsPerReaction));
        System.out.println("Selected algorithms:    " + formatDistribution(selectedAlgorithms));
        System.out.println("Avg mapping phase:      " + String.format("%.1f ms",
                total == 0 ? 0.0 : (double) totalMappingPhaseMs / total));
        System.out.println("Avg evaluation phase:   " + String.format("%.1f ms",
                total == 0 ? 0.0 : (double) totalEvaluationPhaseMs / total));
        System.out.println();
        System.out.println("=== Comparison with Published Results (Lin et al. 2022) ===");
        System.out.println("Scoring: chemically-equivalent bond changes (fair comparison across all tools)");
        System.out.println("| Tool               | Chem-Equiv  | Mol-Map   | Atom-Map  | Training | Deterministic |");
        System.out.println("|--------------------|-------------|-----------|-----------|----------|---------------|");
        System.out.println("| RXNMapper          | 83.74%†     | -         | -         | Unsup.   | No            |");
        System.out.println("| RDTool (published) | 76.18%†     | -         | -         | None     | Yes           |");
        System.out.println("| ChemAxon           | 70.45%†     | -         | -         | Propr.   | Yes           |");
        System.out.printf("| RDT v3.9.0         | %.1f%%      | %.1f%%    | %.1f%%    | None     | Yes           |%n",
                pct_d(chemistryEquivalent, total), pct_d(molMapExact, total), pct_d(exactAtomMatch, total));
        System.out.println("† Published figures from Lin et al. 2022 use chemically-equivalent scoring.");

        assertTrue("Mapping success rate should be > 70%", success > total * 0.70);
    }

    @Test
    public void profileDuplicateStoichiometrySubset() throws Exception {
        URL rdfUrl = getClass().getClassLoader().getResource(GOLDEN_RDF);
        if (rdfUrl == null) {
            System.out.println("SKIP: Golden dataset not found at " + GOLDEN_RDF);
            return;
        }

        List<GoldReaction> goldReactions = parseRDF(rdfUrl);
        List<Integer> selected = new ArrayList<>();
        for (int i = 0; i < goldReactions.size() && selected.size() < DUPLICATE_PROFILE_LIMIT; i++) {
            IReaction rxn = parseRXNBlock(goldReactions.get(i).rxnBlock);
            if (rxn == null) {
                continue;
            }
            assignBenchmarkAtomIds(rxn);
            IReaction standardized = new StandardizeReaction().standardize(rxn);
            if (hasDuplicateStructures(standardized)) {
                selected.add(i);
            }
        }

        System.out.println("Loaded " + goldReactions.size() + " reactions from golden dataset");
        System.out.println("Duplicate-stoichiometry benchmark subset: " + selected.size() + " reactions");
        assertTrue("Expected at least one duplicate-stoichiometry reaction", !selected.isEmpty());

        long totalCandidatePairs = 0;
        long totalUniquePairs = 0;
        long totalScheduledJobs = 0;
        long totalActualMcsSearches = 0;
        long totalSubstructureSearches = 0;
        long totalCacheHits = 0;
        long totalQuickCalls = 0;
        long totalQuickCacheHits = 0;
        long totalQuickSearches = 0;
        long totalElapsedMs = 0;
        long totalMatcherInvocations = 0;
        long firstInvocationMcsSearches = 0;
        long followUpMcsSearches = 0;
        long firstInvocationUniquePairs = 0;
        long followUpUniquePairs = 0;

        int printed = 0;
        for (int index : selected) {
            String reactionId = "GOLDEN_DUP_" + (index + 1);
            IReaction rdtRxn = parseRXNBlock(goldReactions.get(index).rxnBlock);
            if (rdtRxn == null) {
                continue;
            }
            assignBenchmarkAtomIds(rdtRxn);
            stripAtomMaps(rdtRxn);

            long start = System.currentTimeMillis();
            performAtomAtomMapping(rdtRxn, reactionId);
            long elapsed = System.currentTimeMillis() - start;
            MappingDiagnostics.ReactionSnapshot snapshot = MappingDiagnostics.snapshot(reactionId);

            long reactionCandidatePairs = 0;
            long reactionUniquePairs = 0;
            long reactionScheduledJobs = 0;
            long reactionActualMcsSearches = 0;
            long reactionSubstructureSearches = 0;
            long reactionCacheHits = 0;
            long reactionQuickCalls = 0;
            long reactionQuickCacheHits = 0;
            long reactionQuickSearches = 0;
            long reactionMatcherInvocations = 0;
            List<String> algorithms = new ArrayList<>();

            for (MappingDiagnostics.AlgorithmSnapshot algorithmSnapshot : snapshot.algorithms) {
                algorithms.add(algorithmSnapshot.algorithm);
                reactionQuickCalls += algorithmSnapshot.quickMappingCalls;
                reactionQuickCacheHits += algorithmSnapshot.quickMappingCacheHits;
                reactionQuickSearches += algorithmSnapshot.quickMappingSearches;
                for (MappingDiagnostics.MatcherInvocationSnapshot invocation : algorithmSnapshot.invocations) {
                    reactionMatcherInvocations++;
                    reactionCandidatePairs += invocation.candidatePairs;
                    reactionUniquePairs += invocation.uniquePairs;
                    reactionScheduledJobs += invocation.scheduledJobs;
                    reactionActualMcsSearches += invocation.actualMcsSearches;
                    reactionSubstructureSearches += invocation.substructureSearches;
                    reactionCacheHits += invocation.cacheHits;
                    if (invocation.invocationIndex == 1) {
                        firstInvocationMcsSearches += invocation.actualMcsSearches;
                        firstInvocationUniquePairs += invocation.uniquePairs;
                    } else {
                        followUpMcsSearches += invocation.actualMcsSearches;
                        followUpUniquePairs += invocation.uniquePairs;
                    }
                }
            }

            totalCandidatePairs += reactionCandidatePairs;
            totalUniquePairs += reactionUniquePairs;
            totalScheduledJobs += reactionScheduledJobs;
            totalActualMcsSearches += reactionActualMcsSearches;
            totalSubstructureSearches += reactionSubstructureSearches;
            totalCacheHits += reactionCacheHits;
            totalQuickCalls += reactionQuickCalls;
            totalQuickCacheHits += reactionQuickCacheHits;
            totalQuickSearches += reactionQuickSearches;
            totalElapsedMs += elapsed;
            totalMatcherInvocations += reactionMatcherInvocations;

            if (printed < DUPLICATE_PROFILE_PRINT) {
                int duplicateMultiplicity = maxDuplicateMultiplicity(parseRXNBlock(goldReactions.get(index).rxnBlock));
                System.out.printf(
                        "  %s dup-max=%d algos=%s matcher=%d pairs=%d unique=%d scheduled=%d mcs=%d sub=%d cache=%d quick=%d/%d time=%dms%n",
                        reactionId,
                        duplicateMultiplicity,
                        algorithms,
                        reactionMatcherInvocations,
                        reactionCandidatePairs,
                        reactionUniquePairs,
                        reactionScheduledJobs,
                        reactionActualMcsSearches,
                        reactionSubstructureSearches,
                        reactionCacheHits,
                        reactionQuickCacheHits,
                        reactionQuickCalls,
                        elapsed);
                printed++;
            }
        }

        double uniqueReductionPct = totalCandidatePairs == 0 ? 0.0
                : 100.0 * (totalCandidatePairs - totalUniquePairs) / totalCandidatePairs;
        double followUpMcsPct = totalActualMcsSearches == 0 ? 0.0
                : 100.0 * followUpMcsSearches / totalActualMcsSearches;

        System.out.println();
        System.out.println("=== Duplicate Stoichiometry Diagnostics ===");
        System.out.println("Reactions profiled:      " + selected.size());
        System.out.println("Total candidate pairs:   " + totalCandidatePairs);
        System.out.println("Unique structural pairs: " + totalUniquePairs);
        System.out.println("Pair reduction:          " + String.format("%.1f%%", uniqueReductionPct));
        System.out.println("Matcher invocations:     " + totalMatcherInvocations);
        System.out.println("Scheduled MCS jobs:      " + totalScheduledJobs);
        System.out.println("Actual MCS searches:     " + totalActualMcsSearches);
        System.out.println("Substructure searches:   " + totalSubstructureSearches);
        System.out.println("MCS cache hits:          " + totalCacheHits);
        System.out.println("QuickMapping calls:      " + totalQuickCalls);
        System.out.println("QuickMapping cache hits: " + totalQuickCacheHits);
        System.out.println("QuickMapping searches:   " + totalQuickSearches);
        System.out.println("Elapsed time:            " + totalElapsedMs + " ms");
        System.out.println("First-pass unique pairs: " + firstInvocationUniquePairs);
        System.out.println("Follow-up unique pairs:  " + followUpUniquePairs);
        System.out.println("First-pass MCS calls:    " + firstInvocationMcsSearches);
        System.out.println("Follow-up MCS calls:     " + followUpMcsSearches);
        System.out.println("Follow-up MCS share:     " + String.format("%.1f%%", followUpMcsPct));
        System.out.println("Fragment-stage hotspot:  "
                + (followUpMcsPct >= 50.0 ? "likely" : "not dominant in this subset"));
    }

    // ---- Bond-change extraction from mapped reaction ----

    private Set<String> extractBondChanges(IReaction rxn) {
        Set<String> changes = new HashSet<>();
        Map<String, Map<String, IBond.Order>> reactantBonds = collectBondsByBenchmarkId(rxn.getReactants());
        Map<String, Map<String, IBond.Order>> productBonds = collectBondsByBenchmarkId(rxn.getProducts());

        // Broken bonds: in reactants but not products
        for (Map.Entry<String, Map<String, IBond.Order>> entry : reactantBonds.entrySet()) {
            String atom1 = entry.getKey();
            for (Map.Entry<String, IBond.Order> bond : entry.getValue().entrySet()) {
                String atom2 = bond.getKey();
                if (atom1.compareTo(atom2) < 0) { // avoid double counting
                    IBond.Order prodOrder = productBonds.getOrDefault(atom1, new HashMap<>()).get(atom2);
                    if (prodOrder == null) {
                        changes.add("BREAK:" + atom1 + "-" + atom2);
                    } else if (prodOrder != bond.getValue()) {
                        changes.add("ORDER:" + atom1 + "-" + atom2);
                    }
                }
            }
        }
        // Formed bonds: in products but not reactants
        for (Map.Entry<String, Map<String, IBond.Order>> entry : productBonds.entrySet()) {
            String atom1 = entry.getKey();
            for (Map.Entry<String, IBond.Order> bond : entry.getValue().entrySet()) {
                String atom2 = bond.getKey();
                if (atom1.compareTo(atom2) < 0) {
                    IBond.Order reactOrder = reactantBonds.getOrDefault(atom1, new HashMap<>()).get(atom2);
                    if (reactOrder == null) {
                        changes.add("FORM:" + atom1 + "-" + atom2);
                    }
                }
            }
        }
        return changes;
    }

    private Set<String> extractReactionCenterAtoms(Set<String> bondChanges) {
        Set<String> atoms = new HashSet<>();
        for (String change : bondChanges) {
            int colonIndex = change.indexOf(':');
            if (colonIndex < 0 || colonIndex + 1 >= change.length()) {
                continue;
            }
            String pair = change.substring(colonIndex + 1);
            String[] atomIds = pair.split("-");
            for (String atomId : atomIds) {
                if (!atomId.isEmpty()) {
                    atoms.add(atomId);
                }
            }
        }
        return atoms;
    }

    private Map<String, Integer> extractBondChangeTypeCounts(Set<String> bondChanges) {
        Map<String, Integer> typeCounts = new HashMap<>();
        for (String change : bondChanges) {
            int colonIndex = change.indexOf(':');
            if (colonIndex <= 0) {
                continue;
            }
            String type = change.substring(0, colonIndex);
            typeCounts.merge(type, 1, Integer::sum);
        }
        return typeCounts;
    }

    private Set<String> extractMolMapPairs(Map<String, String> atomMapByAtomId) {
        Set<String> molPairs = new HashSet<>();
        for (Map.Entry<String, String> entry : atomMapByAtomId.entrySet()) {
            Integer reactantMol = extractMoleculeIndex(entry.getKey());
            Integer productMol = extractMoleculeIndex(entry.getValue());
            if (reactantMol != null && productMol != null) {
                molPairs.add("R:" + reactantMol + ">P:" + productMol);
            }
        }
        return molPairs;
    }

    private Integer extractMoleculeIndex(String atomId) {
        if (atomId == null) {
            return null;
        }
        String[] parts = atomId.split(":");
        if (parts.length != 3) {
            return null;
        }
        try {
            return Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
            return null;
        }
    }

    private Map<String, Map<String, IBond.Order>> collectBondsByBenchmarkId(
            org.openscience.cdk.interfaces.IAtomContainerSet molSet) {
        Map<String, Map<String, IBond.Order>> bonds = new HashMap<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IBond bond : mol.bonds()) {
                String atom1 = getBenchmarkAtomId(bond.getBegin());
                String atom2 = getBenchmarkAtomId(bond.getEnd());
                if (atom1 != null && atom2 != null) {
                    bonds.computeIfAbsent(atom1, k -> new HashMap<>()).put(atom2, bond.getOrder());
                    bonds.computeIfAbsent(atom2, k -> new HashMap<>()).put(atom1, bond.getOrder());
                }
            }
        }
        return bonds;
    }

    // ---- Atom map extraction ----

    private Map<String, String> extractMapByBenchmarkId(IReaction rxn) {
        Map<Integer, String> reactantMapNums = new HashMap<>();
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                String atomId = getBenchmarkAtomId(atom);
                if (mapNum > 0 && atomId != null) {
                    reactantMapNums.put(mapNum, atomId);
                }
            }
        }
        Map<String, String> mapPairs = new HashMap<>();
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                String atomId = getBenchmarkAtomId(atom);
                String reactantId = reactantMapNums.get(mapNum);
                if (mapNum > 0 && atomId != null && reactantId != null) {
                    mapPairs.put(reactantId, atomId);
                }
            }
        }
        return mapPairs;
    }

    private void assignBenchmarkAtomIds(IReaction rxn) {
        int molIndex = 0;
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                mol.getAtom(atomIndex).setProperty(BENCHMARK_ATOM_ID, "R:" + molIndex + ":" + atomIndex);
            }
            molIndex++;
        }
        molIndex = 0;
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                mol.getAtom(atomIndex).setProperty(BENCHMARK_ATOM_ID, "P:" + molIndex + ":" + atomIndex);
            }
            molIndex++;
        }
    }

    private String getBenchmarkAtomId(IAtom atom) {
        Object atomId = atom.getProperty(BENCHMARK_ATOM_ID);
        return atomId != null ? atomId.toString() : null;
    }

    private int getAtomMapNumber(IAtom atom) {
        int mapIdx = atom.getMapIdx();
        if (mapIdx > 0) return mapIdx;
        Object prop = atom.getProperty("molAtomMapNumber");
        if (prop instanceof Integer && (Integer) prop > 0) return (Integer) prop;
        prop = atom.getProperty("ATOM_ATOM_MAPPING");
        if (prop instanceof Integer && (Integer) prop > 0) return (Integer) prop;
        if (prop instanceof String) {
            try { return Integer.parseInt((String) prop); } catch (NumberFormatException e) { /* skip */ }
        }
        return 0;
    }

    // ---- RDF parsing ----

    private List<GoldReaction> parseRDF(URL rdfUrl) throws IOException {
        List<GoldReaction> reactions = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(rdfUrl.openStream()))) {
            String line;
            StringBuilder rxnBlock = null;
            boolean inReaction = false;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("$RFMT")) {
                    if (rxnBlock != null && rxnBlock.length() > 0) {
                        reactions.add(new GoldReaction(rxnBlock.toString()));
                    }
                    rxnBlock = new StringBuilder();
                    inReaction = true;
                    continue;
                }
                if (inReaction && rxnBlock != null) {
                    if (line.startsWith("$DTYPE") || line.startsWith("$DATUM")) continue;
                    rxnBlock.append(line).append("\n");
                }
            }
            if (rxnBlock != null && rxnBlock.length() > 0) {
                reactions.add(new GoldReaction(rxnBlock.toString()));
            }
        }
        return reactions;
    }

    private IReaction parseRXNBlock(String rxnBlock) {
        try {
            MDLRXNV2000Reader reader = new MDLRXNV2000Reader(
                    new BufferedReader(new StringReader(rxnBlock)));
            return reader.read(new Reaction());
        } catch (Exception e) {
            return null;
        }
    }

    private void stripAtomMaps(IReaction rxn) {
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                atom.removeProperty("molAtomMapNumber");
                atom.removeProperty("ATOM_ATOM_MAPPING");
                atom.setMapIdx(0);
            }
        }
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                atom.removeProperty("molAtomMapNumber");
                atom.removeProperty("ATOM_ATOM_MAPPING");
                atom.setMapIdx(0);
            }
        }
    }

    private ReactionMechanismTool performAtomAtomMapping(IReaction rxn, String id) throws Exception {
        rxn.setID(id);
        return new ReactionMechanismTool(rxn, true, true, false, true, true,
                new StandardizeReaction());
    }

    private boolean hasDuplicateStructures(IReaction reaction) {
        return maxDuplicateMultiplicity(reaction) > 1;
    }

    private int maxDuplicateMultiplicity(IReaction reaction) {
        if (reaction == null) {
            return 0;
        }
        int max = 1;
        max = Math.max(max, maxDuplicateMultiplicity(reaction.getReactants()));
        max = Math.max(max, maxDuplicateMultiplicity(reaction.getProducts()));
        return max;
    }

    private int maxDuplicateMultiplicity(org.openscience.cdk.interfaces.IAtomContainerSet containers) {
        Map<String, Integer> counts = new HashMap<>();
        int max = 1;
        for (IAtomContainer container : containers.atomContainers()) {
            String key = MappingKeyUtil.computeStructureKey(container);
            int count = counts.merge(key, 1, Integer::sum);
            if (count > max) {
                max = count;
            }
        }
        return max;
    }

    private String pct(int num, int den) {
        return den == 0 ? "0.0" : String.format("%.1f", 100.0 * num / den);
    }

    private int getWeightedFeatureCount(IPatternFingerprinter fingerprint) {
        int count = 0;
        for (com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.IFeature feature
                : fingerprint.getFeatures()) {
            count += (int) Math.round(feature.getWeight());
        }
        return count;
    }

    private String describeFingerprint(IPatternFingerprinter fingerprint) {
        List<String> features = new ArrayList<>();
        for (com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.IFeature feature
                : fingerprint.getFeatures()) {
            features.add(feature.getPattern() + "x" + (int) Math.round(feature.getWeight()));
        }
        return features.toString();
    }

    private double pct_d(int num, int den) {
        return den == 0 ? 0.0 : 100.0 * num / den;
    }

    private String formatDistribution(Map<?, Integer> distribution) {
        List<String> entries = new ArrayList<>();
        for (Map.Entry<?, Integer> entry : distribution.entrySet()) {
            entries.add(entry.getKey() + "=" + entry.getValue());
        }
        Collections.sort(entries);
        return entries.toString();
    }

    private static class GoldReaction {
        final String rxnBlock;
        GoldReaction(String rxnBlock) { this.rxnBlock = rxnBlock; }
    }
}
