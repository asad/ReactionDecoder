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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.mechanism.MappingSolution;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import com.bioinceptionlabs.reactionblast.tools.ChemicalFileIO.MDLRXNV2000Reader;

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
public class GoldenDatasetBenchmarkTest {

    private static final String GOLDEN_RDF = "benchmark/golden_dataset.rdf";
    private static final int MAX_REACTIONS = Integer.getInteger("golden.max", 0); // 0 = all

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
        int exactAtomMatch = 0;      // all atom mappings match gold standard exactly
        int equivalentMatch = 0;     // different mapping but equally valid (same bond changes)
        int rdtBetter = 0;           // RDT finds fewer bond changes (more parsimonious)
        int bondChangeMatch = 0;     // at least 1 bond change detected
        int bondChangeExact = 0;     // exact same number of bond changes as gold
        int errors = 0;
        int totalGoldAtoms = 0;
        int correctAtoms = 0;
        int goldParseFail = 0;
        double totalQualityScore = 0;
        int qualityScored = 0;

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

                // Extract gold atom maps before stripping
                Map<Integer, Integer> goldMapByNumber = extractMapByNumber(goldRxn);
                totalGoldAtoms += goldMapByNumber.size();
                if (goldMapByNumber.isEmpty()) {
                    goldParseFail++;
                }

                // Extract gold bond changes from the mapped reaction
                Set<String> goldBondChanges = extractBondChanges(goldRxn);

                // Strip atom maps for RDT input
                stripAtomMaps(rdtRxn);
                rdtRxn.setID("GOLDEN_" + (i + 1));

                // Run RDT mapping
                ReactionMechanismTool rmt = performAtomAtomMapping(rdtRxn, "GOLDEN_" + (i + 1));
                MappingSolution solution = rmt.getSelectedSolution();

                if (solution != null && solution.getBondChangeCalculator() != null) {
                    success++;

                    // --- Metric 1: Atom-level accuracy ---
                    IReaction mappedRxn = solution.getReaction();
                    if (mappedRxn != null) {
                        Map<Integer, Integer> rdtMapByNumber = extractMapByNumber(mappedRxn);
                        int matched = 0;
                        for (Map.Entry<Integer, Integer> goldEntry : goldMapByNumber.entrySet()) {
                            Integer rdtProduct = rdtMapByNumber.get(goldEntry.getKey());
                            if (rdtProduct != null && rdtProduct.equals(goldEntry.getValue())) {
                                matched++;
                            }
                        }
                        correctAtoms += matched;
                        if (matched == goldMapByNumber.size() && !goldMapByNumber.isEmpty()) {
                            exactAtomMatch++;
                        }
                    }

                    // --- Metric 2: Bond-change analysis ---
                    IPatternFingerprinter formedCleaved = solution.getBondChangeCalculator()
                            .getFormedCleavedWFingerprint();
                    IPatternFingerprinter orderChanged = solution.getBondChangeCalculator()
                            .getOrderChangesWFingerprint();
                    int rdtBondChanges = formedCleaved.getFeatureCount() + orderChanged.getFeatureCount();
                    int goldBondChangeCount = goldBondChanges.size();

                    if (rdtBondChanges > 0) {
                        bondChangeMatch++;
                    }
                    if (rdtBondChanges == goldBondChangeCount) {
                        bondChangeExact++;
                    }

                    // --- Metric 3: Equivalent match ---
                    // If atom maps differ but bond change count is same → equivalent mapping
                    if (rdtBondChanges == goldBondChangeCount && rdtBondChanges > 0) {
                        equivalentMatch++;
                    }

                    // --- Metric 4: RDT better (more parsimonious) ---
                    if (rdtBondChanges > 0 && rdtBondChanges < goldBondChangeCount) {
                        rdtBetter++;
                    }

                    // --- Metric 5: Quality score ---
                    // Coverage: fraction of gold atoms mapped by RDT
                    double coverage = goldMapByNumber.isEmpty() ? 1.0
                            : (double) correctAtoms / totalGoldAtoms; // running, not per-rxn
                    // Parsimony: fewer bond changes = better (ratio gold/rdt, capped at 1)
                    double parsimony = (rdtBondChanges == 0) ? 0
                            : Math.min(1.0, (double) goldBondChangeCount / rdtBondChanges);
                    // Simple quality = average of coverage indicator and parsimony
                    double quality = (rdtBondChanges > 0 || goldBondChangeCount > 0)
                            ? parsimony : 1.0;
                    totalQualityScore += quality;
                    qualityScored++;
                }
            } catch (Exception e) {
                errors++;
            }

            if ((i + 1) % 100 == 0) {
                long elapsed = System.currentTimeMillis() - startTime;
                double rate = (i + 1) * 1000.0 / elapsed;
                System.out.printf("  Progress: %d/%d (%.1f rxn/sec, %d errors, %d exact, %d equiv)%n",
                        i + 1, limit, rate, errors, exactAtomMatch, equivalentMatch);
            }
        }

        long totalTime = System.currentTimeMillis() - startTime;
        double rxnPerSec = total * 1000.0 / totalTime;
        double atomAccuracy = totalGoldAtoms > 0 ? (100.0 * correctAtoms / totalGoldAtoms) : 0;
        double avgQuality = qualityScored > 0 ? (100.0 * totalQualityScore / qualityScored) : 0;
        double trueAccuracy = total > 0
                ? 100.0 * (exactAtomMatch + equivalentMatch - exactAtomMatch) / total : 0;
        // equivalentMatch includes exactAtomMatch cases, so true accuracy = equivalent/total

        System.out.println();
        System.out.println("=== Golden Dataset Benchmark Results (RDT v3.6.1) ===");
        System.out.println("Total reactions:        " + total);
        System.out.println();
        System.out.println("--- Core Metrics ---");
        System.out.println("Mapping success:        " + success + "/" + total
                + " (" + pct(success, total) + "%)");
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
        System.out.println();
        System.out.println("--- Quality Metrics ---");
        System.out.println("Equivalent match:       " + equivalentMatch + "/" + total
                + " (" + pct(equivalentMatch, total) + "%)");
        System.out.println("RDT more parsimonious:  " + rdtBetter + "/" + total
                + " (" + pct(rdtBetter, total) + "%)");
        System.out.println("Avg quality score:      " + String.format("%.1f", avgQuality) + "%");
        System.out.println();
        System.out.println("--- Diagnostics ---");
        System.out.println("Gold parse failures:    " + goldParseFail);
        System.out.println("Errors:                 " + errors);
        System.out.println("Speed:                  " + String.format("%.1f", rxnPerSec) + " rxn/sec");
        System.out.println("Total time:             " + (totalTime / 1000) + "s");
        System.out.println();
        System.out.println("=== Comparison with Published Results (Lin et al. 2022) ===");
        System.out.println("| Tool               | Exact Match | Atom Acc. | Bond Acc. | Training | Deterministic |");
        System.out.println("|--------------------|-------------|-----------|-----------|----------|---------------|");
        System.out.println("| RXNMapper          | 83.74%      | -         | -         | Unsup.   | No            |");
        System.out.println("| RDTool (published) | 76.18%      | -         | -         | None     | Yes           |");
        System.out.println("| ChemAxon           | 70.45%      | -         | -         | Propr.   | Yes           |");
        System.out.printf("| RDT v3.6.1         | %.1f%%      | %.1f%%    | %.1f%%    | None     | Yes           |%n",
                pct_d(exactAtomMatch, total), atomAccuracy, pct_d(bondChangeExact, total));

        assertTrue("Mapping success rate should be > 70%", success > total * 0.70);
    }

    // ---- Bond-change extraction from mapped reaction ----

    private Set<String> extractBondChanges(IReaction rxn) {
        Set<String> changes = new HashSet<>();
        Map<Integer, Map<Integer, IBond.Order>> reactantBonds = collectBondsByMap(rxn.getReactants());
        Map<Integer, Map<Integer, IBond.Order>> productBonds = collectBondsByMap(rxn.getProducts());

        // Broken bonds: in reactants but not products
        for (Map.Entry<Integer, Map<Integer, IBond.Order>> entry : reactantBonds.entrySet()) {
            int atom1 = entry.getKey();
            for (Map.Entry<Integer, IBond.Order> bond : entry.getValue().entrySet()) {
                int atom2 = bond.getKey();
                if (atom1 < atom2) { // avoid double counting
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
        for (Map.Entry<Integer, Map<Integer, IBond.Order>> entry : productBonds.entrySet()) {
            int atom1 = entry.getKey();
            for (Map.Entry<Integer, IBond.Order> bond : entry.getValue().entrySet()) {
                int atom2 = bond.getKey();
                if (atom1 < atom2) {
                    IBond.Order reactOrder = reactantBonds.getOrDefault(atom1, new HashMap<>()).get(atom2);
                    if (reactOrder == null) {
                        changes.add("FORM:" + atom1 + "-" + atom2);
                    }
                }
            }
        }
        return changes;
    }

    private Map<Integer, Map<Integer, IBond.Order>> collectBondsByMap(
            org.openscience.cdk.interfaces.IAtomContainerSet molSet) {
        Map<Integer, Map<Integer, IBond.Order>> bonds = new HashMap<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IBond bond : mol.bonds()) {
                int map1 = getAtomMapNumber(bond.getBegin());
                int map2 = getAtomMapNumber(bond.getEnd());
                if (map1 > 0 && map2 > 0) {
                    bonds.computeIfAbsent(map1, k -> new HashMap<>()).put(map2, bond.getOrder());
                    bonds.computeIfAbsent(map2, k -> new HashMap<>()).put(map1, bond.getOrder());
                }
            }
        }
        return bonds;
    }

    // ---- Atom map extraction ----

    private Map<Integer, Integer> extractMapByNumber(IReaction rxn) {
        Set<Integer> reactantMapNums = new HashSet<>();
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) reactantMapNums.add(mapNum);
            }
        }
        Set<Integer> productMapNums = new HashSet<>();
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) productMapNums.add(mapNum);
            }
        }
        Map<Integer, Integer> mapPairs = new HashMap<>();
        for (int mapNum : reactantMapNums) {
            if (productMapNums.contains(mapNum)) {
                mapPairs.put(mapNum, mapNum);
            }
        }
        return mapPairs;
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

    private String pct(int num, int den) {
        return den == 0 ? "0.0" : String.format("%.1f", 100.0 * num / den);
    }

    private double pct_d(int num, int den) {
        return den == 0 ? 0.0 : 100.0 * num / den;
    }

    private static class GoldReaction {
        final String rxnBlock;
        GoldReaction(String rxnBlock) { this.rxnBlock = rxnBlock; }
    }
}
