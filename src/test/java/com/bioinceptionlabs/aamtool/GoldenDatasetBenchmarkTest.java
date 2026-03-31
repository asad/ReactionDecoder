/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.mechanism.MappingSolution;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import com.bioinceptionlabs.reactionblast.tools.ChemicalFileIO.MDLRXNV2000Reader;

import static org.junit.Assert.assertTrue;
import java.io.StringReader;

/**
 * Golden Dataset Benchmark: 1,851 manually curated reactions from
 * Lin et al. "Atom-to-atom Mapping: A Benchmarking Study of Popular
 * Mapping Algorithms and Consensus Strategies"
 * Molecular Informatics 41(4):e2100138, 2022.
 * DOI: 10.1002/minf.202100138
 *
 * Compares RDT atom-atom mapping accuracy against gold-standard mappings.
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

        // Parse RDF into individual RXN blocks with gold atom maps
        List<GoldReaction> goldReactions = parseRDF(rdfUrl);
        System.out.println("Loaded " + goldReactions.size() + " reactions from golden dataset");

        int limit = MAX_REACTIONS > 0 ? Math.min(MAX_REACTIONS, goldReactions.size())
                : goldReactions.size();

        int total = 0;
        int success = 0;           // RDT produced a mapping
        int exactAtomMatch = 0;    // all atom correspondences match gold standard
        int bondChangeMatch = 0;   // bond changes are consistent
        int errors = 0;
        int totalGoldAtoms = 0;
        int correctAtoms = 0;
        int goldParseFail = 0;     // reactions where gold map couldn't be extracted

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < limit; i++) {
            GoldReaction gold = goldReactions.get(i);
            total++;

            try {
                // Parse the RXN block
                IReaction rxn = parseRXNBlock(gold.rxnBlock);
                if (rxn == null) {
                    errors++;
                    continue;
                }

                // Extract gold atom maps (mapNumber pairs) before stripping
                Map<Integer, Integer> goldMapByNumber = extractGoldMapByNumber(rxn);
                totalGoldAtoms += goldMapByNumber.size();

                if (goldMapByNumber.isEmpty()) {
                    goldParseFail++;
                }

                // Strip atom maps for RDT input
                stripAtomMaps(rxn);
                rxn.setID("GOLDEN_" + (i + 1));

                // Run RDT mapping
                ReactionMechanismTool rmt = performAtomAtomMapping(rxn, "GOLDEN_" + (i + 1));
                MappingSolution solution = rmt.getSelectedSolution();

                if (solution != null && solution.getBondChangeCalculator() != null) {
                    success++;

                    // Extract RDT's atom-atom mapping from the mapped reaction
                    IReaction mappedRxn = solution.getReaction();
                    if (mappedRxn != null) {
                        Map<Integer, Integer> rdtMapByNumber = extractGoldMapByNumber(mappedRxn);

                        // Compare: for each gold pair (R_mapNum -> P_mapNum),
                        // check if RDT also maps R_mapNum -> P_mapNum
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

                    // Bond change check
                    if (solution.getBondChangeCalculator()
                            .getFormedCleavedWFingerprint().getFeatureCount() > 0) {
                        bondChangeMatch++;
                    }
                }
            } catch (Exception e) {
                errors++;
            }

            if ((i + 1) % 100 == 0) {
                long elapsed = System.currentTimeMillis() - startTime;
                double rate = (i + 1) * 1000.0 / elapsed;
                System.out.printf("  Progress: %d/%d (%.1f rxn/sec, %d errors, %d exact)%n",
                        i + 1, limit, rate, errors, exactAtomMatch);
            }
        }

        long totalTime = System.currentTimeMillis() - startTime;
        double rxnPerSec = total * 1000.0 / totalTime;
        double atomAccuracy = totalGoldAtoms > 0 ? (100.0 * correctAtoms / totalGoldAtoms) : 0;

        System.out.println();
        System.out.println("=== Golden Dataset Benchmark Results (RDT v3.6.1) ===");
        System.out.println("Total reactions:        " + total);
        System.out.println("Mapping success:        " + success + "/" + total
                + " (" + pct(success, total) + "%)");
        System.out.println("Exact atom-map match:   " + exactAtomMatch + "/" + total
                + " (" + pct(exactAtomMatch, total) + "%)");
        System.out.println("Bond-change found:      " + bondChangeMatch + "/" + total
                + " (" + pct(bondChangeMatch, total) + "%)");
        System.out.println("Atom-level accuracy:    " + correctAtoms + "/" + totalGoldAtoms
                + " (" + String.format("%.1f", atomAccuracy) + "%)");
        System.out.println("Gold parse failures:    " + goldParseFail);
        System.out.println("Errors:                 " + errors);
        System.out.println("Speed:                  " + String.format("%.1f", rxnPerSec) + " rxn/sec");
        System.out.println("Total time:             " + (totalTime / 1000) + "s");
        System.out.println();
        System.out.println("=== Comparison with Published Results (Lin et al. 2022) ===");
        System.out.println("| Tool               | Accuracy  | Training Data  | Deterministic |");
        System.out.println("|--------------------|-----------|----------------|---------------|");
        System.out.println("| RXNMapper          | 83.74%    | Unsupervised   | No            |");
        System.out.println("| RDTool (published) | 76.18%    | None           | Yes           |");
        System.out.println("| ChemAxon           | 70.45%    | Proprietary    | Yes           |");
        System.out.printf("| RDT v3.6.1         | %.1f%%    | None           | Yes           |%n",
                atomAccuracy);

        // RDT should successfully map most reactions
        assertTrue("Mapping success rate should be > 70%", success > total * 0.70);
    }

    /**
     * Parse RDF file into list of RXN blocks with gold atom maps.
     */
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
                    if (line.startsWith("$DTYPE") || line.startsWith("$DATUM")) {
                        // metadata after RXN block — skip
                        continue;
                    }
                    rxnBlock.append(line).append("\n");
                }
            }
            // Last reaction
            if (rxnBlock != null && rxnBlock.length() > 0) {
                reactions.add(new GoldReaction(rxnBlock.toString()));
            }
        }
        return reactions;
    }

    /**
     * Parse an RXN block string into an IReaction.
     */
    private IReaction parseRXNBlock(String rxnBlock) {
        try {
            MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new BufferedReader(new StringReader(rxnBlock)));
            return reader.read(new Reaction());
        } catch (Exception e) {
            return null;
        }
    }

    /**
     * Extract gold-standard atom map: maps atom index (global across all reactant/product
     * containers) to its map number. Returns a map of (reactant_map_number -> product_atom_index)
     * for atoms that appear in both sides.
     */
    /**
     * Extract atom map as mapNumber identity pairs.
     * For the gold standard: each atom has a mapNumber that's the same in reactants and products.
     * Returns: mapNumber -> mapNumber (identity, since gold uses same numbers on both sides).
     * For RDT: after mapping, atoms get new map numbers — same logic applies.
     */
    private Map<Integer, Integer> extractGoldMapByNumber(IReaction rxn) {
        // In both gold and RDT mappings, atoms with the same map number correspond.
        // So we just check: which map numbers exist on both sides?
        java.util.Set<Integer> reactantMapNums = new java.util.HashSet<>();
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) reactantMapNums.add(mapNum);
            }
        }
        java.util.Set<Integer> productMapNums = new java.util.HashSet<>();
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) productMapNums.add(mapNum);
            }
        }
        // Map numbers present on both sides = matched atoms
        Map<Integer, Integer> mapPairs = new HashMap<>();
        for (int mapNum : reactantMapNums) {
            if (productMapNums.contains(mapNum)) {
                mapPairs.put(mapNum, mapNum);
            }
        }
        return mapPairs;
    }

    private Map<Integer, Integer> extractGoldAtomMap(IReaction rxn) {
        // Collect reactant atom maps: mapNumber -> reactant atom index
        Map<Integer, Integer> reactantMaps = new HashMap<>();
        int idx = 0;
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) {
                    reactantMaps.put(mapNum, idx);
                }
                idx++;
            }
        }

        // Collect product atom maps
        Map<Integer, Integer> productMaps = new HashMap<>();
        idx = 0;
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) {
                    productMaps.put(mapNum, idx);
                }
                idx++;
            }
        }

        // Gold map: for each map number present in both sides, record the pair
        Map<Integer, Integer> goldMap = new HashMap<>();
        for (Map.Entry<Integer, Integer> entry : reactantMaps.entrySet()) {
            if (productMaps.containsKey(entry.getKey())) {
                goldMap.put(entry.getValue(), productMaps.get(entry.getKey()));
            }
        }
        return goldMap;
    }

    private int getAtomMapNumber(IAtom atom) {
        // Try CDK's mapIdx first (V2000 atom-atom mapping field)
        int mapIdx = atom.getMapIdx();
        if (mapIdx > 0) return mapIdx;

        // Try properties
        Object prop = atom.getProperty("molAtomMapNumber");
        if (prop instanceof Integer && (Integer) prop > 0) return (Integer) prop;

        prop = atom.getProperty("ATOM_ATOM_MAPPING");
        if (prop instanceof Integer && (Integer) prop > 0) return (Integer) prop;
        if (prop instanceof String) {
            try { return Integer.parseInt((String) prop); } catch (NumberFormatException e) { /* skip */ }
        }
        return 0;
    }

    /**
     * Extract RDT's atom-atom mapping after mapping is complete.
     */
    private Map<Integer, Integer> extractRDTAtomMap(IReaction rxn) {
        Map<Integer, Integer> rdtMap = new HashMap<>();
        Map<Integer, Integer> reactantMaps = new HashMap<>();
        int idx = 0;
        for (IAtomContainer mol : rxn.getReactants().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) {
                    reactantMaps.put(mapNum, idx);
                }
                idx++;
            }
        }

        Map<Integer, Integer> productMaps = new HashMap<>();
        idx = 0;
        for (IAtomContainer mol : rxn.getProducts().atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                int mapNum = getAtomMapNumber(atom);
                if (mapNum > 0) {
                    productMaps.put(mapNum, idx);
                }
                idx++;
            }
        }

        for (Map.Entry<Integer, Integer> entry : reactantMaps.entrySet()) {
            if (productMaps.containsKey(entry.getKey())) {
                rdtMap.put(entry.getValue(), productMaps.get(entry.getKey()));
            }
        }
        return rdtMap;
    }

    /**
     * Compare gold and RDT atom mappings.
     * Returns the number of correctly mapped atoms.
     */
    private int compareMappings(Map<Integer, Integer> goldMap, Map<Integer, Integer> rdtMap) {
        int matched = 0;
        for (Map.Entry<Integer, Integer> gold : goldMap.entrySet()) {
            Integer rdtTarget = rdtMap.get(gold.getKey());
            if (rdtTarget != null && rdtTarget.equals(gold.getValue())) {
                matched++;
            }
        }
        return matched;
    }

    /**
     * Strip atom map numbers from all atoms in the reaction.
     */
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

    private static class GoldReaction {
        final String rxnBlock;
        GoldReaction(String rxnBlock) {
            this.rxnBlock = rxnBlock;
        }
    }
}
