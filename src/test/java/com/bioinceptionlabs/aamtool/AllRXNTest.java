/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.MappingUtility;
import com.bioinceptionlabs.testgroups.FullRegression;

import static com.bioinceptionlabs.reactionblast.tools.TestUtility.*;
import static org.junit.Assert.assertTrue;

/**
 * Comprehensive test: maps ALL 603 RXN files from test resources.
 * Run with: mvn test -Dtest=AllRXNTest
 *
 * Categories:
 * - KEGG (98 files): metabolic reactions from KEGG database
 * - MACiE (358 files): mechanistic enzyme reactions from MACiE
 * - RHEA (55 files): curated reactions from Rhea
 * - BRENDA (5 files): BRENDA enzyme database reactions
 * - Bug (39 files): previously problematic edge cases
 * - Other (48 files): miscellaneous test reactions
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
@Category(FullRegression.class)
public class AllRXNTest extends MappingUtility {

    /**
     * Test ALL KEGG reactions (98 files). Core metabolic reactions.
     */
    @Test
    public void testAllKEGG() throws Exception {
        List<String> files = listRXNFiles("rxn/kegg/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, KEGG_RXN_DIR);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  KEGG no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  KEGG error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("KEGG: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("KEGG mapping success should be > 80%",
                success > files.size() * 0.80);
    }

    /**
     * Test ALL MACiE reactions (358 files). Mechanistic enzyme reactions.
     * These are complex and some may fail — we test coverage rate.
     */
    @Test
    public void testAllMACiE() throws Exception {
        List<String> files = listRXNFiles("rxn/macie/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, MACIE_RXN);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  MACiE no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  MACiE error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("MACiE: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("MACiE mapping success should be > 70%",
                success > files.size() * 0.70);
    }

    /**
     * Test ALL RHEA reactions (55 files). Curated reference reactions.
     */
    @Test
    public void testAllRHEA() throws Exception {
        List<String> files = listRXNFiles("rxn/rhea/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, RHEA_RXN_DIR);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  RHEA no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  RHEA error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("RHEA: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("RHEA mapping success should be > 80%",
                success > files.size() * 0.80);
    }

    /**
     * Test ALL BRENDA reactions (5 files).
     */
    @Test
    public void testAllBRENDA() throws Exception {
        List<String> files = listRXNFiles("rxn/brenda/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, BRENDA_RXN_DIR);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  BRENDA no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  BRENDA error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("BRENDA: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("BRENDA mapping success should be > 60%",
                success > files.size() * 0.60);
    }

    /**
     * Test bug/edge-case reactions (39 files).
     */
    @Test
    public void testAllBugCases() throws Exception {
        List<String> files = listRXNFiles("rxn/bug/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, BUG_RXN_DIR);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  Bug no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  Bug error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("Bug cases: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        // Bug cases are known hard — lower threshold
        assertTrue("Bug case mapping success should be > 50%",
                success > files.size() * 0.50);
    }

    /**
     * Test other/miscellaneous reactions (48 files).
     */
    @Test
    public void testAllOther() throws Exception {
        List<String> files = listRXNFiles("rxn/other/");
        int success = 0, fail = 0;
        for (String rxnFile : files) {
            String id = rxnFile.replace(".rxn", "");
            try {
                ReactionMechanismTool rmt = testReactions(id, OTHER_RXN);
                if (rmt != null && rmt.getSelectedSolution() != null) {
                    success++;
                } else {
                    fail++;
                    System.out.println("  Other no solution: " + id);
                }
            } catch (Exception e) {
                fail++;
                System.out.println("  Other error: " + id + " - " + e.getMessage());
            }
        }
        System.out.println("Other: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("Other mapping success should be > 60%",
                success > files.size() * 0.60);
    }

    /**
     * Known malformed or unsupported test files:
     * - M0354.ov.rxn: mixed RXN + $RDFILE format (corrupted header)
     * - k.rxn: MDL V3000 format (reader only supports V2000)
     * - 200.rxn: atoms/bonds on single line (missing newlines)
     * - Complex.rxn: 31 reactants + 41 products (not a valid single reaction)
     */
    private static final java.util.Set<String> KNOWN_MALFORMED = new java.util.HashSet<>(
            java.util.Arrays.asList("M0354.ov.rxn", "k.rxn", "200.rxn", "Complex.rxn"));

    /**
     * List all .rxn files in a resource directory, excluding known malformed files.
     */
    private List<String> listRXNFiles(String resourceDir) {
        List<String> files = new ArrayList<>();
        try {
            // Read directory listing from classpath
            java.io.File dir = new java.io.File(
                    getClass().getClassLoader().getResource(resourceDir).toURI());
            if (dir.isDirectory()) {
                for (java.io.File f : dir.listFiles()) {
                    if (f.getName().endsWith(".rxn") && !KNOWN_MALFORMED.contains(f.getName())) {
                        files.add(f.getName());
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("Cannot list " + resourceDir + ": " + e.getMessage());
        }
        java.util.Collections.sort(files);
        return files;
    }
}
