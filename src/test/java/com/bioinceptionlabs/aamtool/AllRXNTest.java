/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.aamtool;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import org.openscience.cdk.interfaces.IReaction;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.MappingUtility;

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
                }
            } catch (Exception e) {
                fail++;
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
                }
            } catch (Exception e) {
                fail++;
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
                }
            } catch (Exception e) {
                fail++;
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
                }
            } catch (Exception e) {
                fail++;
            }
        }
        System.out.println("Other: " + success + "/" + files.size()
                + " mapped (" + fail + " failures)");
        assertTrue("Other mapping success should be > 60%",
                success > files.size() * 0.60);
    }

    /**
     * List all .rxn files in a resource directory.
     */
    private List<String> listRXNFiles(String resourceDir) {
        List<String> files = new ArrayList<>();
        try {
            // Read directory listing from classpath
            java.io.File dir = new java.io.File(
                    getClass().getClassLoader().getResource(resourceDir).toURI());
            if (dir.isDirectory()) {
                for (java.io.File f : dir.listFiles()) {
                    if (f.getName().endsWith(".rxn")) {
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
