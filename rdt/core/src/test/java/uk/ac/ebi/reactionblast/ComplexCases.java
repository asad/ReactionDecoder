/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad@ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast;

import static org.junit.Assert.assertEquals;
import org.junit.Test;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import static uk.ac.ebi.reactionblast.TestUtility.KEGG_RXN_DIR;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ComplexCases extends MappingUtility {

    /*
     ************************
     * COMPLEX CASES
     ************************
     */
    /**
     * Complex case, linear to Ring system; Takes longer to Run
     *
     * MIN, fp ID=R03200:Bond Cleaved and Formed (4) [C%C:6.0, C%O:1.0,
     * C-C:3.0,C-H:3.0]
     *
     * MIN, fp ID=R03200:Bond Order Change (2) [C%C*C=C:3.0, C-C*C=C:1.0]
     *
     * BE 3437.0, Fragment 0
     *
     * MIXTURE, fp ID=R03200:Bond Cleaved and Formed (5) [C%C:5.0, C%O:1.0,
     * C-C:2.0, C-H:5.0, H-O:1.0]
     *
     * BE 2750.0, Fragment 0
     *
     * @throws java.lang.Exception
     */
    @Test
    public void R03200() throws Exception {

        String reactionID = "R03200";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR, false, false);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Complex case,
     * Takes longer to Run
     * @BUG
     * TO DO 
     * Functional group regarrangment
     *
     *
     * @throws Exception
     */
//    @Test
//    public void R03775() throws Exception {
//        String reactionID = "R03775";
//        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//
//        /*
//         * Expected Solution
//         * MIXTURE, fp 
//         ID=R07393:Bond Cleaved and Formed (2)
//         {C-H:1.0, H-O:1.0}
//
//         */
//        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
//    }
    /**
     * *******************************************************************
     *
     * TEST CASES FOR MCS OPTIMISATION
     *
     * *******************************************************************
     * @throws java.lang.Exception
     *
     * TO DO, optimise MCS Should break phosphate O-P bond not the C-O
     *
     */
//    @Test
//    public void R04459() throws Exception {
//
//        String reactionID = "R04459";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//
//        /*
//         * Expected Solution
//         * MIN, fp ID=R04459:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
//         * BE 665.0, Fragment 1
//         */
//        assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
//    }

    /*
     * Complex Mapping
     *
     * MAX, fp 
     *  ID=R05030:Bond Cleaved and Formed (5)
     *  [C-N:1.0, C-O:1.0, H-N:1.0, H-O:1.0, O-P:2.0]
     * 
     *  BE 1333.0, Fragment 0
     *    
     *
     *
     * @throws Exception
     */
//    @Test
//    public void R05030() throws Exception {
//
//        String reactionID = "R05030";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
//    }

    /*
     * TO DO, mapping NOT SURE, COMPLEX CASE, 
     * Optimise MCS
     * @BUG
     * 
     *
     * @throws Exception
     */
//    @Test
//    public void R06466() throws Exception {
//
//        String reactionID = "R06466";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//
//        /*
//         * Expected Solution
//         * MIN, fp ID=R06466:Bond Cleaved and Formed (5)  C%C:6;  C%O:1;  C-C:5;  C-H:9;  H-O:1; 
//         * BE 4129.0, Fragment 6
//         */
//        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
//    }

    /*
     * MAX, fp 
     * ID=R08855:Bond Cleaved and Formed (4)
     *  [C-N:2.0, C-O:2.0, H-O:2.0, O-P:2.0]
     * 
     *  BE 1996.0, Fragment 0
     *
     *
     * @throws Exception
     */
//    @Test
//    public void R08855() throws Exception {
//
//        String reactionID = "R08855";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
//    }
    /**
     *
     * @throws Exception
     */
//    @Test
//    public void R09688() throws Exception {
//
//        String reactionID = "R09688";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
//
//    }
    /**
     *
     * @throws Exception
     */
//    @Test
//    public void R09709() throws Exception {
//
//        String reactionID = "R09709";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
//
//    }
    /*
     * MIN, fp 
     * ID=R09907:Bond Cleaved and Formed (5)
     * [C%C:4.0, C%O:1.0, C-C:2.0, C-H:5.0, H-O:1.0]
     * 
     * BE 2409.0, Fragment 0
     * 
     *
     * @throws Exception
     */
//    @Test
//    public void R09907() throws Exception {
//        String reactionID = "R09907";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//
//        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
//    }
    /**
     * Linear to RING SYSTEM
     *
     * @throws java.lang.Exception
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     */
//    @Test
//    public void R09909() throws Exception {
//
//        String reactionID = "R09909";
//        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
//    }
    /**
     *
     * @throws Exception
     */
//    @Test
//    public void LinearToRingTest() throws Exception {
//
//        String reactionID = "LinearToRingTest";
//        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN);
//        IPatternFingerprinter formedCleavedWFingerprint = testReactions
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
//
//    }
}
