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
package aamtool.test;

import java.io.FileNotFoundException;
import junit.framework.Assert;
import static mapping.BaseTest.KEGG_RXN_DIR;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.mechanism.helper.Utility;
import uk.ac.ebi.reactionblast.tools.ReactionSimilarityTool;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MappingTest extends BaseTest {

    @Before
    public void setup() {
        System.out.println("USING : ECBLAST TO MAP FROM RXN DIR: " + KEGG_RXN_DIR);
    }

    /**
     * *******************************************************************
     *
     * TEST CASES FOR WITH BUGS TO DO, mapping NOT SURE, COMPLEX CASE
     *
     * @ BUGS
     * *******************************************************************
     */
    /**
     * @throws java.lang.Exception
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     *
     * @BUG Needs V2000 files system Time taking
     */
    @Test
    public void R06458() throws Exception {
        setup();
        String reactionID = "R06458";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R00959:Bond Cleaved and Formed (2)  H-O:2; O-P:2; 
         * BE 670.0, Fragment 2
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Linear to RING SYSTEM
     *
     * @throws java.lang.Exception
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     */
    @Test
    public void R09909() throws Exception {
        setup();
        String reactionID = "R09909";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }


    /*
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     * Ring rearrangment case
     */
    @Test
    public void R01081() throws Exception {
        setup();
        String reactionID = "R01081";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R01081:Bond Cleaved and Formed (3)  C%O:2;  C-H:2;  H-O:2; 
         * BE 706.0, Fragment 0
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     *
     * It should cleave the bonds between the ring systems (bridge)
     *
     * MIN, fp ID=R01557:Bond Cleaved and Formed (2) [C-O:2.0, H-O:2.0]
     *
     * BE 716.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R01557() throws Exception {
        setup();
        String reactionID = "R01557";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @BUG TO DO, NOT SURE ABOUT THE MAPPINGS, COMPLEX CASE.
     *
     *
     * MIN, fp ID=R02007:Bond Cleaved and Formed (3) [C%C:2.0, C-C:2.0, C-H:4.0]
     *
     * BE 1374.0, Fragment 0
     *
     *
     * @throws Exception
     */
    @Test
    public void R02007() throws Exception {
        setup();
        String reactionID = "R02007";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }


    /*
     TO DO, mapping NOT SURE, COMPLEX CASE, 
     Optimise MCS
     * @BUG
     */
    @Test
    public void R06466() throws Exception {
        setup();
        String reactionID = "R06466";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);

        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R06466:Bond Cleaved and Formed (5)  C%C:6;  C%O:1;  C-C:5;  C-H:9;  H-O:1; 
         * BE 4129.0, Fragment 6
         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * *******************************************************************
     *
     * TEST CASES FOR MCS OPTIMISATION
     *
     * *******************************************************************
     * @throws java.lang.Exception
     */
    /*
     * TO DO, optimise MCS
     * Should break phosphate o-p bond not the c-o
     *
     */
    @Test
    public void R04459() throws Exception {
        setup();
        String reactionID = "R04459";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R04459:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
         * BE 665.0, Fragment 1
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * *******************************************************************
     *
     * SUCESSFUL TEST CASES
     *
     * *******************************************************************
     *
     */
    /*
     *@Ring Breaking
     *
     * MIN, fp 
     ID=R01819:Bond Cleaved and Formed (2)
     [C%O:2.0, C-H:2.0]

     BE 706.0, Fragment 0
     *@throws java.lang.Exception
     */
    @Test
    public void R01819() throws Exception {
        setup();
        String reactionID = "R01819";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *
     *
     * MAX, fp 
     * ID=R09087:Bond Cleaved and Formed (2)
     * [H-O:2.0, O-P:2.0]
     * 
     * BE 670.0, Fragment 0
     *
     * @throws java.lang.Exception
     */
    @Test
    public void R09087() throws Exception {
        setup();
        String reactionID = "R09087";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
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
        setup();
        String reactionID = "R03200";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * 
     * RINGS, fp 
     * ID=R05219:Bond Cleaved and Formed (4)
     * [C-C:2.0, C-O:1.0, C-S:1.0, H-O:2.0]
     * 
     * BE 1322.0, Fragment 0
     * 
     * @throws java.lang.Exception
     */
    @Test
    public void R05219() throws Exception {
        setup();
        String reactionID = "R05219";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Deadlock Test
     *
     *  MIN
     *  ID=R00114:Bond Cleaved and Formed (5)
     *  [C-H:2.0, C-N:2.0, C-O:1.0, C=O:1.0, H-O:1.0]
     * 
     *  BE 1767.0, Fragment 0
     */
    @Test
    public void R00114() throws Exception {
        setup();
        String reactionID = "R00114";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R05217:Bond Cleaved and Formed (5)
     *
     * [C%O:1.0, C-H:1.0, C-O:1.0, H-O:4.0, O=O:1.0]
     *
     * BE 1205.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R05217() throws Exception {
        setup();
        String reactionID = "R05217";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @Long Job
     *
     * MIN, fp ID=R00097:Bond Cleaved and Formed (3) [C-H:1.0, Co-O:2.0,
     * H-O:4.0]
     *
     * BE 0.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00097() throws Exception {
        setup();
        String reactionID = "R00097";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @Flips
     * NAP and NAPH mapping is important, check if P flipping takes place
     */
    @Test
    public void R00023() throws Exception {
        setup();
        String reactionID = "R00023";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @ Carbon
     * NAP and NAPH mapping is important, check if P flipping takes place
     * MIN, fp 
     *ID=R00014:Bond Cleaved and Formed (3)
     *[C-C:2.0, C-H:2.0, H-O:2.0]
     *
     *BE 692.0, Fragment 0
     */
    @Test
    public void R00014() throws Exception {
        setup();
        String reactionID = "R00014";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Complex Mapping
     *
     * MAX, fp 
     *  ID=R05030:Bond Cleaved and Formed (5)
     *  [C-N:1.0, C-O:1.0, H-N:1.0, H-O:1.0, O-P:2.0]
     * 
     *  BE 1333.0, Fragment 0
     *    
     */
    @Test
    public void R05030() throws Exception {
        setup();
        String reactionID = "R05030";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @
     * Phosphate and Water mapping is important, check if P flipping takes place
     * MIN, fp 
     *  ID=R00004:Bond Cleaved and Formed (2)
     *  [H-O:2.0, O-P:2.0]
     * 
     BE 670.0, Fragment 0
     */
    @Test
    public void TestR00004_time() throws Exception {
        setup();
        ReactionMechanismTool testReactions = testReactions("R00004", KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *  MAX, fp 
     *  ID=R00025:Bond Cleaved and Formed (5)
     *  [C=N:1.0, C=O:1.0, H-N:2.0, H-O:3.0, O=O:1.0]
     * 
     *  BE 1908.0, Fragment 0
     */
    @Test
    public void R00025() throws Exception {
        setup();
        String reactionID = "R00025";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }
    /*
     * MAX, fp 
     * ID=R08855:Bond Cleaved and Formed (4)
     *  [C-N:2.0, C-O:2.0, H-O:2.0, O-P:2.0]
     * 
     *  BE 1996.0, Fragment 0
     */

    @Test
    public void R08855() throws Exception {
        setup();
        String reactionID = "R08855";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R03627:Bond Cleaved and Formed (1)
     *  [C-H:2.0]
     * 
     *  BE 0.0, Fragment 0
     */
    @Test
    public void R03627() throws Exception {
        setup();
        String reactionID = "R03627";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * @Ring Breaking and Assimilation
     *
     * MIN, fp ID=R08639:Bond Cleaved and Formed (2) [H-O:2.0, O-P:2.0]
     *
     * BE 670.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R08639() throws Exception {
        setup();
        ReactionMechanismTool testReactions = testReactions("R08639", KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @Ring Breaking and Assimilation
     *
     * MIN, fp ID=R01715:Bond Cleaved and Formed (2) [C-C:1.0, C-O:1.0]
     *
     * BE 704.0, Fragment 0
     *
     * 2 stereo changes, instead of 3
     *
     * @throws Exception
     */
    @Test
    public void R01715() throws Exception {
        setup();
        String reactionID = "R01715";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        IPatternFingerprinter stereoFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        Assert.assertEquals(1, stereoFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R05137:Bond Cleaved and Formed (2)
     *  [C%O:1.0, H-O:2.0]
     * 
     *  BE 353.0, Fragment 0
     * Ring match change missing in the reactant side
     */
    @Test
    public void R05137() throws Exception {
        setup();
        String reactionID = "R05137";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *   MIN, fp 
     *  ID=R05137:Bond Cleaved and Formed (2)
     *   [C%O:1.0, H-O:2.0]
     * 
     *  BE 353.0, Fragment 0
     * 
     */
    @Test
    public void R03673() throws Exception {
        setup();
        String reactionID = "R05137";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R00045:Bond Cleaved and Formed (2)
     *  [H-O:8.0, O=O:1.0]
     * 
     *  BE 494.0, Fragment 0
     * keuekal bond changes mapping misplaced
     */
    @Test
    public void R00045() throws Exception {
        setup();
        String reactionID = "R00045";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        IPatternFingerprinter orderChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        Assert.assertEquals(2, orderChangesWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R01068:Bond Cleaved and Formed (4)
     *  [C%C:1.0, C%O:1.0, C-H:1.0, H-O:3.0]
     * 
     *  BE 694.0, Fragment 0
     */
    @Test
    public void R01068() throws Exception {
        setup();
        String reactionID = "R01068";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R03959:Bond Cleaved and Formed (3)
     *  [C%O:1.0, C-H:1.0, H-O:1.0]
     * 
     *  BE 353.0, Fragment 0
     */
    @Test
    public void R03959() throws Exception {
        setup();
        String reactionID = "R03959";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     ID=R06989:Bond Cleaved and Formed (3)
     [C%O:1.0, C-H:1.0, H-O:1.0]

     BE 353.0, Fragment 0
     */
    @Test
    public void R06989() throws Exception {
        setup();
        String reactionID = "R06989";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }
    /*
     * MIN, fp 
     * ID=R01432:Bond Cleaved and Formed (2)
     * [C%O:2.0, C-H:2.0]
     * 
     *  BE 706.0, Fragment 0
     */

    @Test
    public void R01432() throws Exception {
        setup();
        String reactionID = "R01432";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R00307:Bond Cleaved and Formed (3) [C%O:2.0, C-H:2.0, H-O:2.0]
     *
     * BE 706.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00307() throws Exception {
        setup();
        String reactionID = "R00307";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R03063:Bond Cleaved and Formed (3)
     *  [C-H:2.0, C-N:2.0, C=O:2.0]
     * 
     *   BE 2208.0, Fragment 0
     */
    @Test
    public void R03063() throws Exception {
        setup();
        String reactionID = "R03063";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *  
     *   MIN, fp 
     *   ID=R01015:Bond Cleaved and Formed (2)
     *   [C-H:2.0, H-O:2.0]
     * 
     *   BE 0.0, Fragment 0
     */
    @Test
    public void R01015() throws Exception {
        setup();
        String reactionID = "R01015";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * IMP reaction to test ring re-arrangement Unmapped atoms present
     *
     * MIN, fp ID=R03165:Bond Cleaved and Formed (5) [C%C:2.0, C-C:1.0, C-H:1.0,
     * C-O:1.0, H-O:1.0]
     *
     * BE 1386.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R03165() throws Exception {
        ReactionMechanismTool testReactions = testReactions("R03165", KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R00344:Bond Cleaved and Formed (5) [C-C:1.0, C-H:1.0, C-O:1.0,
     * H-O:2.0, O-P:2.0]
     *
     * BE 1374.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R00344() throws Exception {
        setup();
        String reactionID = "R00344";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Ring exchange test IMP test
     *
     * MIN, fp ID=R04333:Bond Cleaved and Formed (2) [C-O:2.0, H-O:2.0]
     *
     * BE 716.0, Fragment 0
     *
     * @throws Exception
     */
    @Test
    public void R04333() throws Exception {
        setup();
        String reactionID = "R04333";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     * ID=R01717:Bond Cleaved and Formed (1)
     *  [C-O:2.0]
     * 
     *  BE 716.0, Fragment 0
     */
    @Test
    public void R01717() throws Exception {
        setup();
        String reactionID = "R01717";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *   ID=R04165:Bond Cleaved and Formed (2)
     *   [C-H:2.0, H-O:2.0]
     * 
     *   BE 0.0, Fragment 0
     */
    @Test
    public void R04165() throws Exception {
        setup();
        String reactionID = "R04165";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * MIN, fp 
     *  ID=R01394:Bond Cleaved and Formed (2)
     *  [C-H:2.0, H-O:2.0]
     * 
     *  BE 0.0, Fragment 0
     */
    @Test
    public void R01394() throws Exception {
        setup();
        String reactionID = "R01394";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     *  MIN, fp 
     *  ID=R05069:Bond Cleaved and Formed (2)
     *  [C-C:2.0, H-O:2.0]
     *  BE 692.0, Fragment 0
     */
    @Test
    public void R05069() throws Exception {
        setup();
        String reactionID = "R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         *  MIXTURE, fp 
         *  ID=R05069:Bond Cleaved and Formed (2)
         *   {C-C:2; H-O:2;}
         * 
         *   BE 692.0, Fragment 0
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     MIN, fp 
     ID=R07322:Bond Cleaved and Formed (2)
     [C%C:5.0, C-H:2.0]

     BE 1705.0, Fragment 0
     */
    @Test
    public void R07322() throws Exception {
        setup();
        String reactionID = "R07322";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * a single C-N bond and form a single C-C bond expected
     */
    @Test
    public void R03020() throws Exception {
        setup();
        String reactionID = "R03020";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03020:Bond Cleaved and Formed (4)  C-C:1;  C-H:1;  C-N:1;  H-N:1; 
         * BE 651.0, Fragment 2
         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * water phophate bond interaction prefered
     */
    @Test
    public void R03332() throws Exception {
        setup();
        String reactionID = "R03332";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03332:Bond Cleaved and Formed (2)  H-O:2;  O-P:2; 
         * BE 670.0, Fragment 2
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Stereo (R/S)
     */
    @Test
    public void R01903() throws Exception {
        setup();
        String reactionID = "R01903";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R01903:Bond Cleaved and Formed (3)  C%O:1;  C-H:2;  H-O:1; 
         * BE 353.0, Fragment 0
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Stereo (E/Z)
     */
    @Test
    public void R04538() throws Exception {
        setup();
        String reactionID = "R04538";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter fp = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();
        /*
         * Expected Solution
         * MIN, fp ID=R04538:Bond Cleaved and Formed (0)
         * ID=R04538:Bond Stereo Change (1)  C(E/Z):2; 
         * BE 0.0, Fragment 0
         */
        Assert.assertEquals(1, fp.getFeatureCount());
    }

    /*
     * O-P bond changes with water then ATP-ADP prefered big small with water R03187
     */
    @Test
    public void R03187() throws Exception {
        setup();
        String reactionID = "R03187";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03187:Bond Cleaved and Formed (5)  C%N:1;  C-O:1;  H-N:1;  H-O:3;  O-P:2; 
         * BE 1328.0, Fragment 3
         */
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Test for Non Ring RC should differ from RING RC although atoms are same
     * but aromatic and non aromatic bonds should differ
     *
     * @throws Exception
     */
    @Test
    public void NonRingRC() throws Exception {
        setup();
        String reactionID = "R02718";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R02718:Bond Cleaved and Formed (3)  C-O:2;  H-O:2;  O-P:2; 
         * BE 1386.0, Fragment 4
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Test for Non Ring RC should differ from RING RC although atoms are same
     * but aromatic and non aromatic bonds should differ
     *
     * @throws Exception
     */
    @Test
    public void RingRC() throws Exception {
        setup();
        String reactionID = "R01561";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R01561:Bond Cleaved and Formed (4)  C-N:1;  C-O:1;  H-N:1;  H-O:1; 
         * BE 663.0, Fragment 2
         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * MIN, fp ID=R08761:Bond Cleaved and Formed (4) [C-H:2.0, C-O:1.0, H-O:3.0,
     * O=O:1.0]
     *
     * BE 852.0, Fragment 0
     */
    @Test
    public void R08761() throws Exception {
        setup();
        String reactionID = "R08761";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R08765() throws Exception {
        setup();
        String reactionID = "R08765";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R08765:Bond Cleaved and Formed (1)  C%C:2; 
         * BE 682.0, Fragment 0
         */
        Assert.assertEquals(1, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R00011() throws Exception {
        setup();
        String reactionID = "R00011";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R00011:Bond Cleaved and Formed (2)  H-O:2;  O-O:1; 
         * BE 142.0, Fragment 1
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * Image generation failed
     *
     * @throws Exception
     */
    @Test
    public void R03364() throws Exception {
        setup();
        String reactionID = "R03364";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R03364:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
         * BE 665.0, Fragment 1
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     MAX, fp 
     ID=R04558:Bond Cleaved and Formed (6)
     [C-H:1.0, C-N:3.0, C-O:1.0, C=O:1.0, C@N:1.0, H-O:1.0]

     BE 2377.0, Fragment 0
     */
    @Test
    public void R04558() throws Exception {
        setup();
        String reactionID = "R04558";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(6, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Should break phosphate C-N with O-C in RP00084 C00026_C00217 main
     MIXTURE, fp 
     ID=R01148:Bond Cleaved and Formed (3)
     [C-H:2.0, C-N:2.0, C=O:2.0]

     BE 2208.0, Fragment 0
     */
    @Test
    public void R01148() throws Exception {
        setup();
        String reactionID = "R01148";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Expected bond changes
     *
     * R05645	C%O	1 
     * R05645	C(R/S)	10 
     * R05645	C-H	2 
     * R05645	C-O*C=O	1 
     * R05645	H-O	2
     */
    @Test
    public void R05645() throws Exception {
        setup();
        String reactionID = "R05645";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R05645:Bond Cleaved and Formed (3)  C%O:1;  C-H:2;  H-O:2; 
         * BE 353.0, Fragment 0
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     *
     * @throws Exception
     */
    @Test
    public void R06361() throws Exception {
        setup();
        String reactionID = "R06361";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R06361:Bond Cleaved and Formed (5)  C-N:2;  C-O:1;  C-S:1;  H-O:1;  H-S:1; 
         * BE 1240.0, Fragment 4
         */
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    /**
     * IMP N -> O exchange
     *
     * @throws Exception
     */
    @Test
    public void R00093() throws Exception {
        setup();
        String reactionID = "R00093";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R00093:Bond Cleaved and Formed (5)  C-H:2;  C-N:2;  C-O:1;  C=O:1;  H-O:1; 
         * BE 1767.0, Fragment 4 
         */
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void R07635() throws Exception {
        setup();
        String reactionID = "R07635";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp ID=R07635:Bond Cleaved and Formed (3)  H-O:2;  O%P:1;  O-P:1; 
         * BE 665.0, Fragment 
         */
        Assert.assertEquals(3, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     System.out.println("CASE: Condition 11");
     */
    @Test
    public void R05071() throws Exception {
        setup();
        String reactionID = "R05071";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIXTURE, fp ID=R05071:Bond Cleaved and Formed (2)  C-C:2;  H-O:2; 
         * BE 692.0, Fragment 2
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * Water on product should cleave oxygen attached to the ring
     */
    @Test
    public void R02707() throws Exception {
        setup();
        String reactionID = "R02707";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R02707:Bond Cleaved and Formed (2)  C-O:2;  H-O:2; 
         * BE 716.0, Fragment 2
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * phophate should be cleaved - attached to the ring
     */
    @Test
    public void R00959() throws Exception {
        setup();
        String reactionID = "R00959";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MAX, fp ID=R00959:Bond Cleaved and Formed (2)  H-O:2;  O-P:2; 
         * BE 670.0, Fragment 2
         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @Phophate should be cleaved
     * Expected Solution
     * MAX, fp ID=R00959:Bond Cleaved and Formed (2)  H-O:2; O-P:2; 
     * BE 670.0, Fragment 2
     */
    @Test
    public void R01518() throws Exception {
        setup();
        String reactionID = "R01518";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @25405	5.4.99.3
     *  
     * MIXTURE, fp 
     * ID=25405:Bond Cleaved and Formed (2)
     * [C-C:2.0, H-O:2.0]

     * BE 692.0, Fragment 0
     *   
     */
    @Test
    public void Rhea_25405() throws Exception {
        setup();
        String reactionID = "25405";
        ReactionMechanismTool testReactions = testReactions(reactionID, Rhea_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     MAX, fp 
     ID=R07965:Bond Cleaved and Formed (6)
     [C-H:2.0, C-N:1.0, C=O:1.0, H-N:1.0, H-O:2.0, O=O:1.0]

     BE 1598.0, Fragment 0
     * Phosphate mappings
     */
    @Test
    public void R07965() throws Exception {
        setup();
        String reactionID = "R07965";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(6, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     MIN, fp 
     ID=R01569:Bond Cleaved and Formed (2)
     [H-O:2.0, O-P:2.0]

     BE 670.0, Fragment 0
     * Phosphate mappings
     */
    @Test
    public void R01569() throws Exception {
        setup();
        String reactionID = "R01569";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        /*
         Check for the presence of O-P bond
         */
        if (!formedCleavedWFingerprint.getFeatures().contains(new Feature("O-P"))) {
            Assert.fail();
        }
    }

    /*
     MIN, fp 
     ID=R01569:Bond Cleaved and Formed (2)
     [H-O:2.0, O-P:2.0]

     BE 670.0, Fragment 0
     * Phosphate mappings
     */
    @Test
    public void R02918() throws Exception {
        setup();
        String reactionID = "R02918";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
        /*
         Check for the presence of O-P bond
         */
        if (!formedCleavedWFingerprint.getFeatures().contains(new Feature("O-P"))) {
            Assert.fail();
        }
    }

    /*
     MIN, fp 
     ID=R02555:Bond Cleaved and Formed (3)
     [C%O:1.0, C-H:2.0, H-O:1.0]

     BE 353.0, Fragment 0
     */
    @Test
    public void R02555() throws Exception {
        setup();
        String reactionID = "R02555";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void R09688() throws Exception {
        setup();
        String reactionID = "R09688";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void LinearToRingTest() throws Exception {
        setup();
        String reactionID = "LinearToRingTest";
        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R09709() throws Exception {
        setup();
        String reactionID = "R09709";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R00630() throws Exception {
        setup();
        String reactionID = "R00630";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R03979() throws Exception {
        setup();
        String reactionID = "R03979";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R05421() throws Exception {
        setup();
        String reactionID = "R05421";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R10106() throws Exception {
        setup();
        String reactionID = "R10106";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void R01145() throws Exception {
        setup();
        String reactionID = "R01145";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());

    }

    @Test
    public void aromaticity() throws Exception {
        setup();
        String reactionID = "aromaticity_check";
        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());

    }


    /*
     * rxn
     */
    @Test
    public void R02161() throws Exception {
        setup();
        String reactionID = "R02161";
        testReactions(reactionID, KEGG_RXN_DIR);
    }

    /*
     * rxn
     */
    @Test
    public void R00090() throws Exception {
        setup();
        String reactionID = "R00090";
        testReactions(reactionID, KEGG_RXN_DIR);
    }


    /*
     * rxn
     */
    @Test
    public void R00090FP() throws Exception {
        setup();
        String reactionID = "R00090";
        IReaction cdkReaction = readReaction(reactionID, KEGG_RXN_DIR, false);
        for (IAtomContainer ac : cdkReaction.getReactants().atomContainers()) {
            int index = 1;
            for (IAtom a : ac.atoms()) {
                a.setID(String.valueOf(index));
                index++;
            }
            IAtom top = ac.getFirstAtom();
            String circularSMILES = Utility.getCircularSMILES(ac, top, 0, true);
            System.out.println("O: " + circularSMILES);
        }
    }

    @Test
    public void R03627FP() throws Exception {
        setup();
        String reactionID = "R03627";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    @Test
    public void R01194FP() throws Exception {
        setup();
        String reactionID = "R01194";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    @Test
    public void R01188FP() throws Exception {
        setup();
        String reactionID = "R01188";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    @Test
    public void R01123MissingRC() throws FileNotFoundException, Exception {
        setup();
        String reactionID = "R01123";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    /*
     * R01529: Bond change should just report stereo changes
     * 
     */
    @Test
    public void R01529() throws FileNotFoundException, Exception {
        setup();
        String reactionID = "R01529";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    @Test
    public void R01467() throws FileNotFoundException, Exception {
        setup();
        String reactionID = "R01467";
        testRCReactions(reactionID, KEGG_RXN_DIR);
    }

    @Test
    public void R204() throws FileNotFoundException, Exception {
        setup();
        String reactionID = "204";
        testRCReactions(reactionID, Brenda_RXN_DIR);
    }

    @Test
    public void SIM() throws Exception {
        setup();
        String reactionID1 = "R01194";
        BondChangeCalculator testRCReactions1 = map(reactionID1, KEGG_RXN_DIR);
        String reactionID2 = "R01188";
        BondChangeCalculator testRCReactions2 = map(reactionID2, KEGG_RXN_DIR);

        IPatternFingerprinter fp1 = new PatternFingerprinter();
        fp1.add(testRCReactions1.getFormedCleavedWFingerprint());
        fp1.add(testRCReactions1.getOrderChangesWFingerprint());
        fp1.add(testRCReactions1.getStereoChangesWFingerprint());

        IPatternFingerprinter fp2 = new PatternFingerprinter();
        fp2.add(testRCReactions2.getFormedCleavedWFingerprint());
        fp2.add(testRCReactions2.getOrderChangesWFingerprint());
        fp2.add(testRCReactions2.getStereoChangesWFingerprint());

        double similarityBC = ReactionSimilarityTool.getSimilarity(fp1, fp2);
        double similarityRC = ReactionSimilarityTool.getSimilarity(testRCReactions1.getReactionCenterWFingerprint(), testRCReactions2.getReactionCenterWFingerprint());
        System.out.println("RC SIM: " + similarityRC);
        System.out.println("BC SIM: " + similarityBC);
    }

    /*
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R07393() throws Exception {
        String reactionID = "R07393";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIXTURE, fp 
         ID=R07393:Bond Cleaved and Formed (2)
         {C-H:1.0, H-O:1.0}

         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R02996() throws Exception {
        String reactionID = "R02996";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        IPatternFingerprinter orderChangedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getOrderChangesWFingerprint();

        IPatternFingerprinter stereoChangesWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getStereoChangesWFingerprint();

        System.out.println("FP FC" + formedCleavedWFingerprint);
        System.out.println("FP OC" + orderChangedWFingerprint);
        System.out.println("FP ST" + stereoChangesWFingerprint);

        /*
         * Expected Solution
         * MIN, fp 
         ID=R02996:Bond Cleaved and Formed (2)
         [H-O:2.0, O-S:2.0]

         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

    /*
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R03775() throws Exception {
        String reactionID = "R03775";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIXTURE, fp 
         ID=R07393:Bond Cleaved and Formed (2)
         {C-H:1.0, H-O:1.0}

         */
        Assert.assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }


    /*
     * @BUG
     * TO DO 
     * Flip phosphate
     */
    @Test
    public void GeneralBug() throws Exception {
        String reactionID = "reaction2";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIN, fp 
         * ID=reaction:Bond Cleaved and Formed (4)
         * [C-S:1.0, C=O:1.0, O=O:1.0, O=S:1.0]         
         * BE 2087.0, Fragment 0

         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

    @Test
    public void R09907() throws Exception {
        String reactionID = "R09907";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, OTHER_RXN);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * MIN, fp 
         * ID=R09907:Bond Cleaved and Formed (5)
         * [C%C:4.0, C%O:1.0, C-C:2.0, C-H:5.0, H-O:1.0]
         * 
         * BE 2409.0, Fragment 0
         * 
         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }

}
