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
package mapping;

import junit.framework.Assert;
import static mapping.BaseTest.KEGG_RXN_DIR;
import org.junit.Test;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class AlgorithmTest extends BaseTest {

    /*
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R05069() throws Exception {
        String reactionID = "R05069";
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
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R01068() throws Exception {
        String reactionID = "R01068";//"R05069";
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
     * @BUG
     * TO DO 
     * Functional group regarrangment
     */
    @Test
    public void R03200() throws Exception {
        String reactionID = "R03200";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         * MIXTURE, fp 
         ID=R03200:Bond Cleaved and Formed (5)
         {C%C:5; C%O:1; C-C:2; C-H:5; H-O:1;}

         */
        Assert.assertEquals(5, formedCleavedWFingerprint.getFeatureCount());
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
     * linear to ring
     */
    @Test
    public void R06466() throws Exception {
        String reactionID = "R06466";//"R05069";
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
    public void R07965() throws Exception {
        String reactionID = "R07965";//"R05069";
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
        ReactionMechanismTool testReactions = testReactions(reactionID, general_bug);
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
        ReactionMechanismTool testReactions = testReactions(reactionID, general_bug);
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

    @Test
    public void R09709() throws Exception {
        String reactionID = "R09709";//"R05069";
        ReactionMechanismTool testReactions = testReactions(reactionID, KEGG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();

        /*
         * Expected Solution
         *  
         * MIN, fp 
         * ID=R09907:Bond Cleaved and Formed (5)
         * [C%C:4.0, C%O:1.0, C-C:2.0, C-H:5.0, H-O:1.0]
         * 
         *  BE 2409.0, Fragment 0
         * 
         */
        Assert.assertEquals(4, formedCleavedWFingerprint.getFeatureCount());
    }
}
