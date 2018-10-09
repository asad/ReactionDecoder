/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad@ebi.ac.uk>.
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
package uk.ac.ebi.aamtool;

import static org.junit.Assert.assertEquals;
import org.junit.Test;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.MappingUtility;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReArrangments extends MappingUtility {

    /*
     ************************
     * COMPLEX CASES
     ************************
     */
    /**
     * Complex case, Re-arrangementsMAX, fp
     * ID=unbalanced_rearrangement_reactions:Bond Cleaved and Formed (2)
     * [C-C:1.0, C-O:1.0]
     *
     * MAX, fp ID=unbalanced_rearrangement_reactions:Bond Order Change (0) []
     *
     * BE 704.0, Fragment 0
     *
     *
     * @throws java.lang.Exception
     */
    @Test
    public void MolecularRearrangments() throws Exception {

        String reactionID = "unbalanced_rearrangement_reactions";
        ReactionMechanismTool testReactions = testReactions(reactionID, BUG_RXN_DIR);
        IPatternFingerprinter formedCleavedWFingerprint = testReactions
                .getSelectedSolution()
                .getBondChangeCalculator()
                .getFormedCleavedWFingerprint();
        assertEquals(2, formedCleavedWFingerprint.getFeatureCount());
    }

}
