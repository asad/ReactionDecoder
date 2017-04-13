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

import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtom;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import static java.util.logging.Logger.getLogger;
import static uk.ac.ebi.reactionblast.tools.ImageGenerator.TopToBottomReactionLayoutImage;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class RDTTest extends MappingUtility {

    private static final Logger LOG = getLogger(RDTTest.class.getName());

    /*
     */
    @Test
    public void TestCompareIndigoMapping() throws Exception {

        String reactionSM = "[Na+].[CH3:2][C:3](=[O:5])[O-:4].[C:6]([O-:14])(=[O:13])[c:7]1[cH:12][cH:11][cH:10][cH:9][cH:8]1.[Na+]>CN1CCCC1=O>[C:3]([O:4][CH2:8][C:7](=[CH2:6])[CH3:12])(=[O:5])[CH3:2].[O-:14][C:6]([c:7]1[cH:12][cH:11][cH:10][cH:9][cH:8]1)=[O:13]";
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = sp.parseReactionSmiles(reactionSM);
        String rid = "TestUSPTOReaction";
        parseReactionSmiles.setID(rid);
        IPatternFingerprinter formedCleavedWFingerprint = new PatternFingerprinter();
        Map<IAtom, IAtom> mappings = ReactionMappingUtility.getMappings(parseReactionSmiles);
        Set<IBond> bondChanges = ReactionMappingUtility.getBondCleavedFormedChanges(parseReactionSmiles, mappings);
        for (IBond bond : bondChanges) {
            try {
                formedCleavedWFingerprint.add(new Feature(ReactionMappingUtility.getCanonicalisedBondChangePattern(bond), 1.0));
                formedCleavedWFingerprint.setFingerprintID(rid);
            } catch (CDKException ex) {
                Logger.getLogger(USPTOTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        ReactionMechanismTool rmt = new ReactionMechanismTool(parseReactionSmiles, true, new StandardizeReaction());
        MappingSolution s = rmt.getSelectedSolution();

        /*
         * Code for decipt Image generation
         */
        IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator().getReaction();
        TopToBottomReactionLayoutImage(reactionWithCompressUnChangedHydrogens, (rid + s.getAlgorithmID()), "Output");

        //System.out.println("F/C RDT " + s.getBondChangeCalculator().getFormedCleavedWFingerprint());
        System.out.println("INDIGO " + formedCleavedWFingerprint.getFeatures() + ", RMT " + s.getBondChangeCalculator().getFormedCleavedWFingerprint().getFeatures());

    }

}
