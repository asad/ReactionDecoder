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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static mapping.TestUtility.INFORCHEM_RXN;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class AAMTest extends TestUtility {
    private static final Logger LOG = Logger.getLogger(AAMTest.class.getName());

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String RXN_DIR = INFORCHEM_RXN + "rxn" + File.separator + "all" + File.separator;
        /*
         Call AAM Class
         */
        AAMTest aamTest = new AAMTest(RXN_DIR);

    }

    /**
     *
     * @param reactionSet
     * @param remap override existing mappings
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public static List<IReaction> mapReactions(IReactionSet reactionSet, boolean remap) throws FileNotFoundException, Exception {
        List<IReaction> mappedReactionList = new ArrayList<>();
        for (IReaction cdkReaction : reactionSet.reactions()) {

            IReaction mappedReaction = mapReaction(cdkReaction, remap);
            /*
            Add mapped reaction to the list
            */ mappedReactionList.add(mappedReaction);
        }
        return mappedReactionList;
    }

    /**
     *
     * @param cdkReaction reaction for be mapped
     * @param remap override existing mappings
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public static IReaction mapReaction(IReaction cdkReaction, boolean remap) throws FileNotFoundException, Exception {

        String reactionName = cdkReaction.getID();
        IReaction cleanReaction = cleanReaction(cdkReaction, reactionName);
        /*
        * RMT for the reaction mapping
        */
        ReactionMechanismTool rmt = new ReactionMechanismTool(cleanReaction, remap, true, false, new StandardizeReaction());

        /*
        Reaction with hydrogens mapped but unchanged hydrogens suppressed
        */
        //IReaction reactionWithCompressUnChangedHydrogens = rmt.getSelectedSolution().getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        /*
        Reaction with hydrogens mapped
        */
        IReaction mappedReaction = rmt.getSelectedSolution().getReaction();

        /*
        optional step: Renumber the atoms as per mapping
        */
        renumberMappingIDs(mappedReaction);

        return mappedReaction;
    }

    /**
     *
     * @param reaction
     * @param reactionName
     * @return
     * @throws FileNotFoundException
     */
    protected static IReaction cleanReaction(IReaction reaction, String reactionName) throws Exception {
        //write code to fix reactions (Atom type , hydrogens etc.)
        //TO DO
        reaction.setID(reactionName);
        return reaction;
    }

    /**
     *
     * @param RXN_DIR Directory with RXN Files
     */
    public AAMTest(String RXN_DIR) {
        System.out.println("RXN File Directory: " + RXN_DIR);
        /*
        Instance of SMILES with AAM
        */
        SmilesGenerator smilesAAM = SmilesGenerator.generic().withAtomClasses();
        File dir = new File(RXN_DIR);
        File[] files = dir.listFiles();
        for (File file : files) {
            //System.out.println("Iterating: " + file);
            IReaction readReaction;
            try {
                readReaction = readReactionFile(file.getName().split("\\.")[0], RXN_DIR, true, false);
                IReaction mapReaction = mapReaction(readReaction, true);
                System.out.println(" Mapped Reaction SMILES for Reaction " + file.getName() + ": "
                        + smilesAAM.createReactionSMILES(mapReaction));
            } catch (CDKException ex) {
                Logger.getLogger(AAMTest.class.getName()).log(Level.SEVERE, null, ex);
            } catch (Exception ex) {
                Logger.getLogger(AAMTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

}
