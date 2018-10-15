/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.bioinception;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class AAMExample {
    
    private AAMExample() {
    }
    
    public static void main(String[] args) throws CloneNotSupportedException, CDKException, AssertionError, Exception {
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        
        String reactionSM = "CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O";
        String reactionName = "Test";
        
        IReaction cdkReaction = smilesParser.parseReactionSmiles(reactionSM);
        
        IReaction performAtomAtomMapping = AAMExample.performAtomAtomMapping(cdkReaction, reactionName);
        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping));
        //View mapped reaction http://www.simolecule.com/cdkdepict/depict.html
        
    }

    /**
     *
     * @param cdkReaction
     * @param reactionName
     * @return
     * @throws InvalidSmilesException
     * @throws AssertionError
     * @throws Exception
     */
    public static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
        cdkReaction.setID(reactionName);
        /*
         RMT for the reaction mapping
         */
        boolean forceMapping = true;//Overrides any mapping present int the reaction
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers
        StandardizeReaction standardizeReaction = new StandardizeReaction(); //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();//Fetch the AAM Solution
        IReaction reaction = s.getReaction();//Fetch Mapped Reaction
        return reaction;
    }
    
}
