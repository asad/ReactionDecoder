/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package uk.ac.bioinception;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
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
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        String reactionSM
                = "O=C(COCc1ccccc1)N1CCCc2sc(-c3ccc(O[C@@H]4C[C@@H]"
                + "(N5CCCCC5)C4)cc3)nc21>>O=C(CO)N1CCCc2sc"
                + "(-c3ccc(O[C@@H]4C[C@@H](N5CCCCC5)C4)cc3)nc21";

//"CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O";
        //                "CC(C)CCCC(C)CCCC(C)CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])"
        //                + "(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])"
        //                + "=O)N1C=NC2=C1N=CN=C2N.CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
        //                + "COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])"
        //                + "([O-])=O)N1C=NC2=C1N=CN=C2N>>CC(C)CCCC(C)CCCC(C)CCC(=O)C(C)C(=O)"
        //                + "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]"
        //                + "([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N.CC(C)(COP([O-])"
        //                + "(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)"
        //                + "N1C=NC2=C1N=CN=C2N)[C@@H](O)C(=O)NCCC(=O)NCCS";
//                "[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)(O[C@@H]1C[C@@](OC[C@H]2O[C@@H](OC[C@H]3O[C@H]"
//                + "(OP([O-])([O-])=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H]"
//                + "(O)CCCCCCCCCCC)[C@@H]3O)[C@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)"
//                + "[C@@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H]2OP([O-])([O-])=O)"
//                + "(O[C@]([H])([C@H](O)CO)[C@@H]1O[C@H]1O[C@H]([C@@H](O)CO)[C@@H](OP([O-])([O-])=O)"
//                + "[C@H](O[C@H]2O[C@H]([C@@H](O)CO[C@H]3O[C@H]([C@@H](O)CO)[C@@H](O)[C@H](O)[C@@H]3O)"
//                + "[C@@H](OP([O-])([O-])=O)[C@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)"
//                + "[C@@H](O)[C@H](O)[C@H]3O)[C@@H]2O)[C@@H]1O)C([O-])=O)C([O-])=O)[C@H](O)CO.OC[C@H]1O[C@H]"
//                + "(OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=CC(=O)NC2=O)[C@H](O)[C@@H](O)"
//                + "[C@@H]1O>>[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)(O[C@@H]1C[C@@](OC[C@H]2O[C@@H](OC[C@H]3O[C@H]"
//                + "(OP([O-])([O-])=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]3O)"
//                + "[C@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)"
//                + "[C@@H]2OP([O-])([O-])=O)(O[C@]([H])([C@H](O)CO)[C@@H]1O[C@H]1O[C@H]([C@@H](O)CO)[C@@H](OP([O-])([O-])=O)"
//                + "[C@H](O[C@H]2O[C@H]([C@@H](O)CO[C@H]3O[C@H]([C@@H](O)CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H](OP([O-])([O-])=O)"
//                + "[C@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H]"
//                + "(O)[C@H]4O)[C@H]3O)[C@@H]2O)[C@@H]1O)C([O-])=O)C([O-])=O)[C@H](O)CO.[H+].O[C@H]1[C@@H](O)[C@@H](O[C@@H]1COP"
//                + "([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O";
        String reactionName = "Test";

        IReaction cdkReaction = smilesParser.parseReactionSmiles(reactionSM);

        IReaction performAtomAtomMapping = performAtomAtomMapping(cdkReaction, reactionName);
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
        boolean complexMapping = true;//Rings
        StandardizeReaction standardizeReaction = new StandardizeReaction(); //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, complexMapping, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();//Fetch the AAM Solution
        IReaction reaction = s.getReaction();//Fetch Mapped Reaction
        return reaction;
    }

}
