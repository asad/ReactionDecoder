package com.bioinceptionlabs.generic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;
import org.junit.Test;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class StandardizeReactionIdentityTest {

    private final SmilesParser smilesParser
            = new SmilesParser(SilentChemObjectBuilder.getInstance());

    @Test
    public void duplicateComponentsKeepOccurrenceIdentityWhenStandardized() throws Exception {
        IAtomContainer reactantA = smilesParser.parseSmiles("CC");
        IAtomContainer reactantB = smilesParser.parseSmiles("CC");
        IAtomContainer product = smilesParser.parseSmiles("CCCC");

        markAtoms(reactantA, "benchmarkAtomId", "R:0:");
        markAtoms(reactantB, "benchmarkAtomId", "R:1:");
        markAtoms(product, "benchmarkAtomId", "P:0:");

        IReaction reaction = new Reaction();
        reaction.addReactant(reactantA, 1.0);
        reaction.addReactant(reactantB, 1.0);
        reaction.addProduct(product, 1.0);

        IReaction standardized = new StandardizeReaction().standardize(reaction);

        assertEquals(2, standardized.getReactantCount());
        assertEquals(1.0, standardized.getReactantCoefficient(
                standardized.getReactants().getAtomContainer(0)), 0.0);
        assertEquals(1.0, standardized.getReactantCoefficient(
                standardized.getReactants().getAtomContainer(1)), 0.0);

        IAtomContainer standardizedA = standardized.getReactants().getAtomContainer(0);
        IAtomContainer standardizedB = standardized.getReactants().getAtomContainer(1);

        assertEquals("R:0", standardizedA.getProperty(StandardizeReaction.SOURCE_OCCURRENCE_ID));
        assertEquals("R:1", standardizedB.getProperty(StandardizeReaction.SOURCE_OCCURRENCE_ID));
        assertEquals(Boolean.TRUE, standardizedA.getProperty(StandardizeReaction.PRESERVE_OCCURRENCE_IDENTITY));
        assertEquals(Boolean.TRUE, standardizedB.getProperty(StandardizeReaction.PRESERVE_OCCURRENCE_IDENTITY));

        assertNotNull(standardizedA.getAtom(0).getProperty(StandardizeReaction.SOURCE_ATOM_ID));
        assertNotNull(standardizedB.getAtom(0).getProperty(StandardizeReaction.SOURCE_ATOM_ID));
        assertEquals("R:0:0", standardizedA.getAtom(0).getProperty(StandardizeReaction.SOURCE_ATOM_ID));
        assertEquals("R:1:0", standardizedB.getAtom(0).getProperty(StandardizeReaction.SOURCE_ATOM_ID));
    }

    private void markAtoms(IAtomContainer container, String key, String prefix) {
        for (int atomIndex = 0; atomIndex < container.getAtomCount(); atomIndex++) {
            IAtom atom = container.getAtom(atomIndex);
            atom.setProperty(key, prefix + atomIndex);
        }
    }
}
