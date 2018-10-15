/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.tools.labelling;

import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.graph.GraphUtil.toAdjList;
import static org.openscience.cdk.graph.invariant.Canon.label;
import static org.openscience.cdk.graph.invariant.InChINumbersTools.getUSmilesNumbers;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.graphics.direct.DirectReactionDrawer;
import static uk.ac.ebi.reactionblast.tools.labelling.AtomContainerAtomPermutor.permute;

/**
 * Canonically labels (permutes) an atom container according to the SMILES
 * canonicalization algorithm.
 *
 * @author maclean
 *
 */
public class InChiMoleculeLabeller implements ICanonicalMoleculeLabeller {

    private final static ILoggingTool LOGGER
            = createLoggingTool(DirectReactionDrawer.class);

    /**
     *
     * @param container
     * @return
     */
    @Override
    public IAtomContainer getCanonicalMolecule(IAtomContainer container) {
        int[] canonicalPermutation = getCanonicalPermutation(container);
        IAtomContainer permute = permute(
                canonicalPermutation, container);
        if (container.getID() != null) {
            permute.setID(container.getID());
        }

        for (int i : canonicalPermutation) {
            if (container.getAtom(i).getID() != null) {
                permute.getAtom(i).setID(container.getAtom(i).getID());
            }
        }
        return permute;
    }

    /**
     * Given a molecule (possibly disconnected) compute the labels which would
     * order the atoms by increasing canonical labeling.
     *
     * @param container
     * @return the permutation
     */
    @Override
    public int[] getCanonicalPermutation(IAtomContainer container) {
        long[] labels;
        try {
            labels = getUSmilesNumbers(container);
        } catch (CDKException ex) {
            labels = label(container, toAdjList(container));
            LOGGER.error(SEVERE, null, ex);
        }
        int[] permute = new int[labels.length];
        for (int i = 0; i < labels.length; i++) {
            permute[i] = (int) labels[i] - 1;
        }
        return permute;
    }

    /**
     *
     */
    public void getCanonicalPermutation() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
