/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.graph.invariant.InChINumbersTools;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Canonically labels (permutes) an atom container according to the SMILES
 * canonicalization algorithm.
 *
 * @author maclean
 *
 */
public class InChiMoleculeLabeller implements ICanonicalMoleculeLabeller {

    @Override
    public IAtomContainer getCanonicalMolecule(IAtomContainer container) {
        int[] canonicalPermutation = getCanonicalPermutation(container);
        IAtomContainer permute = AtomContainerAtomPermutor.permute(
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
            labels = InChINumbersTools.getUSmilesNumbers(container);
        } catch (CDKException ex) {
            labels = Canon.label(container, GraphUtil.toAdjList(container));
            Logger.getLogger(InChiMoleculeLabeller.class.getName()).log(Level.SEVERE, null, ex);
        }
        int[] permute = new int[labels.length];
        for (int i = 0; i < labels.length; i++) {
            permute[i] = (int) labels[i] - 1;
        }
        return permute;

    }

    public void getCanonicalPermutation() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    private static final Logger LOG = Logger.getLogger(InChiMoleculeLabeller.class.getName());
}
