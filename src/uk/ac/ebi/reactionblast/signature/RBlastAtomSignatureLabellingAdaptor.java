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
package uk.ac.ebi.reactionblast.signature;

import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.tools.labelling.AtomContainerAtomPermutor;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;

/**
 *
 * @author maclean
 *
 */
public class RBlastAtomSignatureLabellingAdaptor
        implements ICanonicalMoleculeLabeller {
    private static final Logger LOG = Logger.getLogger(RBlastAtomSignatureLabellingAdaptor.class.getName());

    private int atomIndex;
    private int height;

    public RBlastAtomSignatureLabellingAdaptor(int atomIndex) {
        this.atomIndex = atomIndex;
        this.height = -1;
    }

    public RBlastAtomSignatureLabellingAdaptor(int atomIndex, int height) {
        this.atomIndex = atomIndex;
        this.height = height;
    }

    @Override
    public IAtomContainer getCanonicalMolecule(IAtomContainer container) {
        return AtomContainerAtomPermutor.permute(
                getCanonicalPermutation(container), container);
    }

    @Override
    public int[] getCanonicalPermutation(IAtomContainer container) {
        RBlastMoleculeSignature molSig
                = new RBlastMoleculeSignature(container);
        if (height == -1) {
            return molSig.getAtomSignature(atomIndex)
                    .getCanonicalLabelling(container.getAtomCount());
        } else {
            return molSig.getAtomSignature(atomIndex, height)
                    .getCanonicalLabelling(container.getAtomCount());
        }
    }
}
