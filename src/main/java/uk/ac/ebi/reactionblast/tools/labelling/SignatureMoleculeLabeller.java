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

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Arrays.sort;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;

/**
 * Canonically labels (permutes) an atom container according to the signature
 * canonicalization algorithm.
 *
 * @author maclean
 *
 */
public class SignatureMoleculeLabeller implements ICanonicalMoleculeLabeller {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(SignatureMoleculeLabeller.class);

    /**
     *
     * @param container
     * @return
     */
    @Override
    public IAtomContainer getCanonicalMolecule(IAtomContainer container) {
        try {
            IAtomContainer permute = container.clone();

            if (container.getID() != null) {
                permute.setID(container.getID());
            }

            for (int i = 0; i < permute.getAtomCount(); i++) {
                if (container.getAtom(i).getID() != null) {
                    permute.getAtom(i).setID(container.getAtom(i).getID());
                }
            }

            int[] canonicalPermutation = getCanonicalPermutation(permute);

            permuteWithoutClone(canonicalPermutation, permute);

            return permute;
        } catch (CloneNotSupportedException ex) {
            LOGGER.error(SEVERE, null, ex);
        }

        return null;
    }

    /**
     *
     * @param container
     * @return
     */
    @Override
    public int[] getCanonicalPermutation(IAtomContainer container) {
        MoleculeSignature molSig = new MoleculeSignature(container);
        return molSig.getCanonicalLabels();
    }

    /*
    This is a very imp code modified by John May
    The idea is to canonicalise the atoms and bonds
     */
    private void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
        int n = atomContainer.getAtomCount();
        IAtom[] permutedAtoms = new IAtom[n];

        for (int i = 0; i < n; i++) {
            IAtom atom = atomContainer.getAtom(i);
            permutedAtoms[p[i]] = atom;
            atom.setProperty("label", p[i]);
        }
        atomContainer.setAtoms(permutedAtoms);

        IBond[] bonds = getBondArray(atomContainer);
        sort(bonds, (IBond o1, IBond o2) -> {
            int u = o1.getAtom(0).getProperty("label");
            int v = o1.getAtom(1).getProperty("label");
            int x = o2.getAtom(0).getProperty("label");
            int y = o2.getAtom(1).getProperty("label");
            int min1 = min(u, v);
            int min2 = min(x, y);
            int max1 = max(u, v);
            int max2 = max(x, y);

            int minCmp = Integer.compare(min1, min2);
            if (minCmp != 0) {
                return minCmp;
            }
            int maxCmp = Integer.compare(max1, max2);
            if (maxCmp != 0) {
                return maxCmp;
            }
            LOGGER.debug("pokemon!");
            throw new InternalError();
        });
        atomContainer.setBonds(bonds);
    }

}
