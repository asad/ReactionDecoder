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
package uk.ac.ebi.reactionblast.signature;

import static java.lang.System.out;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;

/**
 * Finds atoms in a molecule or reaction that have signatures matching one of a
 * set of query signatures.
 *
 * @author maclean
 *
 */
public class SignatureMatcher {

    private int minHeight;

    private int maxHeight;

    /**
     * A signature matcher with min height of 1 and max height of 3.
     */
    public SignatureMatcher() {
        this(1, 3);
    }

    /**
     * A signature matcher with max height of 3.
     *
     * @param minHeight the minimum signature height (inclusive) to check
     */
    public SignatureMatcher(int minHeight) {
        this(minHeight, 3);
    }

    /**
     * A signature matcher that checks signatures between minHeight and
     * maxHeight.
     *
     * @param minHeight
     * @param maxHeight
     */
    public SignatureMatcher(int minHeight, int maxHeight) {
        this.minHeight = minHeight;
        this.maxHeight = maxHeight;
    }

    /**
     * Search a reaction for the list of atoms that are the root atoms for the
     * supplied signatures. To do this, the matcher finds the signature with the
     * greatest height that matches. All atoms in all atom containers in the
     * reaction are checked.
     *
     * @param signatureStrings the signatures to match
     * @param reaction the reaction to search for matching root atoms
     * @return a list of atoms
     */
    public List<IAtom> getMatchingRootAtoms(List<String> signatureStrings, IReaction reaction) {
        List<IAtom> roots = new ArrayList<>();
        for (IAtomContainer atomContainer : getAllAtomContainers(reaction)) {
            getMatchingRootAtoms(roots, signatureStrings, atomContainer);
        }

        return roots;
    }

    /**
     * Search an atom container for the list of atoms that are the root atoms
     * for the supplied signatures. To do this, the matcher finds the signature
     * with the greatest height that matches. All atoms in all atom containers
     * in the reaction are checked.
     *
     * @param signatureStrings the signatures to match
     * @param atomContainer the atom container to search for matching root atoms
     * @return a list of atoms
     */
    public List<IAtom> getMatchingRootAtoms(List<String> signatureStrings, IAtomContainer atomContainer) {
        List<IAtom> roots = new ArrayList<>();
        getMatchingRootAtoms(roots, signatureStrings, atomContainer);
        return roots;
    }

    /**
     * Fill a list with matching root atoms.
     *
     * @param roots
     * @param signatureStrings
     * @param atomContainer
     */
    private void getMatchingRootAtoms(List<IAtom> roots, List<String> signatureStrings, IAtomContainer atomContainer) {
        RBlastMoleculeSignature moleculeSignature = new RBlastMoleculeSignature(atomContainer);
        for (int atomIndex = 0; atomIndex < atomContainer.getAtomCount(); atomIndex++) {
            RBlastAtomSignature matchingAtomSignature = match(atomIndex, signatureStrings, moleculeSignature);
            if (matchingAtomSignature != null) {
                roots.add(atomContainer.getAtom(atomIndex));
            }
        }
    }

    /**
     * Find the highest matching signature.
     *
     * @param atomIndex the atom index
     * @param signatureStrings the list of query signatures
     * @param moleculeSignature the signature object referencing the target mol
     * @return a matching signature or null if none match
     */
    private RBlastAtomSignature match(
            int atomIndex, List<String> signatureStrings, RBlastMoleculeSignature moleculeSignature) {
        for (int height = maxHeight; height >= minHeight; height--) {
            RBlastAtomSignature sigOfHAtI = moleculeSignature.getAtomSignature(atomIndex, height);
            String signatureStringOfHAtI = sigOfHAtI.toCanonicalString();
            for (String querySignatureString : signatureStrings) {
                if (querySignatureString.equals(signatureStringOfHAtI)) {
                    out.println("match of height " + height + " at " + atomIndex + " to " + querySignatureString);
                    return sigOfHAtI;
                }
            }
        }
        return null;
    }

}
