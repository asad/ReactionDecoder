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

import static java.lang.System.out;
import java.util.Iterator;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class AtomContainerAtomPermutor extends Permutor
        implements Iterator<IAtomContainer> {

    private static boolean useA = false;
    private static boolean clone = false;

    /**
     *
     * @param p
     * @param atomContainer
     * @return
     */
    public static IAtomContainer permute(int[] p, IAtomContainer atomContainer) {
//        if (useA) {
//            return permuteA(p, atomContainer);
//        } else if (!clone) {
//            permuteWithoutClone(p, atomContainer);
//            return atomContainer;
//        } else {
//            return permuteB(p, atomContainer);
//        }
        return permuteB(p, atomContainer);
    }

    private static IAtomContainer permuteA(int[] p, IAtomContainer atomContainer) {
        IAtomContainer permutedContainer = null;
        try {
            permutedContainer
                    = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            for (int i = 0; i < p.length; i++) {
                IAtom atom = atomContainer.getAtom(p[i]);
                permutedContainer.addAtom(atom.clone());
            }
            for (IBond bond : atomContainer.bonds()) {
                IBond clonedBond = bond.clone();
                clonedBond.setAtoms(new IAtom[clonedBond.getAtomCount()]);
                int i = 0;
                for (IAtom atom : bond.atoms()) {
                    int index = atomContainer.indexOf(atom);
                    IAtom permutedAtom = permutedContainer.getAtom(p[index]);
                    clonedBond.setAtom(permutedAtom, i++);
                }
                permutedContainer.addBond(clonedBond);
            }

        } catch (CloneNotSupportedException cne) {
            //?
            out.println(cne);
        }

        return permutedContainer;
    }

    private static IAtomContainer permuteB(int[] p, IAtomContainer atomContainer) {
        IAtomContainer permutedContainer = null;
//        System.out.println("permuting " + java.util.Arrays.toString(p));
        try {
            permutedContainer = atomContainer.clone();
            int n = atomContainer.getAtomCount();
            IAtom[] permutedAtoms = new IAtom[n];
            for (int originalIndex = 0; originalIndex < n; originalIndex++) {
                // get the newly cloned atom 
                IAtom atom = permutedContainer.getAtom(originalIndex);

                // permute the index
                int newIndex = p[originalIndex];

                // put the atom in the new place
                permutedAtoms[newIndex] = atom;
            }
            permutedContainer.setAtoms(permutedAtoms);
        } catch (CloneNotSupportedException cne) {
            //?
            out.println(cne);
        }
        return permutedContainer;
    }

    private static void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
        int n = atomContainer.getAtomCount();
        IAtom[] permutedAtoms = new IAtom[n];
        for (int originalIndex = 0; originalIndex < n; originalIndex++) {
            // get the newly cloned atom 
            IAtom atom = atomContainer.getAtom(originalIndex);

            // permute the index
            int newIndex = p[originalIndex];

            // put the atom in the new place
            permutedAtoms[newIndex] = atom;
        }
        atomContainer.setAtoms(permutedAtoms);
    }
    private IAtomContainer original;

    /**
     *
     * @param atomContainer
     */
    public AtomContainerAtomPermutor(IAtomContainer atomContainer) {
        super(atomContainer.getAtomCount());
        original = atomContainer;
    }

    /**
     *
     * @param atomContainer
     * @param useA
     * @param clone
     */
    public AtomContainerAtomPermutor(IAtomContainer atomContainer, boolean useA, boolean clone) {
        this(atomContainer);
        useA = useA;
        clone = clone;
    }

    @Override
    public IAtomContainer next() {
        int[] p = super.getNextPermutation();
        return permute(p, original);
    }

    @Override
    public void remove() {
        // can just increase rank....
    }
}
