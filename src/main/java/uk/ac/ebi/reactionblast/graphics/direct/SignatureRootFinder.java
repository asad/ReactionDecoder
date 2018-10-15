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
package uk.ac.ebi.reactionblast.graphics.direct;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;
import uk.ac.ebi.reactionblast.mapping.helper.RBlastReaction;
import uk.ac.ebi.reactionblast.signature.SignatureMatcher;

/**
 * Given an IReaction, a set of signatures, and (optionally) a set of bond
 * changes this class finds the set of 'root systems' of connected atoms that
 * are roots of one of the signatures.
 *
 * @author maclean
 *
 */
public class SignatureRootFinder {

    /**
     *
     * @param rblReaction
     * @return
     */
    public static Map<IAtomContainer, List<RootSystem>> findRootSystems(
            RBlastReaction rblReaction) {
        // get all the bond changes
        List<IBond> allBondChanges = new ArrayList<>();
        allBondChanges.addAll(rblReaction.getBondsCleavedInReactant());
        allBondChanges.addAll(rblReaction.getBondsFormedInProduct());
        allBondChanges.addAll(rblReaction.getBondsOrderChangedInReactant());
        allBondChanges.addAll(rblReaction.getBondsOrderChangedInProduct());

        // get all the atom stereo changes
        List<IAtom> atomChanges = new ArrayList<>();
        atomChanges.addAll(rblReaction.getAtomStereoProductMap().keySet());
        atomChanges.addAll(rblReaction.getAtomStereoReactantMap().keySet());

        return findRootSystems(
                rblReaction.getReaction(), allBondChanges, atomChanges);
    }

    /**
     *
     * @param reaction
     * @param bondChanges
     * @param atomChanges
     * @return
     */
    public static Map<IAtomContainer, List<RootSystem>> findRootSystems(
            IReaction reaction, List<IBond> bondChanges, List<IAtom> atomChanges) {

        Map<IAtomContainer, List<RootSystem>> rootSystems
                = new HashMap<>();

        // separate bond and atom changes by atomContainer
        for (IAtomContainer atomContainer : getAllAtomContainers(reaction)) {
            List<IBond> bonds = new ArrayList<>();
            for (IBond bond : bondChanges) {
                if (atomContainer.contains(bond)) {
                    bonds.add(bond);
                }
            }
            List<IAtom> atoms = new ArrayList<>();
            for (IAtom atom : atomChanges) {
                if (atomContainer.contains(atom)) {
                    atoms.add(atom);
                }
            }
            rootSystems.put(atomContainer,
                    findRootSystems(
                            atomContainer, bonds, atoms));
        }
        return rootSystems;
    }

    /**
     *
     * @param atomContainer
     * @param bondChanges
     * @param atomChanges
     * @return
     */
    public static List<RootSystem> findRootSystems(
            IAtomContainer atomContainer, List<IBond> bondChanges, List<IAtom> atomChanges) {

        // each bond has a label, to keep track of which root system it is in
        int bSize = bondChanges.size();
        int[] bondSystemLabels = new int[bSize];
        int maxSystemLabel = 1;
        Stack<RootSystem> rootSystems = new Stack<>();
        for (int bondIndex = 0; bondIndex < bondChanges.size(); bondIndex++) {
            IBond bond = bondChanges.get(bondIndex);

            // check what system the bond is in
            if (bondSystemLabels[bondIndex] == 0) {

                // not part of a root system, assign it to one
                int currentSystemLabel = maxSystemLabel;
                for (int rLabel = 1; rLabel <= rootSystems.size(); rLabel++) {
                    RootSystem rootSystem = rootSystems.get(rLabel - 1);
                    if (adjacent(bond, rootSystem, atomContainer)) {
                        currentSystemLabel = rLabel;
                        rootSystem.addRootsFromBond(bond);
                        break;
                    }
                }

                // assign the system label, and increment the max if necessary
                bondSystemLabels[bondIndex] = currentSystemLabel;
                if (currentSystemLabel == maxSystemLabel) {
                    RootSystem system = new RootSystem();
                    system.addRootsFromBond(bond);
                    rootSystems.add(system);
                    maxSystemLabel++;
                }
            } else {
                // already in a root system

            }
        }

        // now do the same for the stereo atoms
        int aSize = atomChanges.size();
        int[] atomSystemLabels = new int[aSize];
        for (int atomIndex = 0; atomIndex < atomChanges.size(); atomIndex++) {
            IAtom atom = atomChanges.get(atomIndex);

            // check what system the atom is in
            if (atomSystemLabels[atomIndex] == 0) {
                // not part of a root system, assign it to one
                int currentSystemLabel = maxSystemLabel;
                for (int rLabel = 1; rLabel <= rootSystems.size(); rLabel++) {
                    RootSystem rootSystem = rootSystems.get(rLabel - 1);
                    if (adjacent(atom, rootSystem, atomContainer)) {
                        currentSystemLabel = rLabel;
                        rootSystem.addRoot(atom);
                        break;
                    }
                }

                // assign the system label, and increment the max if necessary
                atomSystemLabels[atomIndex] = currentSystemLabel;
                if (currentSystemLabel == maxSystemLabel) {
                    RootSystem system = new RootSystem();
                    system.addRoot(atom);
                    rootSystems.add(system);
                    maxSystemLabel++;
                }

            } else {
                // already in a root system

            }
        }
//        System.out.println("RS for " + atomContainer.getID() + " = " + rootSystems);

        // now, join together any root systems that are adjacent
        boolean merging = true;
        while (merging) {
            Stack<RootSystem> mergedRootSystems = new Stack<>();
            if (rootSystems.isEmpty()) {
                merging = false;
                break;
            }
            RootSystem rootSystem = rootSystems.pop();
            boolean hasMerged = false;
            for (RootSystem otherRootSystem : rootSystems) {
                if (adjacent(rootSystem, otherRootSystem, atomContainer)) {
                    mergedRootSystems.add(rootSystem.merge(otherRootSystem));
                    hasMerged = true;
//                    System.out.println("merge " + rootSystem + " and " + otherRootSystem);
                } else {
                    mergedRootSystems.add(otherRootSystem);
                }
            }
            if (hasMerged) {
                merging = true;
            } else {
                mergedRootSystems.add(rootSystem);
                merging = false;
            }
            rootSystems = mergedRootSystems;
        }

        // finally add in the leaves
        for (RootSystem rs : rootSystems) {
            List<IAtom> roots = rs.getRoots();
            for (IAtom root : roots) {
                for (IAtom leaf : atomContainer.getConnectedAtomsList(root)) {
                    if (!roots.contains(leaf)) {
                        rs.addLeaf(leaf);
                    }
                }
            }
        }

//        System.out.println("RS for " + atomContainer.getID() + " = " + rootSystems + " after merging ");
        return rootSystems;
    }

    private static boolean adjacent(RootSystem rsI, RootSystem rsJ, IAtomContainer atomContainer) {
        for (int idxI = 0; idxI < rsI.getRoots().size(); idxI++) {
            IAtom atomI = rsI.getRoots().get(idxI);
            for (int idxJ = 0; idxJ < rsJ.getRoots().size(); idxJ++) {
                IAtom atomJ = rsJ.getRoots().get(idxJ);
                for (IBond bond : atomContainer.bonds()) {
                    if (bond.contains(atomI) && bond.contains(atomJ)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static boolean adjacent(IBond bond, RootSystem rootSystem, IAtomContainer atomContainer) {
        for (IAtom root : rootSystem.getRoots()) {
            for (IBond connectedBond : atomContainer.getConnectedBondsList(root)) {
                if (bond == connectedBond) {
                    return true;
                }

            }
        }
        return false;
    }

    private static boolean adjacent(IAtom atom, RootSystem rootSystem, IAtomContainer atomContainer) {
        // TODO Auto-generated method stub
        return false;
    }

    /**
     *
     * @param reaction
     * @param signatureStrings
     * @return
     */
    public static List<RootSystem> findRootSystems(
            IReaction reaction, List<String> signatureStrings) {
        List<RootSystem> rootSystems = new ArrayList<>();

        // find the root atoms for each container, and connect them
        SignatureMatcher matcher = new SignatureMatcher();
        for (IAtomContainer atomContainer : getAllAtomContainers(reaction)) {
            List<IAtom> roots
                    = matcher.getMatchingRootAtoms(signatureStrings, atomContainer);
            rootSystems.addAll(find(atomContainer, roots));
        }
        return rootSystems;
    }

    private static List<RootSystem> find(IAtomContainer atomContainer,
            List<IAtom> roots) {
        List<RootSystem> rootSystems = new ArrayList<>();

        // the root system labels
        int currentLabel = 1;
        int[] labels = new int[atomContainer.getAtomCount()];

        for (IAtom root : roots) {
            List<IAtom> component = new ArrayList<>();
            dfs(
                    null, root, currentLabel, labels, atomContainer, roots, component);

            // for non-empty components, add all the atoms as roots
            if (component.size() > 0) {
                RootSystem rootSystem = new RootSystem();
                for (IAtom rootAtom : component) {
                    rootSystem.addRoot(rootAtom);

                    // atoms directly connected to the roots are leaves
                    // (unless they are already roots)
                    for (IAtom possibleLeaf : atomContainer.getConnectedAtomsList(rootAtom)) {
                        if (!component.contains(possibleLeaf)) {
                            rootSystem.addLeaf(possibleLeaf);
                        }
                    }
                }
                rootSystems.add(rootSystem);
                currentLabel++;
            }
        }
        return rootSystems;
    }

    /**
     * Depth-first search through an atomContainer to find root systems, which
     * are essentially a specialized type of connected component.
     *
     * @param atomV the current atom
     * @param atomU an atom connected to the current atom
     * @param cLabel the current component label
     * @param labels the component labels for each atom
     * @param atomContainer an atomContainer
     * @param roots the signature roots
     * @param component the component that is being built
     */
    private static void dfs(IAtom atomV, IAtom atomU,
            int cLabel, int[] labels,
            IAtomContainer atomContainer, List<IAtom> roots,
            List<IAtom> component) {
        int uIndex = atomContainer.indexOf(atomU);
        if (atomV == null || (roots.contains(atomU) && labels[uIndex] == 0)) {
            labels[uIndex] = cLabel;
            component.add(atomU);
            for (IAtom atomW : atomContainer.getConnectedAtomsList(atomU)) {
                if (atomW != atomV) {
                    dfs(
                            atomU, atomW, cLabel, labels, atomContainer, roots, component);
                }
            }
        } else {
        }
    }

    private SignatureRootFinder() {
    }
}
