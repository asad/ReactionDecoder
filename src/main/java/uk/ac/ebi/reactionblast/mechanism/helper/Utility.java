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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.Serializable;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.CASE_INSENSITIVE_ORDER;
import java.util.ArrayList;
import static java.util.Arrays.sort;
import java.util.Collection;
import static java.util.Collections.sort;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.silent.RingSet;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import org.openscience.smsd.helper.MoleculeInitializer;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct;
import uk.ac.ebi.reactionblast.signature.RBlastMoleculeSignature;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeDayLight;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class Utility extends MatrixPrinter implements Serializable {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Utility.class);

    /**
     * Used Chemaxon to generate smikrs
     *
     * @param reaction
     * @param remove_AAM
     * @return
     */
    public static String getSMILES(IReaction reaction, boolean remove_AAM) {
        StringBuilder sb = new StringBuilder("");
        try {
            for (IAtomContainer mol : reaction.getReactants().atomContainers()) {
                sb.append(getSMILES(mol, true));
            }

            sb.append(">>");

            for (IAtomContainer mol : reaction.getProducts().atomContainers()) {
                sb.append(getSMILES(mol, true));
            }

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return sb.toString();
    }

    /**
     * Used CDK to generate smiles
     *
     * @param mol
     * @param atom
     * @param level
     * @param remove_AAM
     * @return
     * @throws Exception
     */
    public static String getCircularSMILES(
            IAtomContainer mol, IAtom atom, int level, boolean remove_AAM) throws Exception {
        int refAtom = getAtomIndexByID(mol, atom);
        IAtomContainer fragment = getCircularFragment(mol, refAtom, level);
        String smiles = getSMILES(fragment, remove_AAM);
        return smiles;
    }

    /**
     * Used CDK to generate smiles
     *
     * @param mol
     * @param remove_AAM
     * @return
     */
    public static String getSMILES(
            IAtomContainer mol, boolean remove_AAM) {
        String smiles = "";
        try {
            return new uk.ac.ebi.reactionblast.tools.CDKSMILES(mol, true, remove_AAM).getCanonicalSMILES();
        } catch (CloneNotSupportedException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return smiles;
    }

    /**
     *
     * @param mol
     * @return
     */
    protected static List<IAtom> getAtoms(
            IAtomContainer mol) {
        List<IAtom> atoms = new ArrayList<>(mol.getAtomCount());
        for (IAtom atom : mol.atoms()) {
            atoms.add(atom);
        }
        return atoms;
    }

    /**
     *
     * @param bond
     * @param molset
     * @return
     */
    protected static String getMoleculeID(IBond bond, IAtomContainerSet molset) {
        for (IAtomContainer mol : molset.atomContainers()) {
            if (mol.contains(bond)) {
                return mol.getID();
            }
        }
        return null;
    }

    /**
     *
     * @param atom
     * @param molset
     * @return
     */
    protected static String getMoleculeID(IAtom atom, IAtomContainerSet molset) {
        for (IAtomContainer mol : molset.atomContainers()) {
            if (mol.contains(atom)) {
                return mol.getID();
            }
        }
        return null;
    }

    /**
     *
     * @param bond
     * @param molset
     * @return
     */
    protected static IAtomContainer getAtomContainer(IBond bond, IAtomContainerSet molset) {
        for (IAtomContainer mol : molset.atomContainers()) {
            if (mol.contains(bond)) {
                return mol;
            }
        }
        return null;
    }

    /**
     *
     * @param atom
     * @param molset
     * @return
     */
    protected static IAtomContainer getAtomContainer(IAtom atom, IAtomContainerSet molset) {
        for (IAtomContainer mol : molset.atomContainers()) {
            if (mol.contains(atom)) {
                return mol;
            }
        }
        return null;
    }

    /**
     * Return atom by ID match
     *
     * @param molWithoutH
     * @param refAtom
     * @return
     */
    protected static int getAtomIndexByID(IAtomContainer molWithoutH, IAtom refAtom) {
        for (IAtom atom : molWithoutH.atoms()) {
            if (atom.getID().equalsIgnoreCase(refAtom.getID())) {
                return molWithoutH.indexOf(atom);
            }
        }
        return -1;
    }

    /**
     * Return Signature of height h
     *
     * @param mol
     * @param atom
     * @param height
     * @return
     * @throws CloneNotSupportedException
     */
    protected static String getSignature(IAtomContainer mol, IAtom atom, int height) throws CloneNotSupportedException {
        IAtomContainer molWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(mol);
        int atomIndex = getAtomIndexByID(molWithoutH, atom);
        RBlastMoleculeSignature moleculeSignature = new RBlastMoleculeSignature(molWithoutH);
        moleculeSignature.setUseCharge(true);
        moleculeSignature.setBondSensitive(true);
        moleculeSignature.setUseAromatics(true);
        if (atomIndex >= 0) {
            return moleculeSignature.getAtomSignature(atomIndex, height).toCanonicalString();
        } else {
            return "";
        }
    }

    /**
     *
     * @param atomRCChangesMap
     * @param fragments
     * @throws CloneNotSupportedException
     */
    protected static void setFragmentMatches(SortedMap<String, Integer> atomRCChangesMap, List<IAtomContainer> fragments) throws CloneNotSupportedException {
        for (IAtomContainer fragment : fragments) {
            CountSubstructures countSubstructures = new CountSubstructures(fragment);
            atomRCChangesMap.keySet().stream().forEach((pattern) -> {
                int hit = 0;
                try {
                    hit = countSubstructures.substructureSize(pattern);
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                int val = hit == 0 ? 0 : atomRCChangesMap.get(pattern) + 1;
                atomRCChangesMap.put(pattern, val);
            });
        }
    }

    /**
     *
     * @param atomRCChangesMap
     * @param signatures
     */
    protected static void setReactionCenterMatches(IPatternFingerprinter atomRCChangesMap, List<String> signatures) {
        signatures.stream().forEach((fragment) -> {
            try {
                atomRCChangesMap.add(new Feature(fragment, 1.0));
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        });
    }

    /**
     *
     * @param rid
     * @param mol
     * @param atom
     * @param patternFP
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    protected static void setCircularSignatureFingerprints(String rid, IAtomContainer mol, IAtom atom, Map<Integer, IPatternFingerprinter> patternFP) throws CDKException, CloneNotSupportedException {

        for (int i = 1; i < 5; i++) {
            if (!patternFP.containsKey(i)) {
                IPatternFingerprinter fp = new PatternFingerprinter();
                fp.setFingerprintID(rid + ":" + "Signature: " + i);
                patternFP.put(i, fp);
            }
            String signature = getSignature(mol, atom, i);
            patternFP.get(i).add(new Feature(signature, 1.0));
        }

        if (!patternFP.containsKey(-1)) {
            IPatternFingerprinter fp = new PatternFingerprinter();
            fp.setFingerprintID(rid + ":" + "Signature: " + -1);
            patternFP.put(-1, fp);
        }
        String signature = getSignature(mol, atom, -1);
        patternFP.get(-1).add(new Feature(signature, 1.0));
    }

    /**
     *
     * @param reactBond
     * @param prodBond
     * @return
     */
    protected static String getCanonisedBondChangePattern(IBond reactBond, IBond prodBond) {

        String concatE = getCanonicalisedBondChangePattern(reactBond);
        String concatP = getCanonicalisedBondChangePattern(prodBond);

        List<String> pattern = new ArrayList<>(2);
        pattern.add(0, concatE);
        pattern.add(1, concatP);

        sort(pattern, CASE_INSENSITIVE_ORDER);
        return pattern.get(0).concat("*").concat(pattern.get(1));

    }

    /**
     *
     * @param bond
     * @return
     */
    protected static String getCanonicalisedBondChangePattern(IBond bond) {
        String symbol = getBondOrderSign(bond);
        List<String> atoms = new ArrayList<>(2);
        atoms.add(0, bond.getAtom(0).getSymbol());
        atoms.add(1, bond.getAtom(1).getSymbol());
        sort(atoms, CASE_INSENSITIVE_ORDER);
        String concatenatedSymbols = atoms.get(0).concat(symbol).concat(atoms.get(1));
        return concatenatedSymbols.trim();
    }

    /**
     *
     * @param bond
     * @return
     */
    public static String getBondOrderSign(IBond bond) {
        String bondSymbol = "";
        if (bond.getFlag(ISAROMATIC)) {
            bondSymbol += "@";
        } else if (bond.getFlag(ISINRING)) {
            bondSymbol += "%";
        } else if (bond.getOrder() == SINGLE) {
            bondSymbol += "-";
        } else if (bond.getOrder() == DOUBLE) {
            bondSymbol += "=";
        } else if (bond.getOrder() == TRIPLE) {
            bondSymbol += "#";
        } else if (bond.getOrder() == QUADRUPLE) {
            return "$";
        }

        return bondSymbol;
    }

    /**
     *
     * @param ringBond
     * @param singleRings
     * @return
     */
    public static int getNeighbourBondOrderCountFromRing(IBond ringBond, IRingSet singleRings) {
        int minValue = 9999;
        for (IAtomContainer ring : singleRings.atomContainers()) {
            int value = 0;
            if (ring.contains(ringBond.getAtom(0)) && ring.contains(ringBond.getAtom(1))) {
                for (IBond bond : ring.bonds()) {
                    if (bond.contains(ringBond.getAtom(0)) || bond.contains(ringBond.getAtom(1))) {
                        value += bond.getOrder().numeric();
                    }
                }
            }
            if (value < minValue) {
                minValue = value;
            }
        }
        return minValue;
    }

    /**
     *
     * @param ringBond
     * @param singleRings
     * @return
     */
    public static IRingSet getSmallestRingSet(IBond ringBond, IRingSet singleRings) {
        IRingSet rs = new RingSet();
        for (IAtomContainer ring : singleRings.atomContainers()) {
            if (ring.contains(ringBond.getAtom(0)) && ring.contains(ringBond.getAtom(1))) {
                if (rs.getAtomContainerCount() == 0) {
                    rs.addAtomContainer(ring);
                    continue;
                }
                for (IAtomContainer smallestRing : rs.atomContainers()) {
                    if (ring.getAtomCount() == smallestRing.getAtomCount()) {
                        if (!rs.contains(ring)) {
                            rs.addAtomContainer(ring);
                        }
                    } else if (ring.getAtomCount() < smallestRing.getAtomCount()) {
                        rs.removeAllAtomContainers();
                        rs.addAtomContainer(ring);
                    }
                }
            }
        }
        return rs;
    }

    /**
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    protected static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        MoleculeInitializer.initializeMolecule(atomContainer);
    }

    /**
     *
     * @param mol
     * @param startAtomIndex
     * @param radius
     * @return
     * @throws Exception
     */
    public static IAtomContainer getCircularFragment(IAtomContainer mol, int startAtomIndex, int radius) throws Exception {
        IAtomContainer fragment = cloneWithIDs(mol);
        Set<IAtom> removeList = new HashSet<>();
        Collection<IAtom> solutionSphereList = circularFragment(fragment, startAtomIndex, radius);

        for (IAtom atom : fragment.atoms()) {
            if (!solutionSphereList.contains(atom)) {
                removeList.add(atom);
            }
        }

        for (Iterator<IAtom> it = removeList.iterator(); it.hasNext();) {
            fragment.removeAtom(it.next());
        }

        IAtomContainer canonicalise = canonicalise(fragment);
        aromatizeDayLight(canonicalise);

        return fragment;
    }

    /**
     *
     * @param org_mol
     * @return cloned canonicalized molecule
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public static IAtomContainer canonicalise(IAtomContainer org_mol) throws CloneNotSupportedException, CDKException {

        IAtomContainer cloneMolecule = cloneWithIDs(org_mol);

        int[] p = new int[cloneMolecule.getAtomCount()];
//        /*
//         Signature based canonical Permutations
//         */
//        p = new uk.ac.ebi.reactionblast.tools.labelling.
//        SignatureMoleculeLabeller().getCanonicalPermutation(cloneMolecule);

        /*
        Use the Canonical labelling from the SMILES
        IMP: Suggested by John May
         */
        try {
            unique().create(cloneMolecule, p);
        } catch (CDKException e) {
            e.printStackTrace();
        }

        permuteWithoutClone(p, cloneMolecule);

        /*
        Set the IDs to container
         */
        if (org_mol.getID() != null) {
            cloneMolecule.setID(org_mol.getID());
        }

        return cloneMolecule;
    }

    /*
    This is a very imp code modified by John May
    The idea is to canonicalise the atoms and bonds
     */
    private static void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
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

    /**
     * Performs a breadthFirstSearch in an AtomContainer starting with a
     * particular sphere, which usually consists of one start atom. While
     * searching the graph, the method marks each visited atom. It then puts all
     * the atoms connected to the atoms in the given sphere into a new vector
     * which forms the sphere to search for the next recursive method call. All
     * atoms that have been visited are put into a molecule container. This
     * breadthFirstSearch does thus find the connected graph for a given start
     * atom.
     *
     * @param atomContainer The AtomContainer to be searched
     * @param rootAtom
     * @param max
     * @return
     * @throws CDKException
     */
    public static Collection<IAtom> circularFragment(IAtomContainer atomContainer, int rootAtom, int max) throws CDKException {
        IAtom root = atomContainer.getAtom(rootAtom);
        Set<IAtom> paths = new HashSet<>();
        // list of visited nodes
        LinkedList<IAtom> closedList = new LinkedList<>();
        // list of nodes to visit (sorted)
        LinkedList<IAtom> openList = new LinkedList<>();
        openList.add(root);

        // list of nodes to visit (sorted)
        LinkedList<IAtom> neighbours = new LinkedList<>();

        int level = 0;
        while (!openList.isEmpty()) {
            IAtom currentPath = openList.removeFirst();

            // path found!
            paths.add(currentPath);
            closedList.add(currentPath);

//            System.out.println("level " + (level + "/" + max) + " atom " + currentPath.getSymbol());
            // addBinary neighbors to the open list
            neighbours.addAll(atomContainer.getConnectedAtomsList(currentPath));

            if (openList.isEmpty() && !neighbours.isEmpty() && (max > level || max == -1)) {
                neighbours.stream().filter((a) -> (!closedList.contains(a))).forEach((a) -> {
                    openList.add(a);
                });
                level += 1;
                neighbours.clear();
            }
        }
        return paths;
    }

    /**
     *
     * @param rid
     * @param molOrignal
     * @param atom
     * @param patternFP
     * @throws Exception
     * @throws CloneNotSupportedException
     */
    protected static void setCircularFingerprints(String rid,
            IAtomContainer molOrignal,
            IAtom atom, Map<Integer, IPatternFingerprinter> patternFP)
            throws Exception, CloneNotSupportedException {
        IAtomContainer clone = molOrignal.clone();
        for (int i = 0; i < 3; i++) {
            if (!patternFP.containsKey(i)) {
                IPatternFingerprinter fp = new PatternFingerprinter();
                fp.setFingerprintID(rid + ":" + "Signature: " + i);
                patternFP.put(i, fp);
            }
            String circularSMILES = getCircularSMILES(clone, atom, i, true);
            patternFP.get(i).add(new Feature(circularSMILES, 1.0));
        }
        if (!patternFP.containsKey(-1)) {
            IPatternFingerprinter fp = new PatternFingerprinter();
            fp.setFingerprintID(rid + ":" + "Signature: " + -1);
            patternFP.put(-1, fp);
        }

        String circularSMILES = getCircularSMILES(clone, atom, -1, true);
        patternFP.get(-1).add(new Feature(circularSMILES, 1.0));

    }

    /**
     *
     * @param molOrignal
     * @param type
     * @param atom
     * @return
     * @throws Exception
     * @throws CloneNotSupportedException
     */
    protected static List<ReactionCenterFragment> getCircularReactionPatternFingerprints(IAtomContainer molOrignal,
            IAtom atom,
            EnumSubstrateProduct type)
            throws Exception, CloneNotSupportedException {
        List<ReactionCenterFragment> fragmentsRC = new ArrayList<>();
        IAtomContainer clone = molOrignal.clone();
        for (int i = 0; i < 3; i++) {
            String smiles = getCircularSMILES(clone, atom, i, true);
            ReactionCenterFragment reactionCenterFragment = new ReactionCenterFragment(smiles, i, type);
//            System.out.println(reactionCenterFragment + " smiles " + smiles);
            fragmentsRC.add(reactionCenterFragment);
        }
        String smiles = getCircularSMILES(clone, atom, -1, true);
        ReactionCenterFragment reactionCenterFragment = new ReactionCenterFragment(smiles, -1, type);
        fragmentsRC.add(reactionCenterFragment);
        return fragmentsRC;
    }
}
