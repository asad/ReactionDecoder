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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.Serializable;
import static java.lang.String.CASE_INSENSITIVE_ORDER;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

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
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.interfaces.EnumSubstrateProduct;
import uk.ac.ebi.reactionblast.signature.RBlastMoleculeSignature;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import static java.lang.System.err;
import static org.openscience.cdk.CDKConstants.REACTIVE_CENTER;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.smsd.helper.MoleculeInitializer;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS.BOND_STEREO;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.ATOM_STEREO_CHANGE_INFORMATION;
import static uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS.BOND_CHANGE_INFORMATION;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.E;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.EITHER;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.NONE;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.R;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.S;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.Z;
import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Arrays.sort;
import static java.util.Collections.sort;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.Reaction;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import uk.ac.ebi.reactionblast.mechanism.StereoChange;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeCDK;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class ReactionMappingUtility implements Serializable {

    /**
     *
     * @param reactantAtom
     * @param productAtom
     * @param atomContainerR
     * @param atomContainerP
     * @return
     * @throws Exception
     */
    protected static MoleculeMoleculePair getMolMolPair(
            IAtom reactantAtom,
            IAtom productAtom,
            IAtomContainer atomContainerR,
            IAtomContainer atomContainerP) throws Exception {

        int atomIndexR = getAtomIndexByID(atomContainerR, reactantAtom);

        String signatureR1 = getSignature(atomContainerR, reactantAtom, 1);
        String signatureR2 = getSignature(atomContainerR, reactantAtom, 2);
        String signatureR3 = getSignature(atomContainerR, reactantAtom, 3);
        String signatureR = getSignature(atomContainerR, reactantAtom, -1);

        IAtomContainer fragR1 = getCircularFragment(atomContainerR, atomIndexR, 1);
        IAtomContainer fragR2 = getCircularFragment(atomContainerR, atomIndexR, 2);
        IAtomContainer fragR3 = getCircularFragment(atomContainerR, atomIndexR, 3);
        IAtomContainer fragR = getCircularFragment(atomContainerR, atomIndexR, -1);

        String signatureP1 = getSignature(atomContainerP, productAtom, 1);
        String signatureP2 = getSignature(atomContainerP, productAtom, 2);
        String signatureP3 = getSignature(atomContainerP, productAtom, 3);
        String signatureP = getSignature(atomContainerP, productAtom, -1);

        int atomIndexP = getAtomIndexByID(atomContainerP, productAtom);

        IAtomContainer fragP1 = getCircularFragment(atomContainerP, atomIndexP, 1);
        IAtomContainer fragP2 = getCircularFragment(atomContainerP, atomIndexP, 2);
        IAtomContainer fragP3 = getCircularFragment(atomContainerP, atomIndexP, 3);
        IAtomContainer fragP = getCircularFragment(atomContainerP, atomIndexP, -1);

        IReaction reaction1 = new Reaction();
        reaction1.addReactant(fragR1, 1.0);
        reaction1.addProduct(fragP1, 1.0);

        IReaction reaction2 = new Reaction();
        reaction2.addReactant(fragR2, 1.0);
        reaction2.addProduct(fragP2, 1.0);

        IReaction reaction3 = new Reaction();
        reaction3.addReactant(fragR3, 1.0);
        reaction3.addProduct(fragP3, 1.0);

        IReaction reaction = new Reaction();
        reaction.addReactant(fragR, 1.0);
        reaction.addProduct(fragP, 1.0);

        String smirks;
        String smirks1;
        String smirks2;
        String smirks3;

        smirks = getSMILES(reaction, true);
        smirks1 = getSMILES(reaction1, true);
        smirks2 = getSMILES(reaction2, true);
        smirks3 = getSMILES(reaction3, true);

        String smartsR;
        String smartsR1;
        String smartsR2;
        String smartsR3;

        smartsR = getSMILES(fragR, true);
        smartsR1 = getSMILES(fragR1, true);
        smartsR2 = getSMILES(fragR2, true);
        smartsR3 = getSMILES(fragR3, true);

        String smartsP;
        String smartsP1;
        String smartsP2;
        String smartsP3;

        smartsP = getSMILES(fragP, true);
        smartsP1 = getSMILES(fragP1, true);
        smartsP2 = getSMILES(fragP2, true);
        smartsP3 = getSMILES(fragP3, true);

        ReactantProductPair rrpName = new ReactantProductPair(atomContainerR.getID(), atomContainerP.getID());
        ReactantProductPair rrpSMARTS = new ReactantProductPair(smartsR, smartsP);
        ReactantProductPair rrpSignature = new ReactantProductPair(signatureR, signatureP);

        MoleculeMoleculePair mmp = new MoleculeMoleculePair(rrpName, rrpSMARTS, rrpSignature, smirks);
        mmp.setSignature1(new ReactantProductPair(signatureR1, signatureP1));
        mmp.setSignature2(new ReactantProductPair(signatureR2, signatureP2));
        mmp.setSignature3(new ReactantProductPair(signatureR3, signatureP3));
        mmp.setSmarts1(new ReactantProductPair(smartsR1, smartsP1));
        mmp.setSmarts2(new ReactantProductPair(smartsR2, smartsP2));
        mmp.setSmarts3(new ReactantProductPair(smartsR3, smartsP3));
        mmp.setSmirks1(smirks1);
        mmp.setSmirks2(smirks2);
        mmp.setSmirks3(smirks3);

        return mmp;
    }

    /**
     * Used Chemaxon to generate SMIRKS
     *
     * @param reaction
     * @param remove_AAM
     * @return
     */
    protected static String getSMILES(IReaction reaction, boolean remove_AAM) {
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
            getLogger(ReactionMappingUtility.class.getName()).log(SEVERE, null, ex);
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
    protected static String getCircularSMILES(
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
    protected static String getSMILES(
            IAtomContainer mol, boolean remove_AAM) {
        String smiles = "";
        try {
            return new uk.ac.ebi.reactionblast.tools.CDKSMILES(mol, true, remove_AAM).getCanonicalSMILES();
        } catch (CloneNotSupportedException ex) {
            getLogger(ReactionMappingUtility.class.getName()).log(SEVERE, null, ex);
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
                    getLogger(ReactionMappingUtility.class.getName()).log(SEVERE, null, ex);
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
                getLogger(ReactionMappingUtility.class.getName()).log(SEVERE, null, ex);
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
    protected static String getCanonicalisedBondChangePattern(IBond reactBond, IBond prodBond) {

        String concatE = ReactionMappingUtility.getCanonicalisedBondChangePattern(reactBond);
        String concatP = ReactionMappingUtility.getCanonicalisedBondChangePattern(prodBond);

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
    protected static String getBondOrderSign(IBond bond) {
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
    protected static int getNeighbourBondOrderCountFromRing(IBond ringBond, IRingSet singleRings) {
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
    protected static IRingSet getSmallestRingSet(IBond ringBond, IRingSet singleRings) {
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
    protected static IAtomContainer getCircularFragment(IAtomContainer mol, int startAtomIndex, int radius) throws Exception {
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
        aromatizeCDK(canonicalise);

        return fragment;
    }

    /**
     *
     * @param org_mol
     * @return cloned canonicalized molecule
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    protected static IAtomContainer canonicalise(IAtomContainer org_mol) throws CloneNotSupportedException, CDKException {

        IAtomContainer cloneMolecule = cloneWithIDs(org_mol);

        int[] p = new int[cloneMolecule.getAtomCount()];
//        /*
//         Signature based canonical Permutations
//         */
//        p = new uk.ac.ebi.reactionblast.tools.labelling.
//        SignatureMoleculeLabeller().getCanonicalPermutation(cloneMolecule);

        /*
         * Use the Canonical labelling from the SMILES
         * IMP: Suggested by John May
         */
        unique().create(cloneMolecule, p);

        permuteWithoutClone(p, cloneMolecule);

        /*
         * Set the IDs to container
         */
        if (org_mol.getID() != null) {
            cloneMolecule.setID(org_mol.getID());
        }

        return cloneMolecule;
    }

    /*
     * This is a very imp code modified by John May
     * The idea is to canonicalise the atoms and bonds
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
            err.println("pokemon!");
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
    protected static Collection<IAtom> circularFragment(IAtomContainer atomContainer, int rootAtom, int max) throws CDKException {
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

    protected static Set<IBond> getBondCleavedFormedChanges(IReaction mappedReaction, Map<IAtom, IAtom> mappings) {

        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        Set<IBond> reactantsbonds = getBonds(allReactants);
        Set<IBond> productsbonds = getBonds(allProducts);

        Set<IBond> bondChange = detectBondsCleavedAndFormed(reactantsbonds, productsbonds, mappings);

        return bondChange;
    }

    protected static Map<IBond, IBond> getBondOrderChanges(IReaction mappedReaction, Map<IAtom, IAtom> mappings) {

        IAtomContainerSet allReactants = ExtReactionManipulatorTool.getAllReactants(mappedReaction);
        IAtomContainerSet allProducts = ExtReactionManipulatorTool.getAllProducts(mappedReaction);
        Set<IBond> reactantsbonds = getBonds(allReactants);
        Set<IBond> productsbonds = getBonds(allProducts);

        IRingSet ringsR = getRings(allReactants);
        IRingSet ringsP = getRings(allProducts);

        Map<IBond, IBond> detectBondOrderChanges = detectBondOrderChanges(reactantsbonds, productsbonds, mappings, ringsR, ringsP);

        return detectBondOrderChanges;
    }

    protected static Set<IAtom> getAtomStereoChanges(IReaction mappedReaction, Map<IAtom, IAtom> mappings) {
        Set<IAtom> detectStereoChanges = detectStereoChanges(mappedReaction, mappings);
        return detectStereoChanges;
    }

    private static IRingSet getRings(IAtomContainerSet containerSet) {
        IRingSet rings = new org.openscience.cdk.RingSet();
        for (IAtomContainer container : containerSet.atomContainers()) {
            try {
                /*
                 * set Flag(CDKConstants.ISINRING)
                 */
                initializeMolecule(container);
            } catch (CDKException ex) {
                getLogger(ReactionMappingUtility.class.getName()).log(SEVERE, null, ex);
            }
            IRingSet ssrRings = new SSSRFinder(container).findSSSR();
            rings.add(ssrRings);
        }
        return rings;
    }

    protected static Map<IAtom, IAtom> getMappingsWithouH(IReaction mappedReaction) {
        HashMap<IAtom, IAtom> mappingsWithoutH = new HashMap<>();
        Map<IAtom, IAtom> mappingsFull = getMappings(mappedReaction);
        mappingsFull.entrySet().stream().filter((map) -> (!map.getKey().getSymbol().equals("H"))).forEach((map) -> {
            mappingsWithoutH.put(map.getKey(), map.getValue());
        });
        return mappingsWithoutH;
    }

    protected static Map<IAtom, IAtom> getMappings(IReaction mappedReaction) {
        HashMap<IAtom, IAtom> mappings = new HashMap<>();

        /*
         * Find mapped atoms in reactants
         */
        for (IAtomContainer container : ExtReactionManipulatorTool.getAllReactants(mappedReaction).atomContainers()) {
            for (IAtom a : container.atoms()) {
                Integer mappingNumber = a.getProperty(ATOM_ATOM_MAPPING);
                if (mappingNumber != null) {
                    mappings.put(a, null);
                }
            }
        }
        /*
         * Find mapped atoms in products
         */
        for (IAtomContainer container : ExtReactionManipulatorTool.getAllProducts(mappedReaction).atomContainers()) {
            for (IAtom a : container.atoms()) {
                Integer mappingNumber = a.getProperty(ATOM_ATOM_MAPPING);
                if (mappingNumber != null) {
                    for (IAtom mappedAtom : mappings.keySet()) {
                        Integer storedAtomMappingNumber = mappedAtom.getProperty(ATOM_ATOM_MAPPING);
                        if (mappingNumber.intValue() == storedAtomMappingNumber) {
                            mappings.put(mappedAtom, a);
                            break;
                        }
                    }
                }
            }
        }

        /*
         * Find mapping pairs
         */
        Set<IAtom> unMappedAtomsToBeRemoved = new HashSet<>();
        mappings.entrySet().stream().filter((map) -> (map.getValue() == null)).forEach((map) -> {
            unMappedAtomsToBeRemoved.add(map.getKey());
        });

        /*
         * Removed unpaired atoms
         */
        unMappedAtomsToBeRemoved.stream().forEach((removeAtom) -> {
            mappings.remove(removeAtom);
        });

        return mappings;
    }

    protected static Set<IBond> getBonds(IAtomContainerSet containers) {
        Set<IBond> bonds = new LinkedHashSet<>();
        for (IAtomContainer container : containers.atomContainers()) {
            IBond[] bondArray = ExtAtomContainerManipulator.getBondArray(container);
            bonds.addAll(Arrays.asList(bondArray));
        }
        return bonds;
    }

    /**
     * Reports bond formed or cleaved
     *
     * @param reactantsbonds
     * @param productsbonds
     * @param mappings
     * @return
     */
    private static Set<IBond> detectBondsCleavedAndFormed(Set<IBond> reactantsbonds, Set<IBond> productsbonds, Map<IAtom, IAtom> mappings) {
        Set<IBond> bondChange = new LinkedHashSet<>();

        Set<IBond> cleaved = new LinkedHashSet<>(reactantsbonds);

        reactantsbonds.stream().filter((rb) -> (mappings.containsKey(rb.getAtom(0)) && mappings.containsKey(rb.getAtom(1)))).forEach((rb) -> {
            productsbonds.stream().filter((pb) -> (mappings.containsValue(pb.getAtom(0)) && mappings.containsValue(pb.getAtom(1)))).map((pb) -> {
                if (rb.getAtom(0).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(0).getProperty(ATOM_ATOM_MAPPING))) {
                    if (rb.getAtom(1).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(1).getProperty(ATOM_ATOM_MAPPING))) {
                        cleaved.remove(rb);
                    }
                }
                return pb;
            }).filter((pb) -> (rb.getAtom(1).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(0).getProperty(ATOM_ATOM_MAPPING)))).filter((pb) -> (rb.getAtom(0).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(1).getProperty(ATOM_ATOM_MAPPING)))).forEach((_item) -> {
                cleaved.remove(rb);
            });
        });

        reactantsbonds.stream().filter((bond) -> (!mappings.containsKey(bond.getAtom(0)) && !mappings.containsKey(bond.getAtom(1)))).forEach((bond) -> {
            cleaved.remove(bond);
        });

        Set<IBond> formed = new LinkedHashSet<>(productsbonds);

        productsbonds.stream().filter((pb) -> (mappings.containsValue(pb.getAtom(0)) && mappings.containsValue(pb.getAtom(1)))).forEach((pb) -> {
            reactantsbonds.stream().filter((rb) -> (mappings.containsKey(rb.getAtom(0)) && mappings.containsKey(rb.getAtom(1)))).map((rb) -> {
                if (rb.getAtom(0).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(0).getProperty(ATOM_ATOM_MAPPING))) {
                    if (rb.getAtom(1).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(1).getProperty(ATOM_ATOM_MAPPING))) {
                        formed.remove(pb);
                    }
                }
                return rb;
            }).filter((rb) -> (rb.getAtom(1).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(0).getProperty(ATOM_ATOM_MAPPING)))).filter((rb) -> (rb.getAtom(0).getProperty(ATOM_ATOM_MAPPING).equals(pb.getAtom(1).getProperty(ATOM_ATOM_MAPPING)))).forEach((_item) -> {
                formed.remove(pb);
            });
        });

        productsbonds.stream().filter((bond) -> (!mappings.containsValue(bond.getAtom(0)) && !mappings.containsValue(bond.getAtom(1)))).forEach((bond) -> {
            formed.remove(bond);
        });

        /*
         * set Flags for Bonds cleaved and formed
         */
        cleaved.stream().map((bond) -> {
            bond.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
            return bond;
        }).map((bond) -> {
            bond.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
            return bond;
        }).map((bond) -> {
            bond.getAtom(0).setFlag(REACTIVE_CENTER, true);
            return bond;
        }).map((bond) -> {
            bond.getAtom(1).setFlag(REACTIVE_CENTER, true);
            return bond;
        }).map((bond) -> {
            bond.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED);
            return bond;
        }).forEach((bond) -> {
            bondChange.add(bond);
        });

        /*
         * set Flags for Bonds formed
         */
        formed.stream().map((bond) -> {
            bond.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
            return bond;
        }).map((bond) -> {
            bond.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
            return bond;
        }).map((bond) -> {
            bond.getAtom(0).setFlag(REACTIVE_CENTER, true);
            return bond;
        }).map((bond) -> {
            bond.getAtom(1).setFlag(REACTIVE_CENTER, true);
            return bond;
        }).map((bond) -> {
            bond.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED);
            return bond;
        }).forEach((bond) -> {
            bondChange.add(bond);
        });

        cleaved.clear();
        formed.clear();
        return bondChange;
    }

    private static IAtom getKeyFromValue(IAtom a, Map<IAtom, IAtom> mappings) {
        for (Map.Entry<IAtom, IAtom> map : mappings.entrySet()) {
            if (map.getValue() == a) {
                return map.getKey();
            }
        }
        return null;
    }

    /**
     * Reports bond order changes BOND_ORDER_REDUCED or BOND_ORDER_GAIN Look for
     * BOND_CHANGE_INFORMATION Flag for BOND_ORDER_REDUCED or BOND_ORDER_GAIN
     *
     * @param reactantsbonds
     * @param productsbonds
     * @param mappings
     * @param queryRingSet
     * @param targetRingSet
     * @return
     */
    private static Map<IBond, IBond> detectBondOrderChanges(Set<IBond> reactantsbonds, Set<IBond> productsbonds, Map<IAtom, IAtom> mappings, IRingSet queryRingSet, IRingSet targetRingSet) {
        Map<IBond, IBond> bondChange = new HashMap<>();

        reactantsbonds.stream().filter((rb) -> (mappings.containsKey(rb.getAtom(0)) && mappings.containsKey(rb.getAtom(1)))).forEach((IBond rb) -> {
            productsbonds.stream().filter((pb) -> (pb.contains(mappings.get(rb.getAtom(0))) && pb.contains(mappings.get(rb.getAtom(1))))).map((IBond pb) -> {
                if ((isBondMappingMatch(rb, pb) && !rb.getOrder().equals(pb.getOrder()))) {
                    if (rb.getOrder().numeric() > pb.getOrder().numeric()) {
                        rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_REDUCED);
                        rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_REDUCED);
                        rb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                    } else if (rb.getOrder().numeric() < pb.getOrder().numeric()) {
                        rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                        rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                        rb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                    }

                    rb.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    rb.getAtom(1).setFlag(REACTIVE_CENTER, true);

                    if (pb.getOrder().numeric() > rb.getOrder().numeric()) {
                        pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_REDUCED);
                        pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_REDUCED);
                        pb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                    } else if (pb.getOrder().numeric() < rb.getOrder().numeric()) {
                        pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                        pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                        pb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);
                    }

                    pb.getAtom(0).setFlag(REACTIVE_CENTER, true);
                    pb.getAtom(1).setFlag(REACTIVE_CENTER, true);

                    /*
                    * Store all order changes
                     */
                    bondChange.put(rb, pb);
                }
                return pb;
            }).forEach((pb) -> {
                int kekuleEffect = isAlternateKekuleChange(rb, pb, queryRingSet, targetRingSet);
                if (kekuleEffect == 0) {
                    rb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                    rb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                    rb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);

                    pb.getAtom(0).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                    pb.getAtom(1).setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER_GAIN);
                    pb.setProperty(BOND_CHANGE_INFORMATION, ECBLAST_BOND_CHANGE_FLAGS.BOND_ORDER);

                    /*
                    * Store all order changes
                     */
                    bondChange.put(rb, pb);
                }
            });
        });
        return bondChange;
    }

    /**
     * Reports Stereo changes (R/S) or (E/Z) Look for BOND_CHANGE_INFORMATION
     * flags for BOND_STEREO Changes
     *
     * @param reaction
     * @return
     */
    private static Set<IAtom> detectStereoChanges(IReaction reaction, Map<IAtom, IAtom> mappings) {
        /*
         * Mining Stereo Atom Changes E/Z or R/S only
         */
        Set<IAtom> atomChanges = new LinkedHashSet<>();

        /*
         * Stereo mapping
         */
        Map<IAtom, IStereoAndConformation> chiralityCDK2D = new HashMap<>();
        try {
            chiralityCDK2D = getChirality2D(reaction, mappings);
        } catch (CDKException | CloneNotSupportedException ex) {
            err.println("WARNING: 2D CDK based stereo perception failed");
        }
        /*
         * Generate stereo information
         */
        List<StereoChange> stereogenicCenters = getStereoChanges(mappings, chiralityCDK2D);

        for (StereoChange sc : stereogenicCenters) {
            IAtom atomE = sc.getReactantAtom();
            IAtom atomP = sc.getProductAtom();

            IStereoAndConformation rsb = sc.getReactantAtomStereo();
            IStereoAndConformation psb = sc.getProductAtomStereo();

            if (atomE != null && atomP != null) {
                if (atomE.getSymbol().equals("P") || atomP.getSymbol().equals("P")) {
                    err.println("\nWARNING: The stereo change " + atomE.getSymbol()
                            + " not supported");
                    continue;
                }
                atomE.setFlag(REACTIVE_CENTER, true);
                atomP.setFlag(REACTIVE_CENTER, true);

                if ((sc.getReactantAtomStereo().equals(E)
                        || sc.getProductAtomStereo().equals(Z))
                        || (sc.getReactantAtomStereo().equals(Z)
                        || sc.getProductAtomStereo().equals(E))) {
                    atomE.setProperty(ATOM_STEREO_CHANGE_INFORMATION, rsb);
                    atomP.setProperty(ATOM_STEREO_CHANGE_INFORMATION, psb);

                    atomE.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                    atomP.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);

                } else if ((sc.getReactantAtomStereo().equals(R)
                        || sc.getProductAtomStereo().equals(S))
                        || (sc.getReactantAtomStereo().equals(S)
                        || sc.getProductAtomStereo().equals(R))) {
                    atomE.setProperty(ATOM_STEREO_CHANGE_INFORMATION, rsb);
                    atomP.setProperty(ATOM_STEREO_CHANGE_INFORMATION, psb);

                    atomE.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                    atomP.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO);
                }

                atomChanges.add(atomE);
                atomChanges.add(atomP);
            }
        }
        return atomChanges;

    }

    protected static boolean isBondMappingMatch(IBond a, IBond b) {
        if (isAtomMappingMatch(a.getAtom(1), b.getAtom(0))
                && isAtomMappingMatch(a.getAtom(0), b.getAtom(1))) {
            return true;
        }

        return (isAtomMappingMatch(a.getAtom(0), b.getAtom(0))
                && isAtomMappingMatch(a.getAtom(1), b.getAtom(1)));
    }

    protected static boolean isAtomMappingMatch(IAtom a, IAtom b) {
        return a.getProperty(ATOM_ATOM_MAPPING).equals(b.getProperty(ATOM_ATOM_MAPPING));
    }

    /**
     * Return is the ring has undergone a KekuleChange (0=KekuleChange, -1 or
     * 1=Bond Order Change)
     *
     * @param affectedBondReactants
     * @param affectedBondProducts
     * @param queryRingSet
     * @param targetRingSet
     * @return
     */
    private static int isAlternateKekuleChange(IBond affectedBondReactants, IBond affectedBondProducts, IRingSet queryRingSet, IRingSet targetRingSet) {
        if (affectedBondReactants != null && affectedBondProducts != null) {
            if (affectedBondReactants.getFlag(ISINRING)
                    == affectedBondProducts.getFlag(ISINRING)) {

                if ((!affectedBondReactants.getFlag(ISAROMATIC)
                        && affectedBondProducts.getFlag(ISAROMATIC))
                        || (affectedBondReactants.getFlag(ISAROMATIC)
                        && !affectedBondProducts.getFlag(ISAROMATIC))) {
                    IRingSet smallestRingSetR = getSmallestRingSet(affectedBondReactants, queryRingSet);
                    IRingSet smallestRingSetP = getSmallestRingSet(affectedBondProducts, targetRingSet);
                    int countR = getNeighbourBondOrderCountFromRing(affectedBondReactants, smallestRingSetR);
                    int countP = getNeighbourBondOrderCountFromRing(affectedBondProducts, smallestRingSetP);
                    int sizeR = smallestRingSetR.getAtomContainerCount() > 0 ? smallestRingSetR.atomContainers().iterator().next().getBondCount() : 0;
                    int sizeP = smallestRingSetP.getAtomContainerCount() > 0 ? smallestRingSetP.atomContainers().iterator().next().getBondCount() : 0;

                    IAtom atomR1 = affectedBondReactants.getAtom(0);
                    IAtom atomR2 = affectedBondReactants.getAtom(1);

                    IAtom atomP1 = affectedBondProducts.getAtom(0);
                    IAtom atomP2 = affectedBondProducts.getAtom(1);

                    if (atomR1.getID().equals(atomP1.getID())
                            && atomR2.getID().equals(atomP2.getID())) {
                        if (atomR1.getHybridization() != null
                                && atomR2.getHybridization() != null
                                && atomP1.getHybridization() != null
                                && atomP2.getHybridization() != null) {

                            if (!atomR1.getHybridization().equals(atomP1.getHybridization())
                                    && !atomR2.getHybridization().equals(atomP2.getHybridization())) {
                                if (countR == countP && sizeR == sizeP) {
                                    return 1;
                                } else {
                                    if (!affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                        if (abs(countR - countP) <= 1 && sizeR == sizeP) {
                                            return 1;
                                        }
                                    }
                                    if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                        if (abs(countR - countP) == 0 && sizeR == sizeP) {
                                            return 0;
                                        }
                                    }
                                }
                            } else if (atomR1.getHybridization().equals(atomP1.getHybridization())
                                    && atomR2.getHybridization().equals(atomP2.getHybridization())) {

                                if (!affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (abs(countR - countP) <= 1
                                            && sizeR == sizeP) {
                                        return 0;
                                    }
                                }
                                if (affectedBondReactants.getOrder().equals(affectedBondProducts.getOrder())) {
                                    if (abs(countR - countP) > 1 && sizeR == sizeP) {
                                        return 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return -1;
    }

    /**
     * This Chirality is based on the 2D with stereo code written by John May in
     * our collaboration. Note: Explicit Hydrogens should be added before
     * calling.
     *
     * @param reaction Mapped Reaction
     * @param mappings Atom-Atom Mapping
     * @return
     * @throws CDKException
     * @throws java.lang.CloneNotSupportedException
     */
    private static Map<IAtom, IStereoAndConformation> getChirality2D(IReaction reaction, Map<IAtom, IAtom> mappings) throws CDKException, CloneNotSupportedException {
        Map<IAtom, IStereoAndConformation> chiralityMap = new HashMap<>();

        for (IAtomContainer ac : ExtReactionManipulatorTool.getAllReactants(reaction).atomContainers()) {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
            for (IStereoElement stereoElement : ac.stereoElements()) {
                if (stereoElement instanceof ITetrahedralChirality) {
                    ITetrahedralChirality tetrahedralChirality = ((ITetrahedralChirality) stereoElement);
                    IAtom chiralAtom = tetrahedralChirality.getChiralAtom();
                    CIPTool.CIP_CHIRALITY cipChirality = CIPTool.getCIPChirality(ac, tetrahedralChirality);
                    IStereoAndConformation deduceChirality = deduceChirality(cipChirality);

                    if (mappings.containsKey(chiralAtom)) {
                        chiralityMap.put(chiralAtom, deduceChirality);
                    }
                }
                if (stereoElement instanceof IDoubleBondStereochemistry) {
                    IDoubleBondStereochemistry bondStereochemistry = (IDoubleBondStereochemistry) stereoElement;
                    IBond chiralBond = bondStereochemistry.getStereoBond();
                    CIPTool.CIP_CHIRALITY cipChirality = CIPTool.getCIPChirality(ac, bondStereochemistry);
                    IStereoAndConformation deduceChirality = deduceChirality(cipChirality);
                    if (mappings.containsKey(chiralBond.getAtom(0))) {
                        chiralityMap.put(chiralBond.getAtom(0), deduceChirality);
                    }
                    if (mappings.containsKey(chiralBond.getAtom(1))) {
                        chiralityMap.put(chiralBond.getAtom(1), deduceChirality);
                    }
                }
            }

        }

        for (IAtomContainer ac : ExtReactionManipulatorTool.getAllProducts(reaction).atomContainers()) {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
            for (IStereoElement stereoElement : ac.stereoElements()) {
                if (stereoElement instanceof ITetrahedralChirality) {
                    ITetrahedralChirality tetrahedralChirality = ((ITetrahedralChirality) stereoElement);
                    IAtom chiralAtom = tetrahedralChirality.getChiralAtom();
                    CIPTool.CIP_CHIRALITY cipChirality = CIPTool.getCIPChirality(ac, tetrahedralChirality);
                    IStereoAndConformation deduceChirality = deduceChirality(cipChirality);

                    if (mappings.containsValue(chiralAtom)) {
                        chiralityMap.put(chiralAtom, deduceChirality);
                    }
                }
                if (stereoElement instanceof IDoubleBondStereochemistry) {
                    IDoubleBondStereochemistry bondStereochemistry = (IDoubleBondStereochemistry) stereoElement;
                    IBond chiralBond = bondStereochemistry.getStereoBond();
                    CIPTool.CIP_CHIRALITY cipChirality = CIPTool.getCIPChirality(ac, bondStereochemistry);
                    IStereoAndConformation deduceChirality = deduceChirality(cipChirality);
                    if (mappings.containsValue(chiralBond.getAtom(0))) {
                        chiralityMap.put(chiralBond.getAtom(0), deduceChirality);
                    }
                    if (mappings.containsValue(chiralBond.getAtom(1))) {
                        chiralityMap.put(chiralBond.getAtom(1), deduceChirality);
                    }
                }
            }

        }
        chiralityMap.keySet().stream().forEach((atom) -> {
            atom.setProperty("Stereo", chiralityMap.get(atom));
        });

        return chiralityMap;
    }

    private static IStereoAndConformation deduceChirality(CIPTool.CIP_CHIRALITY cip) {

        IStereoAndConformation stereoType = IStereoAndConformation.NONE;

        switch (cip) {
            case R:
                stereoType = IStereoAndConformation.R;
                break;
            case S:
                stereoType = IStereoAndConformation.S;
                break;
            case E:
                stereoType = IStereoAndConformation.E;
                break;
            case Z:
                stereoType = IStereoAndConformation.Z;
                break;
            case NONE:
                stereoType = IStereoAndConformation.NONE;
                break;
            default:
                break;
        }

        return stereoType;
    }

    protected static IAtom getMappingAtomByID(IAtom atom, IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            if (atom.getProperty(ATOM_ATOM_MAPPING).equals(a.getProperty(ATOM_ATOM_MAPPING))) {
                return a;
            }
        }
        return null;
    }

    /**
     * Generates String from Stereo change types
     *
     * @param atom
     * @return (E/Z) or (R/S) else (NA)
     */
    protected static String getCanonicalisedAtomChangePattern(IAtom atom) {

        if (atom.getProperty("Stereo").equals(IStereoAndConformation.E)) {
            return atom.getSymbol().concat("(E/Z)");
        }

        if (atom.getProperty("Stereo").equals(IStereoAndConformation.Z)) {
            return atom.getSymbol().concat("(E/Z)");
        }

        if (atom.getProperty("Stereo").equals(IStereoAndConformation.R)) {
            return atom.getSymbol().concat("(R/S)");
        }

        if (atom.getProperty("Stereo").equals(IStereoAndConformation.S)) {
            return atom.getSymbol().concat("(R/S)");
        }

        return atom.getSymbol().concat("(NA)");
    }

    /**
     *
     * @param mappings
     * @param chirality2DCDK
     * @return
     */
    private static List<StereoChange> getStereoChanges(Map<IAtom, IAtom> mappings, Map<IAtom, IStereoAndConformation> chirality2DCDK) {

        List<StereoChange> stereoChangeList = new ArrayList<>();

        if (chirality2DCDK.isEmpty()) {
            return stereoChangeList;
        }

        if (mappings.isEmpty()) {
            return stereoChangeList;
        }

        for (Map.Entry<IAtom, IAtom> map : mappings.entrySet()) {
            if (!map.getKey().getSymbol().equalsIgnoreCase("H")
                    && chirality2DCDK.containsKey(map.getKey())
                    && chirality2DCDK.containsKey(map.getValue())) {
                if (isStereogenicChange(chirality2DCDK.get(map.getKey()), chirality2DCDK.get(map.getValue()))) {
                    StereoChange sc = new StereoChange(chirality2DCDK.get(map.getKey()), chirality2DCDK.get(map.getValue()), map.getKey(), map.getValue());
                    stereoChangeList.add(sc);
                }
            }
        }
        return stereoChangeList;
    }

    /**
     * Returns type of stereo changes
     *
     * @param a
     * @param b
     * @return
     */
    private static boolean isStereogenicChange(IStereoAndConformation a, IStereoAndConformation b) {
        if (a.equals(S) && b.equals(NONE)) {
            return true;
        } else if (a.equals(R) && b.equals(NONE)) {
            return true;
        } else if (b.equals(S) && a.equals(NONE)) {
            return true;
        } else if (b.equals(R) && a.equals(NONE)) {
            return true;
        } else if (a.equals(R) && b.equals(S)) {
            return true;
        } else if (a.equals(S) && b.equals(R)) {
            return true;
        } else if (a.equals(S) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(R) && b.equals(EITHER)) {
            return true;
        } else if (b.equals(S) && a.equals(EITHER)) {
            return true;
        } else if (b.equals(R) && a.equals(EITHER)) {
            return true;
        } else if (a.equals(EITHER) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(NONE) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(EITHER) && b.equals(NONE)) {
            return true;
        } else if (a.equals(Z) && b.equals(E)) {
            return true;
        } else if (a.equals(E) && b.equals(Z)) {
            return true;
        }
        return false;
    }

}
