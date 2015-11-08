/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.container;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReaction.Direction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.containers.MolContainer;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.interfaces.IStandardizer;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.AtomContainerSetComparator;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CDKReactionBuilder extends BasicDebugger implements Serializable, IStandardizer {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 19869866609698L;
    private static final Logger LOG = Logger.getLogger(CDKReactionBuilder.class.getName());
    private final IReactionSet reactionSet;
    private int moleculeCounter = 0; //Counter to create Unique Molecules
    private final Map<String, Double> stoichiometryMap;
    private final Map<String, BitSet> fingerprintMap;
    private final Map<String, IAtomContainer> moleculeMap;

    /**
     *
     * @throws java.lang.Exception
     */
    public CDKReactionBuilder() throws Exception {
        reactionSet = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
        stoichiometryMap = Collections.synchronizedMap(new HashMap<String, Double>());
        fingerprintMap = Collections.synchronizedMap(new HashMap<String, BitSet>());
        moleculeMap = Collections.synchronizedMap(new HashMap<String, IAtomContainer>());
    }

    @Override
    public String toString() {
        return "CDKReactionBuilder{" + "reactionSet=" + reactionSet + ", moleculeCounter="
                + moleculeCounter + ", stoichiometryMap=" + stoichiometryMap + ", fingerprintMap="
                + fingerprintMap + ", moleculeMap=" + moleculeMap + '}';
    }

    /**
     *
     * @param reactionSet
     * @throws java.lang.Exception
     */
    public synchronized void standardize(IReactionSet reactionSet) throws Exception {
        for (IReaction reaction : reactionSet.reactions()) {
            IReaction standardizedReaction = standardize(reaction);
            reactionSet.addReaction(standardizedReaction);
        }
    }

    /**
     *
     * @param reaction
     * @return
     * @throws Exception
     */
    @Override
    public synchronized IReaction standardize(IReaction reaction) throws Exception {

        if (DEBUG) {
            String createReactionSMILES = SmilesGenerator.unique().aromatic().createReactionSMILES(reaction);
            System.out.println("createReactionSMILES " + createReactionSMILES);
            System.out.println("standardize reaction module start");
        }

        List<IAtomContainer> _metabolites = new ArrayList<>();
        IReaction standardizedReaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);

        String reactionID = reaction.getID();
        int reactionCounter = 1;
        if (reactionID == null) {
            reactionID = "R" + Long.toString(reactionCounter++);
            reaction.setID(reactionID);
        }

//        System.out.println("****************************");
//        System.out.println("Processing: " + reactionID);
//        System.out.println("****************************");
//        System.out.println("Number of Reactant " + _imoledu.getAtomContainerCount());
//        System.out.println("Number of Product " + _imolpro.getAtomContainerCount());
        _metabolites.clear();

        standardizedReaction.setID(reactionID);

        stoichiometryMap.clear();

        Double tempStoic;

        if (DEBUG) {
            System.out.println("standardize reaction module phase 1");
        }
        for (IAtomContainer mol : reaction.getReactants().atomContainers()) {
            String id = mol.getID() == null || mol.getID().isEmpty() ? null : mol.getID();
            tempStoic = 1.0;
            if (reaction.getReactantCoefficient(mol) > 0) {
                tempStoic = reaction.getReactantCoefficient(mol);
            }

            if (DEBUG) {
                System.out.println("q_mol " + SmilesGenerator.unique().create(mol));
                System.out.println("standardize reaction module phase 1.1");
            }

            IAtomContainer gMol = ExtAtomContainerManipulator.cloneWithIDs(mol);
            if (DEBUG) {
                System.out.println("standardize reaction module phase 1.1.1");
            }
            ExtAtomContainerManipulator.fixDativeBonds(gMol);
            if (DEBUG) {
                System.out.println("standardize reaction module phase 1.1.2");
            }
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(gMol);
            IAtomContainer molWithH = gMol;
            //= ExtAtomContainerManipulator.addExplicitH(gMol);
            ExtAtomContainerManipulator.aromatizeMolecule(molWithH);

            if (DEBUG) {
                System.out.println("standardize reaction module phase 1.2");
            }

            if (id == null) {
                molWithH = setProperty(molWithH);
            } else {
                molWithH.setID(id);
            }

            if (DEBUG) {
                System.out.println("standardize reaction module phase 1.3");
                System.out.println("After Cleanup " + molWithH.getID());
                System.out.println(SmilesGenerator.unique().create(mol));
            }
            if (stoichiometryMap.containsKey(molWithH.getID())) {
                tempStoic += stoichiometryMap.get(molWithH.getID());
                stoichiometryMap.put(molWithH.getID(), tempStoic);
//                System.err.println("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            Collections.sort(_metabolites, comparator);
        } catch (Exception ex) {
            Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }

        setReactantMolecule(standardizedReaction, _metabolites);
        _metabolites.clear();

        if (DEBUG) {
            System.out.println("standardize reaction module phase 2");
        }
        if (DEBUG) {
            System.out.println();
            System.out.println("****************************");
            System.out.println();
        }
        for (IAtomContainer mol : reaction.getProducts().atomContainers()) {
            String id = mol.getID() == null || mol.getID().isEmpty() ? null : mol.getID();
            tempStoic = 1.0;
            if (reaction.getProductCoefficient(mol) > 0) {
                tempStoic = reaction.getProductCoefficient(mol);
            }
            if (DEBUG) {
                System.out.println("standardize reaction module phase 2.1");
                System.out.println("t_mol " + SmilesGenerator.unique().create(mol));
            }
            IAtomContainer gMol = ExtAtomContainerManipulator.cloneWithIDs(mol);
            if (DEBUG) {
                System.out.println("standardize reaction module phase 2.1.1");
                System.out.println("t_mol " + SmilesGenerator.unique().create(gMol));
            }
            ExtAtomContainerManipulator.fixDativeBonds(gMol);
            if (DEBUG) {
                System.out.println("standardize reaction module phase 2.1.2");
                System.out.println("t_mol " + SmilesGenerator.unique().create(gMol));
            }
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(gMol);
            IAtomContainer molWithH = gMol;
            //= ExtAtomContainerManipulator.addExplicitH(gMol);
            ExtAtomContainerManipulator.aromatizeMolecule(molWithH);

            if (DEBUG) {
                System.out.println("standardize reaction module phase 2.2");
                System.out.println("t_mol " + SmilesGenerator.unique().create(molWithH));
            }
            if (id == null) {
                molWithH = setProperty(molWithH);
            } else {
                molWithH.setID(id);
            }
            if (DEBUG) {
                System.out.println("standardize reaction module phase 2.3");
                System.out.println("standardize t_mol " + SmilesGenerator.unique().create(molWithH));
            }
            if (stoichiometryMap.containsKey(molWithH.getID())) {
                tempStoic += stoichiometryMap.get(molWithH.getID());
                stoichiometryMap.put(molWithH.getID(), tempStoic);
//                System.err.println("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            Collections.sort(_metabolites, comparator);
        } catch (Exception ex) {
            Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }

        setProductMolecule(standardizedReaction, _metabolites);
        _metabolites.clear();
        //As per IntEnz 0 for undefined direction, 1 for LR, 2 for RL and 3 for bidirectional
        //As per CDK BIDIRECTION 1, Forward 2, Backward 0

        reactionSet.addReaction(standardizedReaction);

        //BIDIRECTION 1, Forward 2, Backward 0
        if (reaction.getDirection() != null) {
            standardizedReaction.setDirection(reaction.getDirection());
        } else {
            standardizedReaction.setDirection(Direction.BIDIRECTIONAL);
        }
        fingerprintMap.clear();
        moleculeMap.clear();
        stoichiometryMap.clear();

        if (DEBUG) {
            String postCreateReactionSMILES = SmilesGenerator.unique().aromatic().createReactionSMILES(standardizedReaction);
            System.out.println("post CreateReactionSMILES " + postCreateReactionSMILES);
        }

        if (DEBUG) {
            System.out.println("standardize reaction module end");
        }
        return standardizedReaction;
    }

    private IAtomContainer setProperty(IAtomContainer molecule) throws Exception {
        /*
         If ID is NULL or empty please assign it to null
         */
        String molID = molecule.getID() == null || molecule.getID().isEmpty() ? null : molecule.getID();
        if (DEBUG) {
            System.err.println("BEFORE");
            printAtoms(molecule);
        }
        try {
            try {
                if (molecule.getAtomCount() > 0) {
                    IFingerprintGenerator fpr = new FingerprintGenerator();
                    BitSet fingerprint_Present_Mol = fpr.getFingerprint(molecule);
                    //Loop for Unique Mol ID Creation
                    if (!fingerprint_Present_Mol.isEmpty()) {
                        if (!isValuePresent(fingerprint_Present_Mol)) {
                            if (molID == null) {
                                moleculeCounter += 1;
                                int val = moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
                            fingerprintMap.put(molID, fingerprint_Present_Mol);
                            moleculeMap.put(molID, molecule);
                        } else if (isValuePresent(fingerprint_Present_Mol)
                                && isAtomContainerPresent(getMoleculeID(fingerprint_Present_Mol), molecule)) {
                            if (molID == null) {
                                molID = getMoleculeID(fingerprint_Present_Mol);
                                molecule.setID(molID);
                            }
                        } else {
                            if (molID == null) {
                                moleculeCounter += 1;
                                int val = moleculeCounter + 100000;
                                String Temp = Integer.toString(val);
                                molID = Temp.replaceFirst("1", "M");
                                molecule.setID(molID);
                            }
                            fingerprintMap.put(molID, fingerprint_Present_Mol);
                            moleculeMap.put(molID, molecule);
                        }
                    } else {
                        System.err.println("error: Fingerprint can't be generated for this molecules");
                    }
                } else {
                    System.err.println("error: Mol file should contain atleast one atom!");
                }
            } catch (CDKException ex) {
                Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
            } catch (Exception ex) {
                Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
//            System.err.println("After");
//            printAtoms(molecule);
            if (molecule.getID() == null) {
                try {
                    throw new CDKException("Mol ID is NULL");
                } catch (CDKException ex) {
                    Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

        } catch (Exception ex) {
            Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (DEBUG) {
            System.out.println("After Cleanup Mol ID is " + molID);
//        printAtoms(molecule);
        }
        return molecule;
    }

    private void setReactantMolecule(IReaction IR, Collection<IAtomContainer> metabolites) {
//
//        System.out.println("-------------------");
//        System.out.println("Reactants");
//        System.out.println("-------------------");

        Iterator<IAtomContainer> it = metabolites.iterator();
        //System.out.println("Stoic Map Size: " + stoichiometryMap.size());

        while (it.hasNext()) {
            IAtomContainer mol = it.next();
//            System.out.println("Insertion: " + mol.getID() + ", " + stoichiometryMap.get(mol.getID()));
            mol.setProperty("STOICHIOMETRY", stoichiometryMap.get(mol.getID()));
            IR.addReactant(mol, stoichiometryMap.get(mol.getID()));
            //IR.setReactantCoefficient(mol,wt);
//            System.out.println("Reaction: " + mol.getID() + ", " + standardizedReaction.getProductCoefficient(mol));
        }

        metabolites.clear();
        stoichiometryMap.clear();
    }

    private void setProductMolecule(IReaction IR, Collection<IAtomContainer> metabolites) {

//        System.out.println("-------------------");
//        System.out.println("Products");
//        System.out.println("-------------------");
        Iterator<IAtomContainer> it = metabolites.iterator();
        while (it.hasNext()) {
            IAtomContainer mol = it.next();
//
//            System.out.println("Insertion: " + mol.getID() + ", " + stoichiometryMap.get(mol.getID()));
            mol.setProperty("STOICHIOMETRY", stoichiometryMap.get(mol.getID()));
            IR.addProduct(mol, stoichiometryMap.get(mol.getID()));
//         standardizedReaction.setProductCoefficient(mol,stoichiometryMap.get(mol.getID()));
        }

        metabolites.clear();
        stoichiometryMap.clear();
    }

    /**
     *
     * @param value
     * @throws java.io.IOException
     * @return
     */
    private boolean isValuePresent(BitSet value) throws IOException {
        for (BitSet bitset : fingerprintMap.values()) {
            try {
                if (Similarity.getTanimotoSimilarity(value, bitset) == 1.0) {
                    return true;
                }
            } catch (Exception ex) {
                Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return false;
    }

    /**
     *
     * @param bitset
     * @return
     * @throws java.io.IOException
     */
    private String getMoleculeID(BitSet bitset) throws IOException {
        String Key = null;
        for (Map.Entry<String, BitSet> map : fingerprintMap.entrySet()) {
            String key = map.getKey();
            try {
                if (Similarity.getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                    Key = key;
                    break;
                }
            } catch (Exception ex) {
                Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        //System.err.println("Error: Unable to Find AtomContainer ID!!!");
        return Key;
    }

    /**
     *
     * @param key
     * @param molecule
     * @return
     * @throws Exception
     */
    private boolean isAtomContainerPresent(String key, IAtomContainer molecule) throws Exception {
        try {
            boolean flag = moleculeMap.containsKey(key);
            if (flag && molecule.getAtomCount() > 0) {
                IAtomContainer molFromContainer = moleculeMap.get(key);
                return isIdentical(molecule, molFromContainer, true);
            }
        } catch (Exception ex) {
            Logger.getLogger(MolContainer.class.getName()).log(Level.SEVERE, null, ex);
        }
        return false;
    }

    /**
     *
     * @param queryMol_org
     * @param targetMol_org
     * @param removeHydrogen
     * @return
     * @throws Exception
     */
    private boolean isIdentical(IAtomContainer queryMol_org, IAtomContainer targetMol_org, boolean removeHydrogen) throws Exception {

        IAtomContainer queryMol = queryMol_org.clone();
        IAtomContainer targetMol = targetMol_org.clone();

        if (removeHydrogen) {
            queryMol = ExtAtomContainerManipulator.removeHydrogens(queryMol);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(queryMol);
            ExtAtomContainerManipulator.aromatizeMolecule(queryMol);
            targetMol = ExtAtomContainerManipulator.removeHydrogens(targetMol);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(targetMol);
            ExtAtomContainerManipulator.aromatizeMolecule(targetMol);
        }

        if (queryMol.getAtomCount() == 1 && targetMol.getAtomCount() == 1) {
            IAtom a = queryMol.atoms().iterator().next();
            IAtom b = targetMol.atoms().iterator().next();
            return a.getSymbol().equalsIgnoreCase(b.getSymbol())
                    && Objects.equals(a.getFormalCharge(), b.getFormalCharge())
                    && queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount();
        }
        Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
        Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();

        int leftHandAtomCount = 0;

        for (IAtom a : queryMol.atoms()) {
            if (a.getSymbol().equals("H")) {
                continue;
            }
            if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                atomUniqueCounter1.put(a.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                atomUniqueCounter1.put(a.getSymbol(), counter);
            }
            leftHandAtomCount++;
        }
        if (DEBUG) {
            try {
                System.out.println("Q=mol " + SmilesGenerator.generic().create(queryMol));
            } catch (CDKException ex) {
                Logger.getLogger(CDKReactionBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        int rightHandAtomCount = 0;

        for (IAtom b : targetMol.atoms()) {
            if (b.getSymbol().equals("H")) {
                continue;
            }
            if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                atomUniqueCounter2.put(b.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                atomUniqueCounter2.put(b.getSymbol(), counter);
            }
            rightHandAtomCount++;
        }
        if (DEBUG) {
            try {
                System.out.println("T=mol " + SmilesGenerator.generic().create(targetMol));
            } catch (CDKException ex) {
                Logger.getLogger(ReactionMechanismTool.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        if (DEBUG) {
            System.out.println("atomUniqueCounter1 " + leftHandAtomCount);
            System.out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }

        if (leftHandAtomCount != rightHandAtomCount) {
            if (DEBUG) {
                System.err.println();
                System.err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                System.err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            }
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            if (DEBUG) {
                System.err.println();
                System.err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                System.err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            }
            return false;
        }

        if (DEBUG) {
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            boolean flag = queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount();
            System.out.println("is Radical same " + flag);
        }

        return atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())
                ? queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount()
                        ? isSubgraphIdentical(queryMol, targetMol, removeHydrogen) : false : false;
    }

    private boolean isSubgraphIdentical(IAtomContainer qMol, IAtomContainer tMol, boolean removeHydrogen) throws CDKException, IOException, CloneNotSupportedException {
//        System.out.println("Graph matching");

        IAtomContainer mol1 = qMol.clone();
        IAtomContainer mol2 = tMol.clone();

        if (removeHydrogen) {
            mol1 = ExtAtomContainerManipulator.removeHydrogens(mol1);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
            ExtAtomContainerManipulator.aromatizeMolecule(mol1);
            mol2 = ExtAtomContainerManipulator.removeHydrogens(mol2);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
            ExtAtomContainerManipulator.aromatizeMolecule(mol2);
        }
        if (mol1.getAtomCount() != mol2.getAtomCount()) {
            return false;
        }
        Substructure mcs = new Substructure(mol1, mol2, true, true, true, false);
        mcs.setChemFilters(true, true, true);
        return mcs.isSubgraph() && !mcs.isStereoMisMatch();
    }
}
