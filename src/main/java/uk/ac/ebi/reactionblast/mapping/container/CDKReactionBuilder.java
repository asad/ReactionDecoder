/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import static java.util.Collections.sort;
import static java.util.Collections.synchronizedMap;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.tools.AtomContainerSetComparator;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.fixDativeBonds;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogens;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CDKReactionBuilder extends BasicDebugger implements Serializable {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 19869866609698L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(CDKReactionBuilder.class);
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
        reactionSet = getInstance().newInstance(IReactionSet.class);
        stoichiometryMap = synchronizedMap(new HashMap<>());
        fingerprintMap = synchronizedMap(new HashMap<>());
        moleculeMap = synchronizedMap(new HashMap<>());
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
    public synchronized IReaction standardize(IReaction reaction) throws Exception {
        int old_atom_rank_index_reactant = 1;
        int old_atom_rank_index_product = 1;
        if (DEBUG) {
            SmilesGenerator createReactionSMILES = new SmilesGenerator(
                    SmiFlavor.Unique
                    | SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap);
            out.println("createReactionSMILES " + createReactionSMILES.create(reaction));
            out.println("standardize reaction module start");
        }

        List<IAtomContainer> _metabolites = new ArrayList<>();
        IReaction standardizedReaction = getInstance().newInstance(IReaction.class);

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
            out.println("standardize reaction module phase 1");
        }
        for (IAtomContainer mol : reaction.getReactants().atomContainers()) {
            String id = mol.getID() == null || mol.getID().isEmpty() ? null : mol.getID();
            tempStoic = 1.0;
            if (reaction.getReactantCoefficient(mol) > 0) {
                tempStoic = reaction.getReactantCoefficient(mol);
            }

            if (DEBUG) {
                out.println("q_mol " + unique().create(mol));
                out.println("standardize reaction module phase 1.1");
            }

            IAtomContainer gMol = cloneWithIDs(mol);

            /*
             * Set old Atom Index
             */
            for (IAtom a : gMol.atoms()) {
                if (a.getProperties() == null) {
                    a.addProperties(new HashMap<>());
                }
                a.setProperty("OLD_RANK", old_atom_rank_index_reactant++);
            }
            if (DEBUG) {
                out.println("standardize reaction module phase 1.1.1");
            }
            fixDativeBonds(gMol);
            if (DEBUG) {
                out.println("standardize reaction module phase 1.1.2");
            }
            percieveAtomTypesAndConfigureAtoms(gMol);
            IAtomContainer molWithH = gMol;
            //= ExtAtomContainerManipulator.addExplicitH(gMol);
            aromatizeMolecule(molWithH);

            if (DEBUG) {
                out.println("standardize reaction module phase 1.2");
            }

            if (id == null) {
                molWithH = setProperty(molWithH);
            } else {
                molWithH.setID(id);
            }

            if (DEBUG) {
                out.println("standardize reaction module phase 1.3");
                out.println("After Cleanup " + molWithH.getID());
                out.println(unique().create(mol));
            }
            if (stoichiometryMap.containsKey(molWithH.getID())) {
                tempStoic += stoichiometryMap.get(molWithH.getID());
                stoichiometryMap.put(molWithH.getID(), tempStoic);
//                System.LOGGER.debug("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            sort(_metabolites, comparator);
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }

        setReactantMolecule(standardizedReaction, _metabolites);
        _metabolites.clear();

        if (DEBUG) {
            out.println("standardize reaction module phase 2");
        }
        if (DEBUG) {
            out.println();
            out.println("****************************");
            out.println();
        }
        for (IAtomContainer mol : reaction.getProducts().atomContainers()) {
            String id = mol.getID() == null || mol.getID().isEmpty() ? null : mol.getID();
            tempStoic = 1.0;
            if (reaction.getProductCoefficient(mol) > 0) {
                tempStoic = reaction.getProductCoefficient(mol);
            }
            if (DEBUG) {
                out.println("standardize reaction module phase 2.1");
                out.println("t_mol " + unique().create(mol));
            }
            IAtomContainer gMol = cloneWithIDs(mol);

            /*
             * Set old Atom Index
             */
            for (IAtom a : gMol.atoms()) {
                if (a.getProperties() == null) {
                    a.addProperties(new HashMap<>());
                }
                a.setProperty("OLD_RANK", old_atom_rank_index_product++);
            }
            if (DEBUG) {
                out.println("standardize reaction module phase 2.1.1");
                out.println("t_mol " + unique().create(gMol));
            }
            fixDativeBonds(gMol);
            if (DEBUG) {
                out.println("standardize reaction module phase 2.1.2");
                out.println("t_mol " + unique().create(gMol));
            }
            percieveAtomTypesAndConfigureAtoms(gMol);
            IAtomContainer molWithH = gMol;
            //= ExtAtomContainerManipulator.addExplicitH(gMol);
            aromatizeMolecule(molWithH);

            if (DEBUG) {
                out.println("standardize reaction module phase 2.2");
                out.println("t_mol " + unique().create(molWithH));
            }
            if (id == null) {
                molWithH = setProperty(molWithH);
            } else {
                molWithH.setID(id);
            }
            if (DEBUG) {
                out.println("standardize reaction module phase 2.3");
                out.println("standardize t_mol " + unique().create(molWithH));
            }
            if (stoichiometryMap.containsKey(molWithH.getID())) {
                tempStoic += stoichiometryMap.get(molWithH.getID());
                stoichiometryMap.put(molWithH.getID(), tempStoic);
//                System.LOGGER.debug("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            sort(_metabolites, comparator);
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
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
            standardizedReaction.setDirection(BIDIRECTIONAL);
        }
        fingerprintMap.clear();
        moleculeMap.clear();
        stoichiometryMap.clear();

        if (DEBUG) {
            SmilesGenerator smiles = new SmilesGenerator(
                    SmiFlavor.Unique
                    | SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap);
            String postCreateReactionSMILES = smiles.create(standardizedReaction);
            out.println("post CreateReactionSMILES " + postCreateReactionSMILES);
        }

        if (DEBUG) {
            out.println("standardize reaction module end");
        }
        return standardizedReaction;
    }

    private IAtomContainer setProperty(IAtomContainer molecule) throws Exception {
        /*
         If ID is NULL or empty please assign it to null
         */
        String molID = molecule.getID() == null || molecule.getID().isEmpty() ? null : molecule.getID();
        if (DEBUG) {
            LOGGER.debug("BEFORE");
            printAtoms(molecule);
        }
        try {
            try {
                if (molecule.getAtomCount() > 0) {
                    IFingerprintGenerator fpr = new FingerprintGenerator();
                    BitSet fingerprint_Present_Mol = fpr.getFingerprint(molecule);
                    /*
                    Single Atom fingerprints
                     */
                    if (fingerprint_Present_Mol.isEmpty()) {
                        Fingerprinter fingerprinter = new Fingerprinter();
                        fingerprint_Present_Mol = fingerprinter.getBitFingerprint(molecule).asBitSet();
                    }
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
                        LOGGER.debug("error: Fingerprint can't be generated for this molecule " + SmilesGenerator.generic().create(molecule));
                    }
                } else {
                    LOGGER.debug("error: Mol file should contain atleast one atom! " + SmilesGenerator.generic().create(molecule));
                }
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
//            System.LOGGER.debug("After");
//            printAtoms(molecule);
            if (molecule.getID() == null) {
                try {
                    throw new CDKException("Mol ID is NULL");
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        if (DEBUG) {
            out.println("After Cleanup Mol ID is " + molID);
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
    private boolean isValuePresent(BitSet value) throws IOException, Exception {
        for (BitSet bitset : fingerprintMap.values()) {
            if (getTanimotoSimilarity(value, bitset) == 1.0) {
                return true;
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
                if (getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                    Key = key;
                    break;
                }
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        //System.LOGGER.debug("Error: Unable to Find AtomContainer ID!!!");
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
            LOGGER.error(SEVERE, null, ex);
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
            queryMol = removeHydrogens(queryMol);
            percieveAtomTypesAndConfigureAtoms(queryMol);
            aromatizeMolecule(queryMol);
            targetMol = removeHydrogens(targetMol);
            percieveAtomTypesAndConfigureAtoms(targetMol);
            aromatizeMolecule(targetMol);
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
                out.println("Q=mol " + generic().create(queryMol));
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
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
                out.println("T=mol " + generic().create(targetMol));
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + leftHandAtomCount);
            out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }

        if (leftHandAtomCount != rightHandAtomCount) {
            if (DEBUG) {
                LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            }
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            if (DEBUG) {
                LOGGER.warn("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                LOGGER.warn(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            }
            return false;
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            boolean flag = queryMol.getElectronContainerCount() == targetMol.getElectronContainerCount();
            out.println("is Radical same " + flag);
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
            mol1 = removeHydrogens(mol1);
            percieveAtomTypesAndConfigureAtoms(mol1);
            aromatizeMolecule(mol1);
            mol2 = removeHydrogens(mol2);
            percieveAtomTypesAndConfigureAtoms(mol2);
            aromatizeMolecule(mol2);
        }
        if (mol1.getAtomCount() != mol2.getAtomCount()) {
            return false;
        }
        Substructure mcs = new Substructure(mol1, mol2, true, true, true, false);
        mcs.setChemFilters(true, true, true);
        return mcs.isSubgraph() && !mcs.isStereoMisMatch();
    }
}
