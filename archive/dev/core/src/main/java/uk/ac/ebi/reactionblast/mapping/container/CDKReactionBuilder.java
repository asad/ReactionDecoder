/*
 * Copyright (C) 2003-2016 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.lang.System.err;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import static java.util.Collections.synchronizedMap;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.logging.Logger;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.containers.MolContainer;
import uk.ac.ebi.reactionblast.fingerprints.FingerprintGenerator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.AtomContainerSetComparator;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.fixDativeBonds;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogens;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Arrays.sort;
import static java.util.Collections.sort;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import uk.ac.ebi.reactionblast.mapping.helper.MappingHandler;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CDKReactionBuilder extends BasicDebugger implements Serializable {

    private final static boolean DEBUG = false;
    private static final long serialVersionUID = 19869866609698L;
    private static final Logger LOG = getLogger(CDKReactionBuilder.class.getName());
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
        stoichiometryMap = synchronizedMap(new HashMap<String, Double>());
        fingerprintMap = synchronizedMap(new HashMap<String, BitSet>());
        moleculeMap = synchronizedMap(new HashMap<String, IAtomContainer>());
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
        if (DEBUG) {
            String createReactionSMILES = generic().aromatic().createReactionSMILES(reaction);
            out.println("createReactionSMILES " + createReactionSMILES);
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
            int old_atom_rank_index_reactant = 1;
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
//                System.err.println("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            sort(_metabolites, comparator);
        } catch (Exception ex) {
            getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
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
            int old_atom_rank_index_product = 1;
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
//                System.err.println("St Map put: " + mol.getID() + ", St: " + tempStoic);
            } else {
                stoichiometryMap.put(molWithH.getID(), tempStoic);
                _metabolites.add(molWithH);
            }
        }

        try {
            Comparator<IAtomContainer> comparator = new AtomContainerSetComparator();
            sort(_metabolites, comparator);
        } catch (Exception ex) {
            getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
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
            String postCreateReactionSMILES = generic().aromatic().createReactionSMILES(standardizedReaction);
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
            err.println("BEFORE");
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
                        err.println("error: Fingerprint can't be generated for this molecule " + SmilesGenerator.generic().create(molecule));
                    }
                } else {
                    err.println("error: Mol file should contain atleast one atom! " + SmilesGenerator.generic().create(molecule));
                }
            } catch (CDKException ex) {
                getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
            } catch (Exception ex) {
                getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
            }
//            System.err.println("After");
//            printAtoms(molecule);
            if (molecule.getID() == null) {
                try {
                    throw new CDKException("Mol ID is NULL");
                } catch (CDKException ex) {
                    getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
                }
            }

        } catch (Exception ex) {
            getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
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
    private String getMoleculeID(BitSet bitset) throws IOException, Exception {
        for (Map.Entry<String, BitSet> map : fingerprintMap.entrySet()) {
            if (getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                return map.getKey();
            }
        }
        //System.err.println("Error: Unable to Find AtomContainer ID!!!");
        return null;
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
            getLogger(MolContainer.class.getName()).log(SEVERE, null, ex);
        }
        return false;
    }

    /**
     *
     * @param originalQMol
     * @param originalTMol
     * @param removeHydrogen
     * @return
     * @throws Exception
     */
    private boolean isIdentical(IAtomContainer originalQMol, IAtomContainer originalTMol, boolean removeHydrogen) throws Exception {

        IAtomContainer queryMol = originalQMol.clone();
        IAtomContainer targetMol = originalTMol.clone();

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
                getLogger(CDKReactionBuilder.class.getName()).log(SEVERE, null, ex);
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
                getLogger(ReactionMechanismTool.class.getName()).log(SEVERE, null, ex);
            }
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + leftHandAtomCount);
            out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }

        if (leftHandAtomCount != rightHandAtomCount) {
            if (DEBUG) {
                err.println();
                err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
            }
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            if (DEBUG) {
                err.println();
                err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                        + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
                err.println(atomUniqueCounter1 + " =/= " + atomUniqueCounter2);
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

    public static IReaction preprocessStandardizedReaction(IReaction referenceReaction) throws CDKException, IOException, Exception {
        IReaction reactionWithSTOICHIOMETRY = referenceReaction.getBuilder().newInstance(IReaction.class);
        try {
            for (int i = 0; i < referenceReaction.getReactantCount(); i++) {
                IAtomContainer refMol = referenceReaction.getReactants().getAtomContainer(i);
                IAtomContainer mol = cloneWithIDs(refMol);//refMol.clone();
                mol = canonicalLabelling(mol);

                mol.setID(referenceReaction.getReactants().getAtomContainer(i).getID());
                Double st = referenceReaction.getReactantCoefficient(refMol);
                aromatizeMolecule(mol);
                reactionWithSTOICHIOMETRY.addReactant(mol, st);
                /*
                 * Note: This is a work around to help calculate CIPLigand rule
                 * cip.rules.MassNumberRule.getMassNumber
                 */
                for (IAtom atom : mol.atoms()) {
                    int number = atom.getMassNumber() == null ? 11 : atom.getMassNumber();
                    atom.setMassNumber(number);
                }
            }
        } catch (CloneNotSupportedException | CDKException e) {
            getLogger(Reactor.class.getName()).log(SEVERE, null, e);
        }
        try {
            for (int i = 0; i < referenceReaction.getProductCount(); i++) {
                IAtomContainer refMol = referenceReaction.getProducts().getAtomContainer(i);
                IAtomContainer mol = cloneWithIDs(refMol);//refMol.clone();
                mol = canonicalLabelling(mol);

                mol.setID(referenceReaction.getProducts().getAtomContainer(i).getID());
                Double st = referenceReaction.getProductCoefficient(refMol);
                aromatizeMolecule(mol);
                reactionWithSTOICHIOMETRY.addProduct(mol, st);
                /*
                 * Note: This is a work around to help calculate CIPLigand rule
                 * cip.rules.MassNumberRule.getMassNumber
                 */
                for (IAtom atom : mol.atoms()) {
                    int number = atom.getMassNumber() == null ? 11 : atom.getMassNumber();
                    atom.setMassNumber(number);
                }
            }
            reactionWithSTOICHIOMETRY.setID(referenceReaction.getID());
            reactionWithSTOICHIOMETRY.setDirection(referenceReaction.getDirection());
        } catch (CloneNotSupportedException | CDKException e) {
            getLogger(Reactor.class.getName()).log(SEVERE, null, e);
        }

        /*
         * Clean mapping flag and bond change flags
         */
        MappingHandler.cleanMapping(reactionWithSTOICHIOMETRY);
        return expandReaction(reactionWithSTOICHIOMETRY);
    }

    private static IReaction expandReaction(IReaction reactionWithSTOICHIOMETRY) throws CloneNotSupportedException {
        IReaction reactionWithUniqueSTOICHIOMETRY = reactionWithSTOICHIOMETRY.getBuilder().newInstance(IReaction.class);
        for (int i = 0; i < reactionWithSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer _react = reactionWithSTOICHIOMETRY.getReactants().getAtomContainer(i);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getReactantCoefficient(_react);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer _reactDup = cloneWithIDs(_react);
                _reactDup.setID(_react.getID());
                _reactDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addReactant(_reactDup, 1.0);
            }
        }

        for (int j = 0; j < reactionWithSTOICHIOMETRY.getProductCount(); j++) {

            IAtomContainer _prod = reactionWithSTOICHIOMETRY.getProducts().getAtomContainer(j);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getProductCoefficient(_prod);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer prodDup = cloneWithIDs(_prod);
                prodDup.setID(_prod.getID());
                prodDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addProduct(prodDup, 1.0);
            }

        }

        reactionWithUniqueSTOICHIOMETRY.setID(
                reactionWithSTOICHIOMETRY.getID() == null
                        ? "MappedReaction (ecBLAST)"
                        : reactionWithSTOICHIOMETRY.getID());
        reactionWithUniqueSTOICHIOMETRY.setDirection(reactionWithSTOICHIOMETRY.getDirection() == null
                ? BIDIRECTIONAL
                : reactionWithSTOICHIOMETRY.getDirection());
//
//        System.out.println("ExpandedEduct Count: " + reactionWithUniqueSTOICHIOMETRY.getReactantCount()
//                + ", ExpandedProduct Count: " + reactionWithUniqueSTOICHIOMETRY.getProductCount());

        LabelAtoms(reactionWithUniqueSTOICHIOMETRY);

        return reactionWithUniqueSTOICHIOMETRY;
    }

    private static void LabelAtoms(IReaction reactionWithUniqueSTOICHIOMETRY) {
        int new_atom_rank_index_reactant = 1;
        int new_atom_rank_index_product = 1;
        Integer substrateAtomCounter = 1;
        Integer productAtomCounter = 1;

        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (substrateAtomCounter).toString();
                substrateAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
//                System.out.println("EAtom: " + k + " " + atom.getSymbol() + " Rank Atom: " + atom.getProperty("OLD_RANK") + " " + " Id: " + atom.getID());
                if (atom.getProperty("OLD_RANK") != null) {
                    atom.setProperty("NEW_RANK", new_atom_rank_index_reactant++);
                }
            }
        }

        for (int j = 0; j < reactionWithUniqueSTOICHIOMETRY.getProductCount(); j++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(j);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (productAtomCounter).toString();
                productAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
//                System.out.println("PAtom: " + k + " " + atom.getSymbol() + " Id: " + atom.getID());
                if (atom.getProperty("OLD_RANK") != null) {
                    atom.setProperty("NEW_RANK", new_atom_rank_index_product++);
                }
            }
        }

    }

    /**
     * Generates customised canonical labeling
     *
     * @param mol
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    private static IAtomContainer canonicalLabelling(IAtomContainer mol) throws CloneNotSupportedException, CDKException {

        IAtomContainer cloneMolecule = cloneWithIDs(mol);


        /*
        Use the Canonical labelling from the SMILES
        IMP: Suggested by John May
         */
        int[] p = new int[cloneMolecule.getAtomCount()];

        try {
            String smiles = unique().create(cloneMolecule, p);
            if (DEBUG) {
                err.println("smiles " + smiles);
            }
        } catch (CDKException e) {
            getLogger(Reactor.class.getName()).log(SEVERE, null, e);
        }

        permuteWithoutClone(p, cloneMolecule);


        /*
        Generate 2D Diagram without cloning
         */
        if (!has2DCoordinates(cloneMolecule)) {
            try {
                /*
                Clone it else it will loose mol ID
                 */
                StructureDiagramGenerator sdg = new StructureDiagramGenerator();
                sdg.setMolecule(cloneMolecule, false);
                sdg.generateCoordinates();
            } catch (CDKException e) {
            }
        }

        /*
        Set the IDs to -1 very IMP
         */
        for (IAtom atom : cloneMolecule.atoms()) {
            atom.setID("-1");
        }

        /*
        Set the IDs to container
         */
        if (mol.getID() != null) {
            cloneMolecule.setID(mol.getID());
        }

        return cloneMolecule;
    }

    /*
     * Generates customised permutation based canonical labeling of atoms and bonds
     * The idea is to canonicalise the atoms and bonds
     */
    public static void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
        int n = atomContainer.getAtomCount();
        if (DEBUG) {
            err.println("permuting " + java.util.Arrays.toString(p));
        }
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
}
