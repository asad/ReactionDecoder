/*
 * Copyright (C) 2003-2020 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.graph;

import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.getProperty;
import static java.lang.System.nanoTime;
import java.util.ArrayList;
import java.util.Arrays;
import static java.util.Collections.sort;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.mapping.cache.ThreadSafeCache;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MCSThread implements Callable<MCSSolution> {

    private static final boolean DEBUG1 = false;
    private static final boolean DEBUG2 = false;
    private static final boolean DEBUG3 = false;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MCSThread.class);

    private final SmilesGenerator smiles;

    static final String NEW_LINE = getProperty("line.separator");

    /**
     *
     */
    protected final IAtomContainer compound1;

    /**
     *
     */
    protected final IAtomContainer compound2;

    /**
     *
     */
    protected final int queryPosition;

    /**
     *
     */
    protected final int targetPosition;

    /**
     *
     */
    protected final IMappingAlgorithm theory;

    /**
     *
     */
    long startTime;
    private boolean hasRings;
    private int numberOfCyclesEduct;
    private int numberOfCyclesProduct;

    /**
     *
     * @param theory
     * @param queryPosition
     * @param targetPosition
     * @param educt
     * @param product
     * @param bondMatcher
     * @param ringMatcher
     * @param atomMatcher
     * @throws CloneNotSupportedException
     * @throws org.openscience.cdk.exception.CDKException
     */
    MCSThread(IMappingAlgorithm theory, int queryPosition, int targetPosition,
            IAtomContainer educt, IAtomContainer product)
            throws CloneNotSupportedException, CDKException {
        /*
         * create SMILES
         */
        smiles = new SmilesGenerator(
                //SmiFlavor.Unique|
                SmiFlavor.UseAromaticSymbols
                | SmiFlavor.Stereo
                | SmiFlavor.AtomAtomMap);

        this.compound1 = getNewContainerWithIDs(educt);
        this.compound2 = getNewContainerWithIDs(product);
        this.queryPosition = queryPosition;
        this.targetPosition = targetPosition;
        this.theory = theory;
        this.numberOfCyclesEduct = 0;
        this.numberOfCyclesProduct = 0;
    }

    synchronized void printMatch(BaseMapping isomorphism) {
        int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                : isomorphism.getFirstAtomMapping().getCount();

        try {
            System.out.println("Q: " + isomorphism.getQuery().getID()
                    + " T: " + isomorphism.getTarget().getID()
                    + " atoms: " + isomorphism.getQuery().getAtomCount()
                    + " atoms: " + isomorphism.getTarget().getAtomCount()
                    + " overlaps: " + overlap
                    + " mcs " + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
        } catch (CloneNotSupportedException | CDKException ex) {
            LOGGER.error(Level.SEVERE, "Print MCS ", ex.getMessage());
        }
    }

    @Override
    public synchronized MCSSolution call() throws Exception {
        boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

        AtomMatcher am;
        BondMatcher bm;

        if (DEBUG1) {
            System.out.println("in mcsthread call ");
        }

        try {
            /*
             * IMP: Do not perform substructure matching for disconnected molecules
             */
            boolean moleculeConnected = isMoleculeConnected(getCompound1(), getCompound2());
            // boolean moleculeConnected = true;
            /*
                 Check if MCS matching required or not very IMP step
             */
            boolean possibleVFmatch12 = isPossibleSubgraphMatch(getCompound1(), getCompound2());
            if (DEBUG1) {
                System.out.println("VF Matcher 1->2 " + possibleVFmatch12);
            }

            boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
            if (DEBUG1) {
                System.out.println("VF Matcher 2->1 " + possibleVFmatch21);
            }

            if (moleculeConnected && possibleVFmatch12) {
                if (DEBUG1) {
                    System.out.println("Substructure 1");
                    this.startTime = currentTimeMillis();

                }
                IAtomContainer ac1 = duplicate(getCompound1());
                IAtomContainer ac2 = duplicate(getCompound2());

                if (DEBUG1) {
                    System.out.println("---1.1---");
                }
                Substructure substructure;
                am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                substructure = new Substructure(ac1, ac2, am, bm, true);

                if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                    am = AtomBondMatcher.atomMatcher(false, ringFlag);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    if (DEBUG1) {
                        System.out.println("---1.3---");
                    }
                    substructure = new Substructure(ac1, ac2,
                            am, bm, true);
                } else if (moleculeConnected && !substructure.isSubgraph()) {
                    am = AtomBondMatcher.atomMatcher(false, false);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    if (DEBUG1) {
                        System.out.println("---1.2---");
                    }
                    substructure = new Substructure(ac1, ac2, am, bm, true);
                }
                substructure.setChemFilters(true, true, true);
//                    System.out.println("Number of Solutions: " + substructure.getAllAtomMapping());
                if (substructure.isSubgraph()
                        && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
                    if (DEBUG1) {
                        System.out.println("Found Substructure 1");
                    }
                    MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                            substructure.getQuery(), substructure.getTarget(), substructure.getFirstAtomMapping());
                    mcs.setEnergy(substructure.getEnergyScore(0));
                    mcs.setFragmentSize(substructure.getFragmentSize(0));
                    mcs.setStereoScore(substructure.getStereoScore(0));
                    if (DEBUG1) {
                        long stopTime = currentTimeMillis();
                        long time = stopTime - startTime;
                        printMatch(substructure);
                        if (DEBUG1) {
                            System.out.println("\" Time:\" " + time);
                        }
                    }
                    return mcs;
                } else if (DEBUG1) {
                    System.out.println("not a Substructure 1");
                }
            }

            if (moleculeConnected && !possibleVFmatch12 && possibleVFmatch21) {

                if (DEBUG1) {
                    System.out.println("Substructure 2");
                    this.startTime = currentTimeMillis();

                }

                IAtomContainer ac1 = duplicate(getCompound1());
                IAtomContainer ac2 = duplicate(getCompound2());
                Substructure substructure;

                if (DEBUG1) {
                    System.out.println("---2.1---");
                }
                am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                substructure = new Substructure(ac2, ac1, am, bm, true);

                if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                    am = AtomBondMatcher.atomMatcher(false, ringFlag);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    if (DEBUG1) {
                        System.out.println("---2.3---");
                    }
                    substructure = new Substructure(ac2, ac1, am, bm, true);
                } else if (moleculeConnected && !substructure.isSubgraph()) {
                    am = AtomBondMatcher.atomMatcher(false, false);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    if (DEBUG1) {
                        System.out.println("---2.2---");
                    }
                    substructure = new Substructure(ac2, ac1, am, bm, true);
                }
                substructure.setChemFilters(true, true, true);

                if (substructure.isSubgraph()
                        && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {

                    if (DEBUG1) {
                        System.out.println("Found Substructure 2");
                    }
                    AtomAtomMapping aam = new AtomAtomMapping(substructure.getTarget(), substructure.getQuery());
                    Map<IAtom, IAtom> mappings = substructure.getFirstAtomMapping().getMappingsByAtoms();
                    mappings.keySet().stream().forEach((atom1) -> {
                        IAtom atom2 = mappings.get(atom1);
                        aam.put(atom2, atom1);
                    });
                    MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                            substructure.getTarget(), substructure.getQuery(), aam);
                    mcs.setEnergy(substructure.getEnergyScore(0));
                    mcs.setFragmentSize(substructure.getFragmentSize(0));
                    mcs.setStereoScore(substructure.getStereoScore(0));

                    if (DEBUG1) {
                        long stopTime = currentTimeMillis();
                        long time = stopTime - startTime;
                        printMatch(substructure);
                        System.out.println("\" Time:\" " + time);
                    }
                    return mcs;
                } else if (DEBUG1) {
                    System.out.println("not a Substructure 2");
                }
            }

            /*
             * If substructure matches have failed then call MCS
             */
            if (DEBUG1) {
//                /*
//                 * create SMILES
//                 */
//                String createSM1 = null;
//                String createSM2 = null;
//                try {
//                    createSM1 = smiles.create(getCompound1());
//                    createSM2 = smiles.create(getCompound2());
//                } catch (Exception e) {
//                    LOGGER.error(SEVERE, "Error in generating SMILES: ", e.getMessage());
//                }
                System.out.println("==============================================");

                System.out.println("No Substructure found - switching to MCS");
                System.out.println("Q: " + getCompound1().getID()
                        //                        + " molQ: " + createSM1
                        + NEW_LINE
                        + " T: " + getCompound2().getID()
                        //                        + " molT: " + createSM2
                        + NEW_LINE
                        + " atomsE: " + compound1.getAtomCount()
                        + " atomsP: " + compound2.getAtomCount());
                System.out.println("==============================================");

            }
            if (DEBUG1) {
                this.startTime = currentTimeMillis();
            }
            MCSSolution mcs = mcs();
            if (DEBUG1) {
                long stopTime = currentTimeMillis();
                long time = stopTime - startTime;
                System.out.println("\"MCS Time:\" " + time);
            }
            return mcs;

        } catch (CDKException | CloneNotSupportedException ex) {
//            ex.printStackTrace();
            LOGGER.error(SEVERE, "Error in generating MCS Solution: ", ex.getMessage());
        }
        return null;
    }

    private synchronized IAtomContainer getNewContainerWithIDs(IAtomContainer mol)
            throws CDKException, CloneNotSupportedException {
        /*
         Generating SMILES speed ups the mapping process by n-FOLDS
         May be this is CDK inherent bug, which relies on the properties set by the SMILES
         */
        if (mol != null && mol.getAtomCount() > 0) {
            IAtomContainer ac;
            ac = ExtAtomContainerManipulator.cloneWithIDs(mol);
            Aromaticity aromaticity = new Aromaticity(daylight(),
                    Cycles.or(Cycles.all(),
                            Cycles.or(Cycles.relevant(),
                                    Cycles.essential())));
            aromaticity.apply(ac);
            try {
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
                MoleculeInitializer.initializeMolecule(ac);
            } catch (Exception ex) {
                if (DEBUG1) {
                    ex.printStackTrace();
                }
                LOGGER.error(Level.SEVERE, "WARNING: Error in Config. r.mol: ", ex.getMessage());
            }

            for (int i = 0; i < ac.getAtomCount(); i++) {
                String atomID = mol.getAtom(i).getID() == null
                        ? valueOf(i) : mol.getAtom(i).getID();
                ac.getAtom(i).setID(atomID);
            }
            String containerID = mol.getID() == null ? valueOf(nanoTime()) : mol.getID();
            ac.setID(containerID);

            return ac;
        }
        return mol;
    }

    private synchronized boolean isPossibleSubgraphMatch(IAtomContainer q, IAtomContainer t) {
        if (DEBUG1) {
            System.out.println("check isPossibleSubgraphMatch " + q.getID() + "," + t.getID());
        }
        Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
        Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();

        for (IAtom a : q.atoms()) {
            if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                atomUniqueCounter1.put(a.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                atomUniqueCounter1.put(a.getSymbol(), counter);
            }
        }

        for (IAtom b : t.atoms()) {
            if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                atomUniqueCounter2.put(b.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                atomUniqueCounter2.put(b.getSymbol(), counter);
            }
        }

        if (atomUniqueCounter1.size() > atomUniqueCounter2.size()) {
            return false;
        }

        List<String> difference = new LinkedList<>(atomUniqueCounter1.keySet());
        difference.removeAll(atomUniqueCounter2.keySet());

        if (DEBUG2) {
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1.size());
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2.size());
            System.out.println("diff " + difference.size());
        }

        if (difference.isEmpty()) {
            if (!atomUniqueCounter1.keySet().stream().noneMatch((k)
                    -> (atomUniqueCounter1.get(k) > atomUniqueCounter2.get(k)))) {
                return false;
            }
        }

        return difference.isEmpty();
    }

    private synchronized int expectedMaxGraphmatch(IAtomContainer q, IAtomContainer t) {

        /*
         a={c,c,c,o,n}
         b={c,c,c,p}
       
         expectedMaxGraphmatch=3;
         */
        List<String> atomUniqueCounter1 = new ArrayList<>();
        List<String> atomUniqueCounter2 = new ArrayList<>();

        for (IAtom a : q.atoms()) {
            String hyb = a.getHybridization() == UNSET
                    ? a.getSymbol() : a.getAtomTypeName();
            atomUniqueCounter1.add(hyb);
        }

        for (IAtom b : t.atoms()) {
            String hyb = b.getHybridization() == UNSET
                    ? b.getSymbol() : b.getAtomTypeName();
            atomUniqueCounter2.add(hyb);
        }

        sort(atomUniqueCounter1);
        sort(atomUniqueCounter2);

        if (atomUniqueCounter1.isEmpty()) {
            return 0;
        }
        List<String> common = new LinkedList<>(atomUniqueCounter1);
        common.retainAll(atomUniqueCounter2);

        if (DEBUG2) {
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            System.out.println("atomUniqueCounter1 " + atomUniqueCounter1.size());
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2.size());
            System.out.println("Common " + common.size());
        }
        atomUniqueCounter1.clear();
        atomUniqueCounter2.clear();
        return common.size();
    }

    synchronized MCSSolution mcs() throws CDKException, CloneNotSupportedException {

        if (DEBUG3) {
            System.out.println("=============MCS============");
        }
        /*
         * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
         */
        IAtomContainer ac1 = duplicate(getCompound1());
        IAtomContainer ac2 = duplicate(getCompound2());
        Isomorphism isomorphism;
        int expectedMaxGraphmatch = expectedMaxGraphmatch(ac1, ac2);
        //boolean moleculeConnected = isMoleculeConnected(ac1, ac2);
        boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

        if (DEBUG3) {
            System.out.println("Expected matches " + expectedMaxGraphmatch);
        }

        String key;
        MCSSolution mcs;
        AtomMatcher am;
        BondMatcher bm;
        boolean atomType;
        boolean bondMatch;
        boolean ringMatch;
        boolean ringSizeMatch;

        switch (theory) {
            case RINGS:

                atomType = false;
                bondMatch = false;
                ringMatch = ringFlag;
                ringSizeMatch = isHasPerfectRings();

                break;

            case MIN:

                atomType = false;
                bondMatch = false;
                ringMatch = isHasPerfectRings();
                ringSizeMatch = false;

                break;

            case MAX:

                atomType = false;
                bondMatch = true;
                ringMatch = isHasPerfectRings();
                ringSizeMatch = false;

                break;
            default:

                atomType = false;
                bondMatch = false;
                ringMatch = isHasPerfectRings();
                ringSizeMatch = false;

                break;
        }

        am = AtomBondMatcher.atomMatcher(atomType, ringSizeMatch);
        bm = AtomBondMatcher.bondMatcher(bondMatch, ringMatch);

        key = generateUniqueKey(getCompound1().getID(), getCompound2().getID(),
                compound1.getAtomCount(), compound2.getAtomCount(),
                compound1.getBondCount(), compound2.getBondCount(),
                atomType,
                bondMatch,
                ringMatch,
                ringSizeMatch,
                numberOfCyclesEduct,
                numberOfCyclesProduct
        );
        if (ThreadSafeCache.getInstance().containsKey(key)) {
            if (DEBUG3) {
                System.out.println("===={Aladdin} Mapping {Gini}====");
            }
            MCSSolution solution = (MCSSolution) ThreadSafeCache.getInstance().get(key);
            mcs = copyOldSolutionToNew(
                    getQueryPosition(), getTargetPosition(),
                    getCompound1(), getCompound2(),
                    solution);

        } else {
            isomorphism = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, am, bm);
            mcs = addMCSSolution(key, ThreadSafeCache.getInstance(), isomorphism);
        }

        //System.out.println("cache map size " + ThreadSafeCache.getInstance().keySet().size());
        return mcs;

    }

    private synchronized IAtomContainer duplicate(IAtomContainer ac) throws CloneNotSupportedException {
        IAtomContainer a = ac.clone();
        a.setID(ac.getID());

        for (int i = 0; i < a.getAtomCount(); i++) {
            a.getAtom(i).setID(ac.getAtom(i).getID());
        }

        ac.setProperties(ac.getProperties());
        ac.setFlags(ac.getFlags());
        ac.setID(ac.getID());
        ac.notifyChanged();

        return a;
    }

    /**
     * @return the compound1
     */
    synchronized IAtomContainer getCompound1() {
        return compound1;
    }

    /**
     * @return the compound2
     */
    synchronized IAtomContainer getCompound2() {
        return compound2;
    }

    /**
     * @return the queryPosition
     */
    synchronized int getQueryPosition() {
        return queryPosition;
    }

    /**
     * @return the targetPosition
     */
    synchronized int getTargetPosition() {
        return targetPosition;
    }

    synchronized void setHasPerfectRings(boolean ring) {
        this.hasRings = ring;
    }

    /**
     * @return the hasRings
     */
    synchronized boolean isHasPerfectRings() {
        return hasRings;
    }

    /*
     * Check if fragmented container has single atom
     */
    private synchronized boolean isMoleculeConnected(IAtomContainer compound1, IAtomContainer compound2) {
        if (DEBUG1) {
            System.out.println("isMoleculeConnected");
        }
        boolean connected1 = true;

        IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound1);
        if (partitionIntoMolecules.getAtomContainerCount() > 1) {
            connected1 = false;
        }

        boolean connected2 = true;

        partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound2);
        if (partitionIntoMolecules.getAtomContainerCount() > 1) {
            connected2 = false;
        }

        if (DEBUG1) {
            System.out.println("Flag 1 " + connected1);
        }
        if (DEBUG1) {
            System.out.println("Flag 2 " + connected2);
        }
        return connected1 & connected2;
    }

    synchronized void setEductRingCount(int numberOfCyclesEduct) {
        this.numberOfCyclesEduct = numberOfCyclesEduct;
    }

    synchronized void setProductRingCount(int numberOfCyclesProduct) {
        this.numberOfCyclesProduct = numberOfCyclesProduct;
    }

    synchronized String generateUniqueKey(
            String id1, String id2,
            int atomCount1, int atomCount2,
            int bondCount1, int bondCount2,
            boolean atomtypeMatcher,
            boolean bondMatcher,
            boolean ringMatcher,
            boolean hasPerfectRings,
            int numberOfCyclesEduct, int numberOfCyclesProduct) {

        //System.out.println("====generate Unique Key====");
        StringBuilder key = new StringBuilder();
        key.append(id1).append(id2)
                .append(atomCount1)
                .append(atomCount2)
                .append(bondCount1)
                .append(bondCount2)
                .append(atomtypeMatcher)
                .append(bondMatcher)
                .append(ringMatcher)
                .append(hasPerfectRings)
                .append(numberOfCyclesEduct)
                .append(numberOfCyclesProduct);

        try {
            try {
                int[] sm1 = getCircularFP(compound1);
                int[] sm2 = getCircularFP(compound2);
                key.append(Arrays.toString(sm1));
                key.append(Arrays.toString(sm2));
            } catch (Exception ex) {
                LOGGER.error(Level.SEVERE, "Error in Generating Circular FP: ", ex);
            }
            if (DEBUG1) {
                System.out.println("Unique KEY " + key);
            }
        } catch (Exception ex) {
            if (DEBUG1) {
                ex.printStackTrace();
            }
            LOGGER.error(SEVERE, "Error in generating Unique Key: ", ex.getMessage());
        }
        return key.toString();
    }

    private synchronized int[] getCircularFP(IAtomContainer mol) throws CDKException {
        CircularFingerprinter circularFingerprinter = new CircularFingerprinter(6, 1024);
        circularFingerprinter.setPerceiveStereo(true);
        IBitFingerprint bitFingerprint = circularFingerprinter.getBitFingerprint(mol);
        return bitFingerprint.getSetbits();
    }

    /*
     * copy old mapping from the cache to new
     */
    synchronized MCSSolution copyOldSolutionToNew(int queryPosition, int targetPosition,
            IAtomContainer compound1, IAtomContainer compound2, MCSSolution oldSolution) {
        AtomAtomMapping atomAtomMapping = oldSolution.getAtomAtomMapping();
        Map<Integer, Integer> mappingsByIndex = atomAtomMapping.getMappingsByIndex();

        AtomAtomMapping atomAtomMappingNew = new AtomAtomMapping(compound1, compound2);
        mappingsByIndex.entrySet().forEach((m) -> {
            atomAtomMappingNew.put(compound1.getAtom(m.getKey()), compound2.getAtom(m.getValue()));
        });
        MCSSolution mcsSolution = new MCSSolution(queryPosition, targetPosition, compound1, compound2, atomAtomMappingNew);
        mcsSolution.setEnergy(oldSolution.getEnergy());
        mcsSolution.setFragmentSize(oldSolution.getFragmentSize());
        mcsSolution.setStereoScore(oldSolution.getStereoScore());

        return mcsSolution;
    }

    synchronized MCSSolution addMCSSolution(String key, ThreadSafeCache<String, MCSSolution> mappingcache, Isomorphism isomorphism) {

        isomorphism.setChemFilters(true, true, true);
        if (DEBUG3) {
            try {
                System.out.println("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", "
                        + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
            } catch (CloneNotSupportedException | CDKException e) {
                LOGGER.error(SEVERE, "Error in computing MCS ", e.getMessage());
            }
        }
        /*
         * In case of Complete subgraph, don't use Energy filter
         *
         */
        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                isomorphism.getQuery(), isomorphism.getTarget(), isomorphism.getFirstAtomMapping());
        mcs.setEnergy(isomorphism.getEnergyScore(0));
        mcs.setFragmentSize(isomorphism.getFragmentSize(0));
        mcs.setStereoScore(isomorphism.getStereoScore(0));
        if (DEBUG1) {
            long stopTime = currentTimeMillis();
            long time = stopTime - startTime;
            printMatch(isomorphism);
            System.out.println("\" Time:\" " + time);
        }
        if (!mappingcache.containsKey(key)) {
            if (DEBUG3) {
                System.out.println("Key " + key);
                try {
                    System.out.println("mcs size " + mcs.getAtomAtomMapping().getCount());
                    System.out.println("mcs map " + mcs.getAtomAtomMapping().getMappingsByIndex());
                    System.out.println("mcs " + mcs.getAtomAtomMapping().getCommonFragmentAsSMILES());
                    System.out.println("\n\n\n ");
                } catch (CloneNotSupportedException | CDKException ex) {
                    LOGGER.error(SEVERE, "Unable to create SMILES ", ex.getMessage());
                }
            }
            mappingcache.put(key, mcs);
        }
        return mcs;
    }
}
