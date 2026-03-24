/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping.graph;

import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.getProperty;
import static java.lang.System.nanoTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
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
import com.bioinceptionlabs.reactionblast.mapping.cache.ThreadSafeCache;
import com.bioinceptionlabs.reactionblast.mapping.interfaces.IMappingAlgorithm;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class MCSThread implements Callable<MCSSolution> {

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

    void printMatch(BaseMapping isomorphism) {
        int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                : isomorphism.getFirstAtomMapping().getCount();

        try {
            LOGGER.debug("Q: " + isomorphism.getQuery().getID()
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
    public MCSSolution call() throws Exception {
        boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

        AtomMatcher am;
        BondMatcher bm;

        LOGGER.debug("in mcsthread call ");

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
            LOGGER.debug("VF Matcher 1->2 " + possibleVFmatch12);

            boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
            LOGGER.debug("VF Matcher 2->1 " + possibleVFmatch21);

            if (moleculeConnected && possibleVFmatch12) {
                LOGGER.debug("Substructure 1");
                this.startTime = currentTimeMillis();

                IAtomContainer ac1 = duplicate(getCompound1());
                IAtomContainer ac2 = duplicate(getCompound2());

                LOGGER.debug("---1.1---");
                Substructure substructure;
                am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                substructure = new Substructure(ac1, ac2, am, bm, true);

                if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                    am = AtomBondMatcher.atomMatcher(false, ringFlag);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    LOGGER.debug("---1.3---");
                    substructure = new Substructure(ac1, ac2,
                            am, bm, true);
                } else if (moleculeConnected && !substructure.isSubgraph()) {
                    am = AtomBondMatcher.atomMatcher(false, false);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    LOGGER.debug("---1.2---");
                    substructure = new Substructure(ac1, ac2, am, bm, true);
                }
                substructure.setChemFilters(true, true, true);
//                    System.out.println("Number of Solutions: " + substructure.getAllAtomMapping());
                if (substructure.isSubgraph()
                        && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
                    LOGGER.debug("Found Substructure 1");
                    MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                            substructure.getQuery(), substructure.getTarget(), substructure.getFirstAtomMapping());
                    mcs.setEnergy(substructure.getEnergyScore(0));
                    mcs.setFragmentSize(substructure.getFragmentSize(0));
                    mcs.setStereoScore(substructure.getStereoScore(0));
                    long stopTime = currentTimeMillis();
                    long time = stopTime - startTime;
                    printMatch(substructure);
                    LOGGER.debug("\" Time:\" " + time);
                    return mcs;
                } else {
                    LOGGER.debug("not a Substructure 1");
                }
            }

            if (moleculeConnected && !possibleVFmatch12 && possibleVFmatch21) {

                LOGGER.debug("Substructure 2");
                this.startTime = currentTimeMillis();

                IAtomContainer ac1 = duplicate(getCompound1());
                IAtomContainer ac2 = duplicate(getCompound2());
                Substructure substructure;

                LOGGER.debug("---2.1---");
                am = AtomBondMatcher.atomMatcher(true, isHasPerfectRings());
                bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                substructure = new Substructure(ac2, ac1, am, bm, true);

                if (!substructure.isSubgraph() && !theory.equals(IMappingAlgorithm.RINGS)) {
                    am = AtomBondMatcher.atomMatcher(false, ringFlag);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    LOGGER.debug("---2.3---");
                    substructure = new Substructure(ac2, ac1, am, bm, true);
                } else if (moleculeConnected && !substructure.isSubgraph()) {
                    am = AtomBondMatcher.atomMatcher(false, false);
                    bm = AtomBondMatcher.bondMatcher(false, isHasPerfectRings());

                    LOGGER.debug("---2.2---");
                    substructure = new Substructure(ac2, ac1, am, bm, true);
                }
                substructure.setChemFilters(true, true, true);

                if (substructure.isSubgraph()
                        && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {

                    LOGGER.debug("Found Substructure 2");
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

                    long stopTime = currentTimeMillis();
                    long time = stopTime - startTime;
                    printMatch(substructure);
                    LOGGER.debug("\" Time:\" " + time);
                    return mcs;
                } else {
                    LOGGER.debug("not a Substructure 2");
                }
            }

            /*
             * If substructure matches have failed then call MCS
             */
            LOGGER.debug("==============================================");
            LOGGER.debug("No Substructure found - switching to MCS");
            LOGGER.debug("Q: " + getCompound1().getID()
                    + NEW_LINE
                    + " T: " + getCompound2().getID()
                    + NEW_LINE
                    + " atomsE: " + compound1.getAtomCount()
                    + " atomsP: " + compound2.getAtomCount());
            LOGGER.debug("==============================================");
            this.startTime = currentTimeMillis();
            MCSSolution mcs = mcs();
            long stopTime = currentTimeMillis();
            long time = stopTime - startTime;
            LOGGER.debug("\"MCS Time:\" " + time);
            return mcs;

        } catch (CDKException | CloneNotSupportedException ex) {
//            ex.printStackTrace();
            LOGGER.error(SEVERE, "Error in generating MCS Solution: ", ex.getMessage());
        }
        return null;
    }

    /**
     * Clone molecule preserving IDs.
     * Aromaticity and atom-type perception already done in GraphMatcher
     * before MCSThread is created — do NOT repeat here (was ~25% of total time).
     */
    private IAtomContainer getNewContainerWithIDs(IAtomContainer mol)
            throws CDKException, CloneNotSupportedException {
        if (mol != null && mol.getAtomCount() > 0) {
            IAtomContainer ac = ExtAtomContainerManipulator.cloneWithIDs(mol);

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

    private boolean isPossibleSubgraphMatch(IAtomContainer q, IAtomContainer t) {
        LOGGER.debug("check isPossibleSubgraphMatch " + q.getID() + "," + t.getID());
        Map<String, Integer> atomCount1 = new HashMap<>();
        Map<String, Integer> atomCount2 = new HashMap<>();

        for (IAtom a : q.atoms()) {
            atomCount1.merge(a.getSymbol(), 1, Integer::sum);
        }

        for (IAtom b : t.atoms()) {
            atomCount2.merge(b.getSymbol(), 1, Integer::sum);
        }

        if (atomCount1.size() > atomCount2.size()) {
            return false;
        }

        // Check all atom types in query exist in target with sufficient count
        for (Map.Entry<String, Integer> entry : atomCount1.entrySet()) {
            Integer targetCount = atomCount2.get(entry.getKey());
            if (targetCount == null || entry.getValue() > targetCount) {
                return false;
            }
        }

        return true;
    }

    private int expectedMaxGraphmatch(IAtomContainer q, IAtomContainer t) {
        /*
         a={c,c,c,o,n}
         b={c,c,c,p}
         expectedMaxGraphmatch=3;
         */
        Map<String, Integer> countQ = new HashMap<>();
        Map<String, Integer> countT = new HashMap<>();

        for (IAtom a : q.atoms()) {
            String hyb = a.getHybridization() == UNSET
                    ? a.getSymbol() : a.getAtomTypeName();
            countQ.merge(hyb, 1, Integer::sum);
        }

        for (IAtom b : t.atoms()) {
            String hyb = b.getHybridization() == UNSET
                    ? b.getSymbol() : b.getAtomTypeName();
            countT.merge(hyb, 1, Integer::sum);
        }

        if (countQ.isEmpty()) {
            return 0;
        }

        // Multiset intersection: min of counts for each common type
        int common = 0;
        for (Map.Entry<String, Integer> entry : countQ.entrySet()) {
            Integer tCount = countT.get(entry.getKey());
            if (tCount != null) {
                common += Math.min(entry.getValue(), tCount);
            }
        }

        return common;
    }

    MCSSolution mcs() throws CDKException, CloneNotSupportedException {

        LOGGER.debug("=============MCS============");
        /*
         * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
         */
        IAtomContainer ac1 = duplicate(getCompound1());
        IAtomContainer ac2 = duplicate(getCompound2());

        // Guard: cannot compute MCS on empty molecules
        if (ac1 == null || ac2 == null || ac1.getAtomCount() == 0 || ac2.getAtomCount() == 0) {
            return null;
        }
        Isomorphism isomorphism;
        int expectedMaxGraphmatch = expectedMaxGraphmatch(ac1, ac2);
        //boolean moleculeConnected = isMoleculeConnected(ac1, ac2);
        boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

        LOGGER.debug("Expected matches " + expectedMaxGraphmatch);

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
            LOGGER.debug("===={Aladdin} Mapping {Gini}====");
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

    private IAtomContainer duplicate(IAtomContainer ac) throws CloneNotSupportedException {
        IAtomContainer a = ac.clone();
        a.setID(ac.getID());
        a.setProperties(ac.getProperties());

        for (int i = 0; i < a.getAtomCount(); i++) {
            a.getAtom(i).setID(ac.getAtom(i).getID());
        }

        // Fix aromatic bond consistency: if a bond is aromatic but its atoms
        // are not flagged aromatic, downgrade the bond to SINGLE to prevent
        // "Aromatic bond connects non-aromatic atoms" errors in SMSD
        for (IBond bond : a.bonds()) {
            if (bond.isAromatic()) {
                IAtom begin = bond.getBegin();
                IAtom end = bond.getEnd();
                if ((begin != null && !begin.isAromatic())
                        || (end != null && !end.isAromatic())) {
                    bond.setIsAromatic(false);
                    if (bond.getOrder() == null || bond.getOrder() == IBond.Order.UNSET) {
                        bond.setOrder(IBond.Order.SINGLE);
                    }
                }
            }
        }

        return a;
    }

    /**
     * @return the compound1
     */
    IAtomContainer getCompound1() {
        return compound1;
    }

    /**
     * @return the compound2
     */
    IAtomContainer getCompound2() {
        return compound2;
    }

    /**
     * @return the queryPosition
     */
    int getQueryPosition() {
        return queryPosition;
    }

    /**
     * @return the targetPosition
     */
    int getTargetPosition() {
        return targetPosition;
    }

    void setHasPerfectRings(boolean ring) {
        this.hasRings = ring;
    }

    /**
     * @return the hasRings
     */
    boolean isHasPerfectRings() {
        return hasRings;
    }

    /*
     * Check if fragmented container has single atom
     */
    private boolean isMoleculeConnected(IAtomContainer compound1, IAtomContainer compound2) {
        LOGGER.debug("isMoleculeConnected");
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

        return connected1 && connected2;
    }

    void setEductRingCount(int numberOfCyclesEduct) {
        this.numberOfCyclesEduct = numberOfCyclesEduct;
    }

    void setProductRingCount(int numberOfCyclesProduct) {
        this.numberOfCyclesProduct = numberOfCyclesProduct;
    }

    private static final String CACHED_SMILES = "CACHED_CANONICAL_SMILES";
    private static final SmilesGenerator CANONICAL_SMIGEN = new SmilesGenerator(
            SmiFlavor.Canonical | SmiFlavor.Stereo);

    /**
     * Generate a unique cache key based on canonical SMILES (structure-based,
     * not ID-based). This enables cross-reaction cache hits: if the same
     * molecule pair appears in different reactions, the MCS result is reused.
     */
    String generateUniqueKey(
            String id1, String id2,
            int atomCount1, int atomCount2,
            int bondCount1, int bondCount2,
            boolean atomtypeMatcher,
            boolean bondMatcher,
            boolean ringMatcher,
            boolean hasPerfectRings,
            int numberOfCyclesEduct, int numberOfCyclesProduct) {

        StringBuilder key = new StringBuilder();

        // Use canonical SMILES as the molecular identity (not mol IDs)
        String smi1 = getCanonicalSmiles(compound1);
        String smi2 = getCanonicalSmiles(compound2);
        key.append(smi1).append(">>").append(smi2);

        // Append matcher flags that affect the MCS result
        key.append('|')
                .append(atomtypeMatcher ? '1' : '0')
                .append(bondMatcher ? '1' : '0')
                .append(ringMatcher ? '1' : '0')
                .append(hasPerfectRings ? '1' : '0');

        return key.toString();
    }

    /**
     * Get or compute canonical SMILES for a molecule. Cached on the molecule
     * to avoid recomputation across calls within the same reaction.
     */
    private String getCanonicalSmiles(IAtomContainer mol) {
        String cached = mol.getProperty(CACHED_SMILES);
        if (cached != null) {
            return cached;
        }
        try {
            cached = CANONICAL_SMIGEN.create(mol);
        } catch (CDKException e) {
            // Fallback: use atom/bond counts + fingerprint hash
            cached = mol.getAtomCount() + ":" + mol.getBondCount() + ":" + mol.hashCode();
        }
        mol.setProperty(CACHED_SMILES, cached);
        return cached;
    }

    /*
     * copy old mapping from the cache to new
     */
    MCSSolution copyOldSolutionToNew(int queryPosition, int targetPosition,
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

    MCSSolution addMCSSolution(String key, ThreadSafeCache<String, MCSSolution> mappingcache, Isomorphism isomorphism) {

        isomorphism.setChemFilters(true, true, true);
        try {
            LOGGER.debug("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", "
                    + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
        } catch (CloneNotSupportedException | CDKException e) {
            LOGGER.error(SEVERE, "Error in computing MCS ", e.getMessage());
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
        long stopTime = currentTimeMillis();
        long time = stopTime - startTime;
        printMatch(isomorphism);
        LOGGER.debug("\" Time:\" " + time);
        if (!mappingcache.containsKey(key)) {
            LOGGER.debug("Key " + key);
            try {
                LOGGER.debug("mcs size " + mcs.getAtomAtomMapping().getCount());
                LOGGER.debug("mcs map " + mcs.getAtomAtomMapping().getMappingsByIndex());
                LOGGER.debug("mcs " + mcs.getAtomAtomMapping().getCommonFragmentAsSMILES());
                LOGGER.debug("\n\n\n ");
            } catch (CloneNotSupportedException | CDKException ex) {
                LOGGER.error(SEVERE, "Unable to create SMILES ", ex.getMessage());
            }
            mappingcache.put(key, mcs);
        }
        return mcs;
    }
}
