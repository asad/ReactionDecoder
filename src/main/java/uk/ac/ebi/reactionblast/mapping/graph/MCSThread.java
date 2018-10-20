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
package uk.ac.ebi.reactionblast.mapping.graph;

import static java.lang.String.valueOf;
import static java.lang.System.currentTimeMillis;
import static java.lang.System.nanoTime;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.exception.CDKException;
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
import org.openscience.smsd.interfaces.Algorithm;
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

    private boolean stereoFlag;
    private boolean fragmentFlag;
    private boolean energyFlag;

    private SmilesGenerator smiles;

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
    protected final boolean bondMatcher;

    /**
     *
     */
    protected final boolean ringMatcher;

    /**
     *
     */
    protected final IMappingAlgorithm theory;

    /**
     *
     */
    protected final boolean atomMatcher;

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
            IAtomContainer educt, IAtomContainer product,
            boolean bondMatcher, boolean ringMatcher, boolean atomMatcher)
            throws CloneNotSupportedException, CDKException {
        this.stereoFlag = true;
        this.fragmentFlag = true;
        this.energyFlag = true;
        this.compound1 = getNewContainerWithIDs(educt);
        this.compound2 = getNewContainerWithIDs(product);
        this.queryPosition = queryPosition;
        this.targetPosition = targetPosition;
        this.bondMatcher = bondMatcher;
        this.ringMatcher = ringMatcher;
        this.theory = theory;
        this.atomMatcher = atomMatcher;
        this.numberOfCyclesEduct = 0;
        this.numberOfCyclesProduct = 0;

        Aromaticity aromaticity = new Aromaticity(daylight(),
                Cycles.or(Cycles.all(),
                        Cycles.or(Cycles.relevant(),
                                Cycles.essential())));
        aromaticity.apply(this.compound1);
        aromaticity.apply(this.compound2);

        /*
         * create SMILES
         */
        smiles = new SmilesGenerator(
                //                SmiFlavor.Unique |
                SmiFlavor.Stereo
                | SmiFlavor.AtomAtomMap);

    }

    synchronized void printMatch(BaseMapping isomorphism) {
        int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                : isomorphism.getFirstAtomMapping().getCount();

        System.out.println("Q: " + isomorphism.getQuery().getID()
                + " T: " + isomorphism.getTarget().getID()
                + " atoms: " + isomorphism.getQuery().getAtomCount()
                + " atoms: " + isomorphism.getTarget().getAtomCount()
                + " overlaps " + overlap);
    }

    @Override
    public synchronized MCSSolution call() throws Exception {
        try {
            if (!theory.equals(IMappingAlgorithm.RINGS)) {
                if (DEBUG1) {
                    String createSM1 = null;
                    String createSM2 = null;
                    try {
                        createSM1 = smiles.create(this.compound1);
                        createSM2 = smiles.create(this.compound2);
                    } catch (CDKException e) {
                        LOGGER.error(SEVERE, null, e);
                    }
                    System.out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + createSM1
                            + " molT: " + createSM2
                            + " atomsQ: " + compound1.getAtomCount()
                            + " atomsT: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + " isHasPerfectRings: " + isHasPerfectRings()
                            + "]");
                }

                /*
                 * IMP: Do not perform substructure matching for disconnected molecules
                 */
                boolean moleculeConnected = isMoleculeConnected(getCompound1(), getCompound2());

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

                if (moleculeConnected && possibleVFmatch12
                        && getCompound1().getAtomCount() <= getCompound2().getAtomCount()
                        && getCompound1().getBondCount() <= getCompound2().getBondCount()) {
                    if (DEBUG1) {
                        System.out.println("Substructure 1");
                        this.startTime = currentTimeMillis();

                    }
                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());
                    Substructure substructure;
                    substructure = new Substructure(ac1, ac2,
                            true, false, isHasPerfectRings(), true);
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac1, ac2,
                                false, false, isHasPerfectRings(), true);
                    }
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac1, ac2,
                                false, false, isHasPerfectRings(), false);
                    }
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
//                    System.out.println("Number of Solutions: " + substructure.getAllAtomMapping());
                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
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
                            System.out.println("\" Time:\" " + time);
                        }
                        return mcs;
                    } else if (DEBUG1) {
                        System.out.println("not a Substructure 1");
                    }
                } else if (moleculeConnected && possibleVFmatch21) {

                    if (DEBUG1) {
                        System.out.println("Substructure 2");
                        this.startTime = currentTimeMillis();

                    }

                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());
                    Substructure substructure;
                    substructure = new Substructure(ac2, ac1,
                            true, false, isHasPerfectRings(), true);
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac2, ac1,
                                false, false, isHasPerfectRings(), true);
                    }
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac2, ac1,
                                false, false, isHasPerfectRings(), false);
                    }
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);

                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {

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
                if (DEBUG1) {

                    /*
                     * create SMILES
                     */
                    smiles = new SmilesGenerator(SmiFlavor.Unique
                            | SmiFlavor.Stereo
                            | SmiFlavor.AtomAtomMap);
                    String createSM1 = null;
                    String createSM2 = null;
                    try {
                        createSM1 = smiles.create(this.compound1);
                        createSM2 = smiles.create(this.compound2);
                    } catch (CDKException e) {
                        LOGGER.error(SEVERE, null, e);
                    }
                    System.out.println("No Substructure found - switching to MCS");
                    System.out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + createSM1
                            + " molT: " + createSM2
                            + " atomsE: " + compound1.getAtomCount()
                            + " atomsP: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + " isHasPerfectRings: " + isHasPerfectRings()
                            + "]");
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
            } else {
                if (DEBUG1) {

                    /*
                     * create SMILES
                     */
                    smiles = new SmilesGenerator(SmiFlavor.Unique
                            | SmiFlavor.Stereo
                            | SmiFlavor.AtomAtomMap);
                    String createSM1 = null;
                    String createSM2 = null;
                    try {
                        createSM1 = smiles.create(this.compound1);
                        createSM2 = smiles.create(this.compound2);
                    } catch (CDKException e) {
                        LOGGER.error(SEVERE, null, e);
                    }
                    System.out.println("No Substructure found - switching to MCS");
                    System.out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + createSM1
                            + " molT: " + createSM2
                            + " atomsE: " + compound1.getAtomCount()
                            + " atomsP: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + " isHasPerfectRings: " + isHasPerfectRings()
                            + "]");
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
            }

        } catch (CDKException | CloneNotSupportedException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return null;
    }

    private IAtomContainer getNewContainerWithIDs(IAtomContainer mol)
            throws CDKException, CloneNotSupportedException {
        /*
         Generating SMILES speed ups the mapping process by n-FOLDS
         May be this is CDK inherent bug, which relies on the properties set by the SMILES
         */
        IAtomContainer ac;
        ac = mol.clone();
        for (int i = 0; i < ac.getAtomCount(); i++) {
            String atomID = mol.getAtom(i).getID() == null
                    ? valueOf(i) : mol.getAtom(i).getID();
            ac.getAtom(i).setID(atomID);
        }
        String containerID = mol.getID() == null ? valueOf(nanoTime()) : mol.getID();
        ac.setID(containerID);

        return ac;
    }

    synchronized void setStereoFlag(boolean b) {
        this.stereoFlag = b;
    }

    /**
     * @param fragmentFilterFlag the fragmentFlag to set
     */
    synchronized void setFragmentFilterFlag(boolean fragmentFilterFlag) {
        this.fragmentFlag = fragmentFilterFlag;
    }

    /**
     * @param energyFlag the energyFlag to set
     */
    synchronized void setEnergyFlag(boolean energyFlag) {
        this.energyFlag = energyFlag;
    }

    private boolean isPossibleSubgraphMatch(IAtomContainer q, IAtomContainer t) {

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
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
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

    private int expectedMaxGraphmatch(IAtomContainer q, IAtomContainer t) {

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
            System.out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            System.out.println("common " + common.size());
        }
        atomUniqueCounter1.clear();
        atomUniqueCounter2.clear();
        return common.size();
    }

    synchronized MCSSolution mcs() throws CDKException, CloneNotSupportedException {
        /*
         * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
         */
        IAtomContainer ac1 = duplicate(getCompound1());
        IAtomContainer ac2 = duplicate(getCompound2());
        Isomorphism isomorphism;
        int expectedMaxGraphmatch = expectedMaxGraphmatch(ac1, ac2);
        boolean moleculeConnected = isMoleculeConnected(ac1, ac2);
        boolean ringFlag = this.numberOfCyclesEduct > 0 && this.numberOfCyclesProduct > 0;

        if (DEBUG3) {
            System.out.println("Expected matches " + expectedMaxGraphmatch);
        }

        isomorphism = new Isomorphism(ac1, ac2, Algorithm.DEFAULT,
                false, isHasPerfectRings(), false);

        isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
        if (DEBUG3) {
            try {
                System.out.println("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", "
                        + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
            } catch (CloneNotSupportedException | CDKException e) {
                LOGGER.error(SEVERE, "Error in computing MCS ", e);
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
        return mcs;
    }

    private IAtomContainer duplicate(IAtomContainer ac) throws CloneNotSupportedException {
        IAtomContainer a = ac.clone();
        a.setID(ac.getID());
        for (int i = 0; i < a.getAtomCount(); i++) {
            a.getAtom(i).setID(ac.getAtom(i).getID());
        }
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
        if (DEBUG1) {
            System.out.println("isMoleculeConnected");
        }
        boolean connected1 = true;

        IAtomContainerSet partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound1);
        for (IAtomContainer a : partitionIntoMolecules.atomContainers()) {
            if (DEBUG1) {
                System.out.println("QContainer size " + a.getAtomCount());
            }
            if (a.getAtomCount() == 1) {
                connected1 = false;
            }
        }

        boolean connected2 = true;

        partitionIntoMolecules = ConnectivityChecker.partitionIntoMolecules(compound2);
        for (IAtomContainer a : partitionIntoMolecules.atomContainers()) {
            if (DEBUG1) {
                System.out.println("TContainer size " + a.getAtomCount());
            }
            if (a.getAtomCount() == 1) {
                connected2 = false;
            }
        }

        if (DEBUG1) {
            System.out.println("Flag 1 " + connected1);
        }
        if (DEBUG1) {
            System.out.println("Flag 2 " + connected2);
        }
        return connected1 & connected2;
    }

    void setEductRingCount(int numberOfCyclesEduct) {
        this.numberOfCyclesEduct = numberOfCyclesEduct;
    }

    void setProductRingCount(int numberOfCyclesProduct) {
        this.numberOfCyclesProduct = numberOfCyclesProduct;
    }
}
