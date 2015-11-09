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
package uk.ac.ebi.reactionblast.mapping.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import uk.ac.ebi.reactionblast.tools.labelling.SmilesMoleculeLabeller;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public final class MCSThread implements Callable<MCSSolution> {

    private boolean stereoFlag;
    private boolean fragmentFlag;
    private boolean energyFlag;

    private final static boolean DEBUG1 = false;
    private final static boolean DEBUG2 = false;

    private SmilesGenerator smiles;
    private Aromaticity aromaticity;

    protected final IAtomContainer compound1;
    protected final IAtomContainer compound2;
    protected final int queryPosition;
    protected final int targetPosition;
    protected final boolean bondMatcher;
    protected final boolean ringMatcher;
    protected final IMappingAlgorithm theory;
    protected final boolean atomMatcher;
    private final ICanonicalMoleculeLabeller labeller;

    final long startTime;
    private boolean hasRings;
    private Integer eductCount;
    private Integer productCount;

    /**
     *
     * @param theory
     * @param educt
     * @param product
     * @param bondMatcher
     * @param ringMatcher
     * @throws CloneNotSupportedException
     */
    MCSThread(IMappingAlgorithm theory, int queryPosition, int targetPosition,
            IAtomContainer educt, IAtomContainer product,
            boolean bondMatcher, boolean ringMatcher, boolean atomMatcher)
            throws CloneNotSupportedException, CDKException {
        this.stereoFlag = true;
        this.fragmentFlag = true;
        this.energyFlag = true;

        this.startTime = System.currentTimeMillis();

        this.compound1 = getNewContainerWithIDs(educt);
        this.compound2 = getNewContainerWithIDs(product);
        this.queryPosition = queryPosition;
        this.targetPosition = targetPosition;
        this.bondMatcher = bondMatcher;
        this.ringMatcher = ringMatcher;
        this.theory = theory;
        this.atomMatcher = atomMatcher;
        this.labeller = new SmilesMoleculeLabeller();

        if (DEBUG1) {
            aromaticity = new Aromaticity(ElectronDonation.daylight(),
                    Cycles.or(Cycles.all(), Cycles.relevant()));
            smiles = SmilesGenerator.unique().aromatic();
        }

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
                    System.out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + smiles.create(compound1)
                            + " molT: " + smiles.create(compound2)
                            + " atoms: " + compound1.getAtomCount()
                            + " atoms: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + "]");
                }

                /*
                 Check if MCS matching required or not very IMP step
                 */
                boolean possibleVFmatch12 = isPossibleSubgraphMatch(getCompound1(), getCompound2());
                if (DEBUG1) {
                    System.out.println("VF Matcher " + possibleVFmatch12);
                    System.out.println("Compound1 " + "A: " + getCompound1().getAtomCount()
                            + " B: " + getCompound1().getBondCount());

                    System.out.println("Compound2 " + "A: " + getCompound2().getAtomCount()
                            + " B: " + getCompound2().getBondCount());
                }

                boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
                if (DEBUG1) {
                    System.out.println("VF Matcher " + possibleVFmatch21);
                }

                if (possibleVFmatch12
                        && getCompound1().getAtomCount() <= getCompound2().getAtomCount()
                        && getCompound1().getBondCount() <= getCompound2().getBondCount()) {
                    if (DEBUG1) {
                        System.out.println("Substructure 5");
                    }
                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());
                    Substructure substructure = new Substructure(ac1, ac2,
                            false, false, isHasPerfectRings(), true);
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
                        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                                substructure.getQuery(), substructure.getTarget(), substructure.getFirstAtomMapping());
                        mcs.setEnergy(substructure.getEnergyScore(0));
                        mcs.setFragmentSize(substructure.getFragmentSize(0));
                        mcs.setStereoScore(substructure.getStereoScore(0));
                        if (DEBUG1) {
                            long stopTime = System.currentTimeMillis();
                            long time = stopTime - startTime;
                            printMatch(substructure);
                            System.out.println("\" Time:\" " + time);
                        }
                        return mcs;
                    } else {
                        if (DEBUG1) {
                            System.out.println("not a Substructure 5");
                        }
                    }
                } else if (possibleVFmatch21) {

                    if (DEBUG1) {
                        System.out.println("Substructure 6");
                    }

                    IAtomContainer ac1 = duplicate(getCompound1());
                    IAtomContainer ac2 = duplicate(getCompound2());
                    BaseMapping substructure = new Substructure(ac2, ac1,
                            false, false, isHasPerfectRings(), true);
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);

                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {
                        AtomAtomMapping aam = new AtomAtomMapping(substructure.getTarget(), substructure.getQuery());
                        Map<IAtom, IAtom> mappings = substructure.getFirstAtomMapping().getMappingsByAtoms();
                        for (IAtom atom1 : mappings.keySet()) {
                            IAtom atom2 = mappings.get(atom1);
                            aam.put(atom2, atom1);
                        }
                        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                                substructure.getTarget(), substructure.getQuery(), aam);
                        mcs.setEnergy(substructure.getEnergyScore(0));
                        mcs.setFragmentSize(substructure.getFragmentSize(0));
                        mcs.setStereoScore(substructure.getStereoScore(0));

                        if (DEBUG1) {
                            long stopTime = System.currentTimeMillis();
                            long time = stopTime - startTime;
                            printMatch(substructure);
                            System.out.println("\" Time:\" " + time);
                        }
                        return mcs;
                    } else {
                        if (DEBUG1) {
                            System.out.println("not a Substructure 6");
                        }
                    }
                }

                if (DEBUG1) {
                    System.out.println("calling mcs");
                    System.out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + smiles.create(compound1)
                            + " molT: " + smiles.create(compound2)
                            + " atoms: " + compound1.getAtomCount()
                            + " atoms: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + "]");
                }
                MCSSolution mcs = mcs();
                return mcs;
            } else {
                MCSSolution mcs = mcs();
                return mcs;
            }

        } catch (CDKException | CloneNotSupportedException ex) {
            Logger.getLogger(MCSThread.class.getName()).log(Level.SEVERE, null, ex);
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
                    ? String.valueOf(i) : mol.getAtom(i).getID();
            ac.getAtom(i).setID(atomID);
        }
        String containerID = mol.getID() == null ? String.valueOf(System.nanoTime()) : mol.getID();
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
            for (String k : atomUniqueCounter1.keySet()) {
                if (atomUniqueCounter1.get(k) > atomUniqueCounter2.get(k)) {
                    return false;
                }
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
            String hyb = a.getHybridization() == CDKConstants.UNSET
                    ? a.getSymbol() : a.getAtomTypeName();
            atomUniqueCounter1.add(hyb);
        }

        for (IAtom b : t.atoms()) {
            String hyb = b.getHybridization() == CDKConstants.UNSET
                    ? b.getSymbol() : b.getAtomTypeName();
            atomUniqueCounter2.add(hyb);
        }

        Collections.sort(atomUniqueCounter1);
        Collections.sort(atomUniqueCounter2);

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

    synchronized MCSSolution mcs() {
        try {
            /*
             * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
             */
            Isomorphism isomorphism;
            int expectedMaxGraphmatch = expectedMaxGraphmatch(getCompound1(), getCompound2());

            if (eductCount == 1 && productCount == 1) {
                /*
                 This handles large aliphatics to ring system (ex: R09907)
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT,
                        false, isHasPerfectRings(), false);
            } else if (expectedMaxGraphmatch > 30) {
                /*
                 This handles large aliphatics to ring system (ex: R06466)
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.VFLibMCS,
                        false, ringMatcher, !isHasPerfectRings());
            } else {
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT,
                        false, isHasPerfectRings(), !isHasPerfectRings());
            }

            isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
            if (DEBUG1) {
                System.out.println("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", " + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
            }      /*
             * In case of Complete subgraph, don't use Energy filter
             *
             */

            MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                    isomorphism.getQuery(), isomorphism.getTarget(), isomorphism.getFirstAtomMapping());
            mcs.setEnergy(isomorphism.getEnergyScore(0));
            mcs.setFragmentSize(isomorphism.getFragmentSize(0));
            mcs.setStereoScore(isomorphism.getStereoScore(0));

            if (DEBUG1) {
                long stopTime = System.currentTimeMillis();
                long time = stopTime - startTime;
                printMatch(isomorphism);
                System.out.println("\" Time:\" " + time);

            }
            return mcs;
        } catch (Exception e) {
            Logger.getLogger(MCSThread.class.getName()).log(Level.SEVERE, "Error in computing MCS ", e);
        }
        return null;
    }

    synchronized MCSSolution combimcs(boolean stereoFlag, boolean fragmentFlag,
            boolean energyFlag) throws CloneNotSupportedException, CDKException {
        double energy = 0.0d;
        int fragmentSize = 0;
        int stereoScore = 0;
        IAtomContainer ac1 = duplicate(getCompound1());
        IAtomContainer ac2 = duplicate(getCompound2());
        /*
         * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
         */
        Isomorphism isomorphism;

        isomorphism = new Isomorphism(ac1, ac2, Algorithm.DEFAULT, true, ringMatcher, true);
        isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);

        Map<IAtom, IAtom> acceptedSolution = new HashMap<>();

        if (!isomorphism.getFirstAtomMapping().isEmpty()) {
            for (IAtom a : isomorphism.getFirstAtomMapping().getMappingsByAtoms().keySet()) {
                IAtom refA = getAtomByID(getCompound1(), a);
                IAtom refB = getAtomByID(getCompound2(),
                        isomorphism.getFirstAtomMapping().getMappingsByAtoms().get(a));
                acceptedSolution.put(refA, refB);
            }
            ac1 = reduceAtomContainer(ac1, isomorphism.getFirstAtomMapping().getMappingsByAtoms().keySet());
            ac2 = reduceAtomContainer(ac2, isomorphism.getFirstAtomMapping().getMappingsByAtoms().values());
            energy += isomorphism.getEnergyScore(0);
            fragmentSize += isomorphism.getFragmentSize(0);
            stereoScore += isomorphism.getStereoScore(0);
        }

        if (ac1.getAtomCount() > 0 && ac2.getAtomCount() > 0) {

            if (DEBUG2) {
                System.out.println(smiles.create(getCompound1()) + "ac1 reduced by "
                        + (getCompound1().getAtomCount() - ac1.getAtomCount())
                        + ", " + smiles.create(getCompound2()) + " ac2 reduced by "
                        + (getCompound2().getAtomCount() - ac2.getAtomCount()));
            }
            isomorphism = new Isomorphism(ac1, ac2, Algorithm.VFLibMCS, false, ringMatcher, true);
            isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
            List<AtomAtomMapping> allAtomMapping = isomorphism.getAllAtomMapping();
            int solIndex = 0;
            for (AtomAtomMapping s : allAtomMapping) {
                boolean stitchingFeasible = isStitchingFeasible(getCompound1(),
                        getCompound2(), acceptedSolution, s);
                if (stitchingFeasible) {
                    for (IAtom a : s.getMappingsByAtoms().keySet()) {
                        IAtom refA = getAtomByID(getCompound1(), a);
                        IAtom refB = getAtomByID(getCompound2(), s.getMappingsByAtoms().get(a));
                        acceptedSolution.put(refA, refB);
                    }
                    energy += isomorphism.getEnergyScore(solIndex);
                    fragmentSize += isomorphism.getFragmentSize(solIndex);
                    stereoScore += isomorphism.getStereoScore(solIndex);
                    ac1 = reduceAtomContainer(ac1, s.getMappingsByAtoms().keySet());
                    ac2 = reduceAtomContainer(ac2, s.getMappingsByAtoms().values());
                    break;
                }
                solIndex++;
            }

            /*
             * In case of Complete subgraph, don't use Energy filter
             *
             */
            AtomAtomMapping combi = new AtomAtomMapping(getCompound1(), getCompound2());

            for (IAtom a : acceptedSolution.keySet()) {
                IAtom b = acceptedSolution.get(a);
                combi.put(a, b);
            }
            MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(), isomorphism.getQuery(), isomorphism.getTarget(), combi);
            mcs.setEnergy(energy);
            mcs.setFragmentSize(fragmentSize);
            mcs.setStereoScore(stereoScore);

        }

        /*
         * In case of Complete subgraph, don't use Energy filter
         *
         */
        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(), isomorphism.getQuery(), isomorphism.getTarget(), isomorphism.getFirstAtomMapping());
        mcs.setEnergy(isomorphism.getEnergyScore(0));
        mcs.setFragmentSize(isomorphism.getFragmentSize(0));
        mcs.setStereoScore(isomorphism.getStereoScore(0));

        if (DEBUG1) {
            long stopTime = System.currentTimeMillis();
            long time = stopTime - startTime;
            printMatch(isomorphism);
            System.out.println("\" Time:\" " + time);

        }

        return mcs;
    }

    private IAtomContainer reduceAtomContainer(IAtomContainer ac,
            Collection<IAtom> keySet) throws CloneNotSupportedException {
        IAtomContainer ac_new = ac.getBuilder().newInstance(IAtomContainer.class);
        HashMap<IAtom, IAtom> ref_new_atom = new HashMap<>();
        for (int i = 0; i < ac.getAtomCount(); i++) {
            IAtom ref = ac.getAtom(i);
            if (!keySet.contains(ref)) {
                IAtom a = ref.clone();
                a.setID(ref.getID());
                ac_new.addAtom(a);
                ref_new_atom.put(ref, a);
            }
        }

        for (IBond b : ac.bonds()) {
            if (keySet.contains(b.getAtom(0)) || keySet.contains(b.getAtom(1))) {
                continue;
            }
            IAtom a1 = ref_new_atom.get(b.getAtom(0));
            IAtom a2 = ref_new_atom.get(b.getAtom(1));
            IBond.Order order = b.getOrder();
            ac_new.addBond(ac_new.getAtomNumber(a1), ac_new.getAtomNumber(a2), order);
        }

        ac_new.setID(ac.getID());
        return ac_new;
    }

    private synchronized IAtom getAtomByID(IAtomContainer ac, IAtom atom) {
        if (atom.getID() == null) {
            return null;
        }
        for (IAtom a : ac.atoms()) {
            if (a.getID().equals(atom.getID())) {
                return a;
            }
        }
        return null;
    }

    private boolean isStitchingFeasible(IAtomContainer compound1, IAtomContainer compound2,
            Map<IAtom, IAtom> map, AtomAtomMapping mapping) {

        boolean t1 = false;
        boolean t2 = false;
        for (IAtom a : map.keySet()) {
            IAtom refAtomA = getAtomByID(compound1, a);
            for (IAtom b : mapping.getMappingsByAtoms().keySet()) {
                IAtom refAtomB = getAtomByID(compound1, b);
                IBond bond = compound1.getBond(refAtomA, refAtomB);
                if (bond != null) {
                    t1 = true;
                    break;
                }
            }
        }

        for (IAtom a : map.values()) {
            IAtom refAtomA = getAtomByID(compound2, a);
            for (IAtom b : mapping.getMappingsByAtoms().values()) {
                IAtom refAtomB = getAtomByID(compound2, b);
                IBond bond = compound2.getBond(refAtomA, refAtomB);
                if (bond != null) {
                    t2 = true;
                    break;
                }
            }
        }
        return t1 && t2;
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

    void setEductCount(Integer eductCount) {
        this.eductCount = eductCount;
    }

    void setProductCount(Integer productCount) {
        this.productCount = productCount;
    }
}
