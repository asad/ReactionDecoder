/*
 * Copyright (C) 2003-2017 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.sort;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import static org.openscience.cdk.graph.Cycles.all;
import static org.openscience.cdk.graph.Cycles.or;
import static org.openscience.cdk.graph.Cycles.relevant;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;
import static org.openscience.smsd.interfaces.Algorithm.DEFAULT;
import static org.openscience.smsd.interfaces.Algorithm.VFLibMCS;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import uk.ac.ebi.reactionblast.tools.labelling.SmilesMoleculeLabeller;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MCSThread implements Callable<MCSSolution> {

    private static final boolean DEBUG1 = false;
    private static final boolean DEBUG2 = false;
    private static final boolean DEBUG3 = false;
    private static final Logger LOG = getLogger(MCSThread.class.getName());

    private boolean stereoFlag;
    private boolean fragmentFlag;
    private boolean energyFlag;

    private SmilesGenerator smiles;
    private Aromaticity aromaticity;

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
    private final ICanonicalMoleculeLabeller labeller;

    final long startTime;
    private boolean hasRings;
    private Integer eductCount;
    private Integer productCount;

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

        this.startTime = currentTimeMillis();

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
            aromaticity = new Aromaticity(daylight(),
                    or(all(), relevant()));
            smiles = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        }

    }

    synchronized void printMatch(BaseMapping isomorphism) {
        int overlap = isomorphism.getFirstAtomMapping().isEmpty() ? 0
                : isomorphism.getFirstAtomMapping().getCount();

        out.println("Q: " + isomorphism.getQuery().getID()
                + " T: " + isomorphism.getTarget().getID()
                + " atoms: " + isomorphism.getQuery().getAtomCount()
                + " atoms: " + isomorphism.getTarget().getAtomCount()
                + " overlaps " + overlap);
    }

    @Override
    public synchronized MCSSolution call() throws Exception {
        try {
            if (!theory.equals(RINGS)) {
                if (DEBUG1) {
                    out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + smiles.create(compound1)
                            + " molT: " + smiles.create(compound2)
                            + " atomsE: " + compound1.getAtomCount()
                            + " atomsP: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + " isHasPerfectRings: " + isHasPerfectRings()
                            + "]");
                }

                /*
                 Check if MCS matching required or not very IMP step
                 */
                boolean possibleVFmatch12 = isPossibleSubgraphMatch(getCompound1(), getCompound2());
                if (DEBUG1) {
                    out.println("VF Matcher 1->2 " + possibleVFmatch12);
                }

                boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
                if (DEBUG1) {
                    out.println("VF Matcher 2->1 " + possibleVFmatch21);
                }

                if (possibleVFmatch12
                        && getCompound1().getAtomCount() <= getCompound2().getAtomCount()
                        && getCompound1().getBondCount() <= getCompound2().getBondCount()) {
                    if (DEBUG1) {
                        out.println("Substructure 5");
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
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
//                    System.out.println("Number of Solutions: " + substructure.getAllAtomMapping());
                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac1.getAtomCount()) {
                        MCSSolution mcs = new MCSSolution(getQueryPosition(), getTargetPosition(),
                                substructure.getQuery(), substructure.getTarget(), substructure.getFirstAtomMapping());
                        mcs.setEnergy(substructure.getEnergyScore(0));
                        mcs.setFragmentSize(substructure.getFragmentSize(0));
                        mcs.setStereoScore(substructure.getStereoScore(0));
                        if (DEBUG1) {
                            long stopTime = currentTimeMillis();
                            long time = stopTime - startTime;
                            printMatch(substructure);
                            out.println("\" Time:\" " + time);
                        }
                        return mcs;
                    } else if (DEBUG1) {
                        out.println("not a Substructure 5");
                    }
                } else if (possibleVFmatch21) {

                    if (DEBUG1) {
                        out.println("Substructure 6");
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
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);

                    if (substructure.isSubgraph() && substructure.getFirstAtomMapping().getCount() == ac2.getAtomCount()) {
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
                            out.println("\" Time:\" " + time);
                        }
                        return mcs;
                    } else if (DEBUG1) {
                        out.println("not a Substructure 6");
                    }
                }

                if (DEBUG1) {
                    out.println("calling mcs");
                    out.println("Q: " + getCompound1().getID()
                            + " T: " + getCompound2().getID()
                            + " molQ: " + smiles.create(compound1)
                            + " molT: " + smiles.create(compound2)
                            + " atomsQ: " + compound1.getAtomCount()
                            + " atomsT: " + compound2.getAtomCount()
                            + " [bonds: " + bondMatcher
                            + " rings: " + ringMatcher
                            + " isHasPerfectRings: " + isHasPerfectRings()
                            + "]");
                }
                MCSSolution mcs = mcs();
                return mcs;
            } else {
                MCSSolution mcs = mcs();
                return mcs;
            }

        } catch (CDKException | CloneNotSupportedException ex) {
            getLogger(MCSThread.class.getName()).log(SEVERE, null, ex);
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
            out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            out.println("diff " + difference.size());
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
            out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            out.println("atomUniqueCounter2 " + atomUniqueCounter2);
            out.println("common " + common.size());
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

            if (getCompound1().getAtomCount() == 1
                    || getCompound2().getAtomCount() == 1) {
                if (DEBUG3) {
                    System.out.println("CASE 1");
                }
                /*
                 * This handles large aliphatics to ring system (ex: R09907)
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT,
                        false, isHasPerfectRings(), false);
            } else if (expectedMaxGraphmatch >= 30
                    && ConnectivityChecker.isConnected(getCompound1())) {
                if (DEBUG3) {
                    System.out.println("CASE 2");
                }
                /*
                 * Although the bond changes are set to true but its only used by filters
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.MCSPlus,
                        false, isHasPerfectRings(), !isHasPerfectRings());
            } else if (expectedMaxGraphmatch < 30) {
                if (DEBUG3) {
                    System.out.println("CASE 3");
                }
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.CDKMCS,
                        false, isHasPerfectRings(), !isHasPerfectRings());
            } else {
                if (DEBUG3) {
                    System.out.println("CASE 4");
                }
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT,
                        false, isHasPerfectRings(), !isHasPerfectRings());
            }

            isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
            if (DEBUG3) {
                out.println("MCS " + isomorphism.getFirstAtomMapping().getCount() + ", " + isomorphism.getFirstAtomMapping().getCommonFragmentAsSMILES());
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
                out.println("\" Time:\" " + time);

            }
            return mcs;
        } catch (CloneNotSupportedException | CDKException e) {
            getLogger(MCSThread.class.getName()).log(SEVERE, "Error in computing MCS ", e);
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

        isomorphism = new Isomorphism(ac1, ac2, DEFAULT, true, ringMatcher, true);
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
                out.println(smiles.create(getCompound1()) + "ac1 reduced by "
                        + (getCompound1().getAtomCount() - ac1.getAtomCount())
                        + ", " + smiles.create(getCompound2()) + " ac2 reduced by "
                        + (getCompound2().getAtomCount() - ac2.getAtomCount()));
            }
            isomorphism = new Isomorphism(ac1, ac2, VFLibMCS, false, ringMatcher, true);
            isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
            List<AtomAtomMapping> allAtomMapping = isomorphism.getAllAtomMapping();
            int solIndex = 0;
            for (AtomAtomMapping s : allAtomMapping) {
                boolean stitchingFeasible = isStitchingFeasible(getCompound1(),
                        getCompound2(), acceptedSolution, s);
                if (stitchingFeasible) {
                    s.getMappingsByAtoms().keySet().stream().forEach((a) -> {
                        IAtom refA = getAtomByID(getCompound1(), a);
                        IAtom refB = getAtomByID(getCompound2(), s.getMappingsByAtoms().get(a));
                        acceptedSolution.put(refA, refB);
                    });
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

            acceptedSolution.keySet().stream().forEach((a) -> {
                IAtom b = acceptedSolution.get(a);
                combi.put(a, b);
            });
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
            long stopTime = currentTimeMillis();
            long time = stopTime - startTime;
            printMatch(isomorphism);
            out.println("\" Time:\" " + time);

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
