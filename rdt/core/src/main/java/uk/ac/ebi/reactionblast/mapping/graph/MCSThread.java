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

import static java.lang.System.currentTimeMillis;
import static java.lang.System.nanoTime;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.graph.Cycles.or;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import uk.ac.ebi.reactionblast.tools.labelling.SmilesMoleculeLabeller;
import org.openscience.smsd.interfaces.Algorithm;
import static java.lang.String.valueOf;
import static java.util.Collections.sort;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.graph.Cycles.all;
import static org.openscience.cdk.graph.Cycles.relevant;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MCSThread implements Callable<MCSSolution> {

    private static final boolean DEBUG1 = false;
    private static final boolean DEBUG2 = false;
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
            aromaticity = new Aromaticity(daylight(), or(all(), relevant()));
            smiles = unique().aromatic();
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
                            + " isHasPerfectRings: " + isPerfectRings()
                            + "]");
                }

                /*
                 Check if MCS matching required or not very IMP step
                 */
                boolean possibleVFmatch12 = isPossibleSubgraphMatch(getCompound1(), getCompound2());
                if (DEBUG1) {
                    out.println("VF Matcher " + possibleVFmatch12);
                }

                boolean possibleVFmatch21 = isPossibleSubgraphMatch(getCompound2(), getCompound1());
                if (DEBUG1) {
                    out.println("VF Matcher " + possibleVFmatch21);
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
                    substructure = new Substructure(ac1, ac2, true, false, isPerfectRings(), true);
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac1, ac2, false, false, isPerfectRings(), true);
                    }
                    substructure.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
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
                    substructure = new Substructure(ac2, ac1, true, false, isPerfectRings(), true);
                    if (!substructure.isSubgraph()) {
                        substructure = new Substructure(ac2, ac1, false, false, isPerfectRings(), true);
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
                            + " isHasPerfectRings: " + isPerfectRings()
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
            if (!atomUniqueCounter1.keySet().stream().noneMatch((k) -> (atomUniqueCounter1.get(k) > atomUniqueCounter2.get(k)))) {
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
//        System.out.println("MCS called");
        try {
            /*
             * 0: default Isomorphism, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
             */
            Isomorphism isomorphism;
            int expectedMaxGraphmatch = expectedMaxGraphmatch(getCompound1(), getCompound2());

            if (eductCount == 1 && productCount == 1) {
                /*
                 * This handles large aliphatics to ring system (ex: R09907)
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT, false, isPerfectRings(), false);
            } else if (expectedMaxGraphmatch > 30) {
//                System.out.println("CDK MCS called");
//                System.out.println("getCompound1() " + getCompound1().getAtomCount());
//                System.out.println("getCompound2() " + getCompound2().getAtomCount());
//                System.out.println("isPerfectRings() " + isPerfectRings());
                /*
                 * This handles large aliphatics to ring system (ex: R06466)
                 */
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.CDKMCS, false, isPerfectRings(), !isPerfectRings());
            } else {
//                System.out.println("Default called");
//                System.out.println("getCompound1() " + getCompound1().getAtomCount());
//                System.out.println("getCompound2() " + getCompound2().getAtomCount());
                isomorphism = new Isomorphism(getCompound1(), getCompound2(), Algorithm.DEFAULT, false, isPerfectRings(), !isPerfectRings());
            }

            isomorphism.setChemFilters(stereoFlag, fragmentFlag, energyFlag);
            if (DEBUG1) {
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
    boolean isPerfectRings() {
        return hasRings;
    }

    void setEductCount(Integer eductCount) {
        this.eductCount = eductCount;
    }

    void setProductCount(Integer productCount) {
        this.productCount = productCount;
    }
}
