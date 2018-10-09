/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received sourceAtom copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.rgraph;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 * This algorithm derives from the algorithm described in [Tonnelier, C. and
 * Jauffret, Ph. and Hanser, Th. and Jauffret, Ph. and Kaufmann, G., Machine
 * Learning of generic reactions: 3. An efficient algorithm for maximal common
 * substructure determination, Tetrahedron Comput. Methodol., 1990, 3:351-358]
 * and modified in the thesis of T. Hanser [Unknown BibTeXML type: HAN93].
 *
 *
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class CDKRMapHandler {

    public CDKRMapHandler() {
        this.timeout = false;
    }

    /**
     * Returns source molecule
     *
     * @return the source
     */
    public synchronized IAtomContainer getSource() {
        return source;
    }

    /**
     * Set source molecule
     *
     * @param aSource the source to set
     */
    public synchronized void setSource(IAtomContainer aSource) {
        source = aSource;
    }

    /**
     * Returns target molecule
     *
     * @return the target
     */
    public synchronized IAtomContainer getTarget() {
        return target;
    }

    /**
     * Set target molecule
     *
     * @param aTarget the target to set
     */
    public synchronized void setTarget(IAtomContainer aTarget) {
        target = aTarget;
    }
    private List<Map<Integer, Integer>> mappings;
    private IAtomContainer source;
    private IAtomContainer target;
    private boolean timeout;

    /**
     * This function calculates all the possible combinations of MCS
     *
     * @param molecule1
     * @param molecule2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return List
     * @throws CDKException
     */
    public synchronized List<Map<Integer, Integer>> calculateOverlapsAndReduce(IAtomContainer molecule1,
            IAtomContainer molecule2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        setSource(molecule1);
        setTarget(molecule2);
        List<Map<Integer, Integer>> solution = new ArrayList<>();
        setMappings(solution);

        if ((getSource().getAtomCount() == 1) || (getTarget().getAtomCount() == 1)) {
            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(getSource(), getTarget());
            this.setTimeout(CDKMCS.isTimeout());
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifySingleAtomsMatchedParts(overlaps, getSource(), getTarget());

            }

        } else {
            List<List<CDKRMap>> overlaps = CDKMCS.search(getSource(), getTarget(), new BitSet(), new BitSet(), true, true, shouldMatchBonds, shouldMatchRings, matchAtomType);
            this.setTimeout(CDKMCS.isTimeout());
            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);
            while (!allMaxOverlaps.empty()) {
//                System.out.println("source: " + source.getAtomCount() + ", target: " + target.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                List<List<CDKRMap>> maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), getSource(), getTarget());
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, getSource(), getTarget());
//                identifyMatchedParts(allMaxOverlaps.peek(), source, target);
                allMaxOverlaps.pop();
            }
        }
        return solution;
    }

    /**
     * This function calculates all the possible combinations of MCS
     *
     * @param molecule1
     * @param molecule2
     * @return List
     * @throws CDKException
     */
    public synchronized List<Map<Integer, Integer>> calculateOverlapsAndReduce(IAtomContainer molecule1,
            IQueryAtomContainer molecule2) throws CDKException {
        setSource(molecule1);
        setTarget(molecule2);
        List<Map<Integer, Integer>> solution = new ArrayList<>();
        setMappings(solution);

        if ((getSource().getAtomCount() == 1) || (getTarget().getAtomCount() == 1)) {
            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(getSource(), getTarget());
            this.setTimeout(CDKMCS.isTimeout());
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifySingleAtomsMatchedParts(overlaps, getSource(), (IQueryAtomContainer) getTarget());

            }

        } else {
            List<List<CDKRMap>> overlaps = CDKMCS.search(getSource(), (IQueryAtomContainer) getTarget(), new BitSet(), new BitSet(), true, true, true, true, true);
            this.setTimeout(CDKMCS.isTimeout());
            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);
            while (!allMaxOverlaps.empty()) {
//                System.out.println("source: " + source.getAtomCount() + ", target: " + target.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                List<List<CDKRMap>> maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), getSource(), (IQueryAtomContainer) getTarget());
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, getSource(), (IQueryAtomContainer) getTarget());
//                identifyMatchedParts(allMaxOverlaps.peek(), source, target);
                allMaxOverlaps.pop();
            }
        }
        return solution;
    }

    /**
     * This function calculates only one solution (exact) because we are looking
     * at the molecules which are exactly same in terms of the bonds and atoms
     * determined by the Fingerprint
     *
     * @param Molecule1
     * @param Molecule2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @throws CDKException
     */
    public synchronized void calculateOverlapsAndReduceExactMatch(
            IAtomContainer Molecule1,
            IAtomContainer Molecule2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws CDKException {

        setSource(Molecule1);
        setTarget(Molecule2);

        setMappings(new ArrayList<>());

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);
        if ((getSource().getAtomCount() == 1) || (getTarget().getAtomCount() == 1)) {

            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(getSource(), getTarget());
            this.setTimeout(CDKMCS.isTimeout());
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                identifySingleAtomsMatchedParts(overlaps, getSource(), getTarget());
            }

        } else {

            List<List<CDKRMap>> overlaps
                    = CDKMCS.search(getSource(), getTarget(), new BitSet(), new BitSet(), true, true,
                            shouldMatchBonds, shouldMatchRings, matchAtomType);
            this.setTimeout(CDKMCS.isTimeout());
            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);

            while (!allMaxOverlaps.empty()) {
                List<List<CDKRMap>> maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), getSource(), getTarget());
                identifyMatchedParts(maxOverlapsAtoms, getSource(), getTarget());
                allMaxOverlaps.pop();
            }
        }
    }

    /**
     * This function calculates only one solution (exact) because we are looking
     * at the molecules which are exactly same in terms of the bonds and atoms
     * determined by the Fingerprint
     *
     * @param Molecule1
     * @param Molecule2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return List
     * @throws CDKException
     */
    public synchronized List<Map<Integer, Integer>> calculateSubGraphs(IAtomContainer Molecule1,
            IAtomContainer Molecule2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws CDKException {

        setSource(Molecule1);
        setTarget(Molecule2);

        List<Map<Integer, Integer>> solutions = new ArrayList<>();
        setMappings(solutions);

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);
        if ((getSource().getAtomCount() == 1) || (getTarget().getAtomCount() == 1)) {

            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(getSource(), getTarget());
            this.setTimeout(CDKMCS.isTimeout());
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                identifySingleAtomsMatchedParts(overlaps, getSource(), getTarget());
            }

        } else {

            List<List<CDKRMap>> overlaps
                    = CDKMCS.getSubgraphMaps(getSource(), getTarget(),
                            shouldMatchBonds,
                            shouldMatchRings,
                            matchAtomType);
            this.setTimeout(CDKMCS.isTimeout());
            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);

            while (!allMaxOverlaps.empty()) {
                List<List<CDKRMap>> maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), getSource(), getTarget());
                identifyMatchedParts(maxOverlapsAtoms, getSource(), getTarget());
                allMaxOverlaps.pop();
            }
        }
        return solutions;
    }

    /**
     * This function calculates only one solution (exact) because we are looking
     * at the molecules which are exactly same in terms of the bonds and atoms
     * determined by the Fingerprint
     *
     * @param Molecule1
     * @param Molecule2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return List
     * @throws CDKException
     */
    public synchronized List<Map<Integer, Integer>> calculateIsomorphs(IAtomContainer Molecule1,
            IAtomContainer Molecule2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws CDKException {

        setSource(Molecule1);
        setTarget(Molecule2);
        List<Map<Integer, Integer>> solutions = new ArrayList<>();
        setMappings(solutions);

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);
        if ((getSource().getAtomCount() == 1) || (getTarget().getAtomCount() == 1)) {

            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(getSource(), getTarget());
            this.setTimeout(CDKMCS.isTimeout());
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                identifySingleAtomsMatchedParts(overlaps, getSource(), getTarget());
            }

        } else {

            List<List<CDKRMap>> overlaps
                    = CDKMCS.getIsomorphMaps(getSource(), getTarget(), shouldMatchBonds, shouldMatchRings, matchAtomType);
            this.setTimeout(CDKMCS.isTimeout());
            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);

            while (!allMaxOverlaps.empty()) {
                List<List<CDKRMap>> maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), getSource(), getTarget());
                identifyMatchedParts(maxOverlapsAtoms, getSource(), getTarget());
                allMaxOverlaps.pop();
            }
        }
        return solutions;
    }

    /**
     *
     * @param overlaps
     * @return removed List
     */
    protected synchronized List<List<CDKRMap>> removeSubGraph(List<List<CDKRMap>> overlaps) {

        List<List<CDKRMap>> reducedList = new ArrayList<>(overlaps);

        for (int i = 0; i < overlaps.size(); i++) {
            List<CDKRMap> graphI = overlaps.get(i);

            for (int j = i + 1; j < overlaps.size(); j++) {
                List<CDKRMap> graphJ = overlaps.get(j);

                // Gi included in Gj or Gj included in Gi then
                // reduce the irrelevant solution
                if (graphI.size() != graphJ.size()) {
                    if (isSubgraph(graphJ, graphI)) {
                        reducedList.remove(graphI);
                    } else if (isSubgraph(graphI, graphJ)) {
                        reducedList.remove(graphJ);
                    }
                }

            }
        }
        return reducedList;
    }

    /**
     *
     * @param overlaps
     * @return List removed
     */
    protected synchronized List<CDKRMap> removeRedundantMappingsForSingleAtomCase(List<CDKRMap> overlaps) {
        List<CDKRMap> reducedList = Collections.synchronizedList(new ArrayList<CDKRMap>());
        reducedList.add(overlaps.get(0));
        //reducedList.add(overlaps.get(1));
        return reducedList;
    }

    /**
     * This makes sourceAtom map1 of matching atoms out of sourceAtom map1 of
     * matching bonds as produced by the get(Subgraph|Ismorphism)Map methods.
     *
     * @param rMapList The list produced by the getMap method.
     * @param graph1 first molecule. Must not be an IQueryAtomContainer.
     * @param graph2 second molecule. May be an IQueryAtomContainer.
     * @return The mapping found projected on graph1. This is sourceAtom List of
     * CDKRMap objects containing Ids of matching atoms.
     */
    private synchronized List<List<CDKRMap>> makeAtomsMapOfBondsMap(List<CDKRMap> rMapList, IAtomContainer graph1, IAtomContainer graph2) {
        if (rMapList == null) {
            return (null);
        }
        List<List<CDKRMap>> result;
        if (rMapList.size() == 1) {
            result = makeAtomsMapOfBondsMapSingleBond(rMapList, graph1, graph2);
        } else {
            List<CDKRMap> resultLocal = new ArrayList<>();
            for (CDKRMap rMapList2 : rMapList) {
                IBond qBond = graph1.getBond(rMapList2.getId1());
                IBond tBond = graph2.getBond(rMapList2.getId2());
                IAtom[] qAtoms = BondManipulator.getAtomArray(qBond);
                IAtom[] tAtoms = BondManipulator.getAtomArray(tBond);
                for (int j = 0; j < 2; j++) {
                    List<IBond> bondsConnectedToAtom1j = graph1.getConnectedBondsList(qAtoms[j]);
                    for (IBond bondsConnectedToAtom1j1 : bondsConnectedToAtom1j) {
                        if (bondsConnectedToAtom1j1 != qBond) {
                            IBond testBond = bondsConnectedToAtom1j1;
                            for (CDKRMap rMapList1 : rMapList) {
                                IBond testBond2;
                                if ((rMapList1).getId1() == graph1.getBondNumber(testBond)) {
                                    testBond2 = graph2.getBond((rMapList1).getId2());
                                    for (int n = 0; n < 2; n++) {
                                        List<IBond> bondsToTest = graph2.getConnectedBondsList(tAtoms[n]);
                                        if (bondsToTest.contains(testBond2)) {
                                            CDKRMap map1;
                                            if (j == n) {
                                                map1 = new CDKRMap(graph1.indexOf(qAtoms[0]), graph2.indexOf(tAtoms[0]));
                                            } else {
                                                map1 = new CDKRMap(graph1.indexOf(qAtoms[1]), graph2.indexOf(tAtoms[0]));
                                            }
                                            if (!resultLocal.contains(map1)) {
                                                resultLocal.add(map1);
                                            }
                                            CDKRMap map2;
                                            if (j == n) {
                                                map2 = new CDKRMap(graph1.indexOf(qAtoms[1]), graph2.indexOf(tAtoms[1]));
                                            } else {
                                                map2 = new CDKRMap(graph1.indexOf(qAtoms[0]), graph2.indexOf(tAtoms[1]));
                                            }
                                            if (!resultLocal.contains(map2)) {
                                                resultLocal.add(map2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            result = new ArrayList<>();
            result.add(resultLocal);
        }
        return result;
    }

    /**
     * This makes atom map1 of matching atoms out of atom map1 of matching bonds
     * as produced by the get(Subgraph|Ismorphism)Map methods. Added by Asad
     * since CDK one doesn't pick up the correct changes
     *
     * @param list The list produced by the getMap method.
     * @param sourceGraph first molecule. Must not be an IQueryAtomContainer.
     * @param targetGraph second molecule. May be an IQueryAtomContainer.
     * @return The mapping found projected on sourceGraph. This is atom List of
     * CDKRMap objects containing Ids of matching atoms.
     */
    private synchronized List<List<CDKRMap>> makeAtomsMapOfBondsMapSingleBond(List<CDKRMap> list, IAtomContainer sourceGraph, IAtomContainer targetGraph) {
        if (list == null) {
            return null;
        }
        Map<IBond, IBond> bondMap = new HashMap<>(list.size());
        for (CDKRMap solBondMap : list) {
            int id1 = solBondMap.getId1();
            int id2 = solBondMap.getId2();
            IBond qBond = sourceGraph.getBond(id1);
            IBond tBond = targetGraph.getBond(id2);
            bondMap.put(qBond, tBond);
        }
        List<CDKRMap> result1 = new ArrayList<>();
        List<CDKRMap> result2 = new ArrayList<>();
        for (IBond qbond : sourceGraph.bonds()) {
            if (bondMap.containsKey(qbond)) {
                IBond tbond = bondMap.get(qbond);
                CDKRMap map00 = null;
                CDKRMap map01 = null;
                CDKRMap map10 = null;
                CDKRMap map11 = null;

                if ((qbond.getAtom(0).getSymbol().equals(tbond.getAtom(0).getSymbol()))
                        && (qbond.getAtom(1).getSymbol().equals(tbond.getAtom(1).getSymbol()))) {
                    map00 = new CDKRMap(sourceGraph.indexOf(qbond.getAtom(0)),
                            targetGraph.indexOf(tbond.getAtom(0)));
                    map11 = new CDKRMap(sourceGraph.indexOf(qbond.getAtom(1)),
                            targetGraph.indexOf(tbond.getAtom(1)));
                    if (!result1.contains(map00)) {
                        result1.add(map00);
                    }
                    if (!result1.contains(map11)) {
                        result1.add(map11);
                    }
                }
                if ((qbond.getAtom(0).getSymbol().equals(tbond.getAtom(1).getSymbol()))
                        && (qbond.getAtom(1).getSymbol().equals(tbond.getAtom(0).getSymbol()))) {
                    map01 = new CDKRMap(sourceGraph.indexOf(qbond.getAtom(0)),
                            targetGraph.indexOf(tbond.getAtom(1)));
                    map10 = new CDKRMap(sourceGraph.indexOf(qbond.getAtom(1)),
                            targetGraph.indexOf(tbond.getAtom(0)));
                    if (!result2.contains(map01)) {
                        result2.add(map01);
                    }
                    if (!result2.contains(map10)) {
                        result2.add(map10);
                    }
                }
            }
        }
        List<List<CDKRMap>> result = new ArrayList<>();
        if (result1.size() == result2.size()) {
            result.add(result1);
            result.add(result2);
        } else if (result1.size() > result2.size()) {
            result.add(result1);
        } else {
            result.add(result2);
        }
        return result;
    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected synchronized List getMaximum(List overlaps) {
        List list = null;
        int count = 0;
        for (Object o : overlaps) {
            List arrayList = (ArrayList) o;
            if (arrayList.size() > count) {
                list = arrayList;
                count = arrayList.size();
            }

        }
        return list;
    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected synchronized Stack<List<CDKRMap>> getAllMaximum(List<List<CDKRMap>> overlaps) {

        Stack<List<CDKRMap>> allMaximumMappings = null;

        int count = -1;

        for (List<CDKRMap> arrayList : overlaps) {
            //System.out.println("O size" + sourceAtom.size());

            if (arrayList.size() > count) {

                List<CDKRMap> list = new ArrayList<>(arrayList);
                count = arrayList.size();
                allMaximumMappings = new Stack<>();
                allMaximumMappings.push(list);
            } else if (arrayList.size() == count) {
                List<CDKRMap> list = new ArrayList<>(arrayList);
                count = arrayList.size();
                allMaximumMappings.push(list);
            }

        }
        return allMaximumMappings;
    }

    /**
     *
     * @param list
     * @param source
     * @param target
     */
    protected synchronized void identifyMatchedParts(List<List<CDKRMap>> list, IAtomContainer source, IAtomContainer target) {

//        List<IAtom> array1 = new ArrayList<IAtom>();
//        List<IAtom> array2 = new ArrayList<IAtom>();

        /*
         * We have serial numbers of the bonds/Atoms to delete
         * Now we will collect the actual bond/Atoms rather than
         * serial number for deletion. RonP flag check whether reactant is
         * mapped on product or Vise Versa
         *
         */
        for (List<CDKRMap> rMap : list) {
            Map<Integer, Integer> atomNumbersFromContainer = new TreeMap<>();
            for (CDKRMap rmap : rMap) {
                IAtom sourceAtom = source.getAtom(rmap.getId1());
                IAtom targetAtom = target.getAtom(rmap.getId2());

//                array1.add(sourceAtom);
//                array2.add(targetAtom);
                int indexI = source.indexOf(sourceAtom);
                int indexJ = target.indexOf(targetAtom);

                atomNumbersFromContainer.put(indexI, indexJ);
            }
            /*Added the Mapping Numbers to the FinalMapping*
             */
            getMappings().add(atomNumbersFromContainer);
        }
    }

    /**
     *
     * @param list
     * @param source
     * @param target
     */
    protected synchronized void identifySingleAtomsMatchedParts(List<CDKRMap> list,
            IAtomContainer source,
            IAtomContainer target) {

//        List<IAtom> array1 = new ArrayList<>();
//        List<IAtom> array2 = new ArrayList<>();

        /* We have serial numbers of the bonds/Atoms to delete
         * Now we will collect the actual bond/Atoms rather than
         * serial number for deletion. RonP flag check whether reactant is
         * mapped on product or Vise Versa
         */
        TreeMap<Integer, Integer> atomNumbersFromContainer = new TreeMap<>();

        for (CDKRMap rmap : list) {
            //System.err.print("Map " + o.getClass());

            IAtom sAtom = source.getAtom(rmap.getId1());
            IAtom tAtom = target.getAtom(rmap.getId2());

//            array1.add(sAtom);
//            array2.add(tAtom);
            int indexI = source.indexOf(sAtom);
            int indexJ = target.indexOf(tAtom);

            atomNumbersFromContainer.put(indexI, indexJ);

            /*Added the Mapping Numbers to the FinalMapping*
             */
            getMappings().add(atomNumbersFromContainer);

        }
    }

    /**
     *
     * @param rmaps1
     * @param rmaps2
     * @return true if condition meet else false
     */
    protected synchronized boolean isSubgraph(List<CDKRMap> rmaps1, List<CDKRMap> rmaps2) {
        //System.out.println("Entering isSubgraph.");
        List<CDKRMap> rmaps2clone = (List<CDKRMap>) ((ArrayList<CDKRMap>) rmaps2).clone();
        for (CDKRMap rmap1 : rmaps1) {
            boolean found = false;
            for (int i = 0; i < rmaps2clone.size(); ++i) {
                CDKRMap rmap2 = rmaps2clone.get(i);
                if (isSameRMap(rmap1, rmap2)) {
                    rmaps2clone.remove(i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }

        }
        return true;
    }

    /**
     *
     * @param sourceRMap sourceAtom
     * @param targetRMap targetAtom
     * @return true if condition meet else false
     */
    protected synchronized boolean isSameRMap(CDKRMap sourceRMap, CDKRMap targetRMap) {
        return sourceRMap.getId1() == targetRMap.getId1()
                && sourceRMap.getId2() == targetRMap.getId2();
    }

    /**
     * Returns mapping solutions
     *
     * @return the mappings
     */
    public synchronized List<Map<Integer, Integer>> getMappings() {
        return mappings;
    }

    /**
     * Set mapping solutions
     *
     * @param mappings the mappings to set
     */
    public synchronized void setMappings(List<Map<Integer, Integer>> mappings) {
        this.mappings = mappings;
    }

    /**
     * @return the timeout
     */
    public boolean isTimeout() {
        return timeout;
    }

    /**
     * @param timeout the timeout to set
     */
    public void setTimeout(boolean timeout) {
        this.timeout = timeout;
    }
}
