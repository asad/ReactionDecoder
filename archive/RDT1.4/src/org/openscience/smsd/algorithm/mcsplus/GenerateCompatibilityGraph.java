/*
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ebi.ac.uk>
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
 * You should have received iIndex copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.smsd.algorithm.matchers.DefaultMatcher;
import org.openscience.smsd.helper.LabelContainer;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class GenerateCompatibilityGraph implements Serializable {

    private static final long serialVersionUID = 96986606860861L;
    private List<Integer> compGraphNodes = null;
    private List<Integer> compGraphNodesCZero = null;
    private final List<Integer> cEdges;
    private final List<Integer> dEdges;
    private int cEdgesSize = 0;
    private int dEdgesSize = 0;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchBonds;
    private final boolean shouldMatchRings;
    private final boolean matchAtomType;

    /**
     * Generates a compatibility graph between two molecules
     *
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @throws java.io.IOException
     */
    public GenerateCompatibilityGraph(
            IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws IOException {
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.matchAtomType = matchAtomType;
        this.source = source;
        this.target = target;
        compGraphNodes = new ArrayList<>();
        compGraphNodesCZero = new ArrayList<>();
        cEdges = Collections.synchronizedList(new ArrayList<Integer>());
        dEdges = Collections.synchronizedList(new ArrayList<Integer>());

        /*
         Generate all possible graphs when no ring match or atom type is required
         */
        /*
         Modification for AAM only
         */
        if ((!shouldMatchBonds || !matchAtomType)
                && source.getAtomCount() > 30 && target.getAtomCount() > 30) {
            compatibilityGraphNodesIfCEdgeIsZero();
            compatibilityGraphCEdgeZero();
            clearCompGraphNodesCZero();
        } else {
//        System.out.println("compatibilityGraphNodes ");
            compatibilityGraphNodes();
//        System.out.println("compatibilityGraph ");
            compatibilityGraph();
//        System.out.println("c-edges " + getCEgdes().size());
//        System.out.println("d-edges " + getDEgdes().size());

            if (getCEdgesSize() == 0) {
                clearCompGraphNodes();

                clearCEgdes();
                clearDEgdes();

                resetCEdgesSize();
                resetDEdgesSize();

                compatibilityGraphNodesIfCEdgeIsZero();
                compatibilityGraphCEdgeZero();
                clearCompGraphNodesCZero();
            }
        }
    }

    private Map<IAtom, List<String>> labelAtomsBySymbol(IAtomContainer atomCont) {
        Map<IAtom, List<String>> label_list = new HashMap<>();

        for (int i = 0; i < atomCont.getAtomCount(); i++) {
            List<String> label = new ArrayList<>(7);
            for (int a = 0; a < 7; a++) {
                label.add(a, "Z9");
            }

            IAtom refAtom = atomCont.getAtom(i);
            /*
             * Important Step: Discriminate between source atom types
             */
            String referenceAtom;
            if (refAtom instanceof IQueryAtom) {
                referenceAtom = ((IQueryAtom) refAtom).getSymbol() == null ? "*" : ((IQueryAtom) refAtom).getSymbol();
//                System.out.println("referenceAtom " + referenceAtom);
            } else if (!(refAtom instanceof IQueryAtom) && this.matchAtomType) {
                referenceAtom = refAtom.getAtomTypeName() == null ? refAtom.getSymbol() : refAtom.getAtomTypeName();
            } else {
                referenceAtom = refAtom.getSymbol();
            }
            label.set(0, referenceAtom);
            List<IAtom> connAtoms = atomCont.getConnectedAtomsList(refAtom);

            int counter = 1;

            for (IAtom negAtom : connAtoms) {
                String neighbouringAtom;
                if (refAtom instanceof IQueryAtom) {
                    neighbouringAtom = ((IQueryAtom) negAtom).getSymbol() == null ? "*" : ((IQueryAtom) negAtom).getSymbol();
//                    System.out.println("neighbouringAtom " + neighbouringAtom);
                } else if (!(negAtom instanceof IQueryAtom) && this.matchAtomType) {
                    neighbouringAtom = negAtom.getAtomTypeName() == null ? negAtom.getSymbol() : negAtom.getAtomTypeName();
                } else {
                    neighbouringAtom = negAtom.getSymbol();
                }
                label.set(counter, neighbouringAtom);
                counter += 1;
            }
//            System.out.println("label " + label);
            bubbleSort(label);
            label_list.put(refAtom, label);
        }
        return label_list;
    }

    private void bubbleSort(List<String> num) {
        int j;
        boolean flag = true;   // set flag to true to begin first pass
        String temp;   //holding variable

        while (flag) {
            flag = false;    //set flag to false awaiting a possible swap
            for (j = 0; j < (num.size() - 1); j++) {
                if (num.get(j).compareTo(num.get(j + 1)) > 0) // change to < for descending sort
                {
                    temp = num.get(j);                //swap elements
                    num.set(j, num.get(j + 1));
                    num.set(j + 1, temp);
                    flag = true;              //shows a swap occurred  
                }
            }
        }
    }

    /**
     * Generate Compatibility Graph Nodes
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphNodes() throws IOException {

        compGraphNodes.clear();

        Set<Edge> edges = new HashSet<>();

        int nodeCount = 1;
        Map<IAtom, List<String>> labelAtomsBySymbolA = labelAtomsBySymbol(source);
        Map<IAtom, List<String>> labelAtomsBySymbolB = labelAtomsBySymbol(target);

        for (Map.Entry<IAtom, List<String>> labelA : labelAtomsBySymbolA.entrySet()) {
//            System.err.println("labelA.getValue() " + labelA.getValue());
            for (Map.Entry<IAtom, List<String>> labelB : labelAtomsBySymbolB.entrySet()) {
                IAtom atom = labelA.getKey();
                if (((atom instanceof IQueryAtom) && ((IQueryAtom) atom).matches(labelB.getKey()))
                        || (!(atom instanceof IQueryAtom) && atom.getSymbol().equals(labelB.getKey().getSymbol()))) {
//                        System.err.println("labelB.getValue() " + labelB.getValue());
                    int atomNumberI = source.indexOf(labelA.getKey());
                    int atomNumberJ = target.indexOf(labelB.getKey());
                    Edge e = new Edge(atomNumberI, atomNumberJ);
                    if (!edges.contains(e)) {
                        edges.add(e);
                        compGraphNodes.add(atomNumberI);
                        compGraphNodes.add(atomNumberJ);
                        compGraphNodes.add(nodeCount);
                        nodeCount += 1;
                    }
                }
            }
        }
        return 0;
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraph() throws IOException {
        int comp_graph_nodes_List_size = compGraphNodes.size();
//        System.out.println("Source atom count " + source.getAtomCount());
//        System.out.println("target atom count " + target.getAtomCount());
//        System.out.println("Expected " + (source.getAtomCount() * target.getAtomCount())
//                + " Found Compatibilty: " + ((compGraphNodes.size() / 3) * 2));
//        System.out.println("compGraphNodes " + compGraphNodes);
        for (int a = 0; a < comp_graph_nodes_List_size; a += 3) {
            for (int b = a; b < comp_graph_nodes_List_size; b += 3) {
                if ((a != b)
                        && (!Objects.equals(compGraphNodes.get(a), compGraphNodes.get(b)))
                        && (!Objects.equals(compGraphNodes.get(a + 1), compGraphNodes.get(b + 1)))) {

                    IBond reactantBond;
                    IBond productBond;

//                    System.out.println("a " + compGraphNodes.get(a) + " b " + compGraphNodes.get(b));
                    //exists a bond in molecule 2, so that molecule 1 pair is connected?
                    reactantBond = source.getBond(source.getAtom(compGraphNodes.get(a)), source.getAtom(compGraphNodes.get(b)));
                    productBond = target.getBond(target.getAtom(compGraphNodes.get(a + 1)), target.getAtom(compGraphNodes.get(b + 1)));

                    if (reactantBond != null && productBond != null) {
                        addEdges(reactantBond, productBond, a, b);
                    } else if (reactantBond == null && productBond == null) {
                        dEdges.add((a / 3) + 1);
                        dEdges.add((b / 3) + 1);
                    }
                }
            }
        }
        cEdgesSize = cEdges.size();
        dEdgesSize = dEdges.size();
        return 0;
    }

    private void addEdges(IBond reactantBond, IBond productBond, int iIndex, int jIndex) {

        if (!isMatchBond() && !isMatchRings() && !matchAtomType) {
            if (isRawMatch(reactantBond, productBond)) {
                cEdges.add((iIndex / 3) + 1);
                cEdges.add((jIndex / 3) + 1);
            }
        } else if (isMatchFeasible(reactantBond, productBond, isMatchBond(), isMatchRings(), matchAtomType)) {
            cEdges.add((iIndex / 3) + 1);
            cEdges.add((jIndex / 3) + 1);
        } else {
            dEdges.add((iIndex / 3) + 1);
            dEdges.add((jIndex / 3) + 1);
        }
    }

    private boolean isRawMatch(IBond reactantBond, IBond productBond) {
        if (reactantBond.getAtom(0).getSymbol().equals(productBond.getAtom(0).getSymbol())
                && reactantBond.getAtom(1).getSymbol().equals(productBond.getAtom(1).getSymbol())) {
            if (reactantBond.getAtom(0).getHybridization().equals(productBond.getAtom(0).getHybridization())
                    && reactantBond.getAtom(1).getHybridization().equals(productBond.getAtom(1).getHybridization())) {
                return true;
            }
        } else if (reactantBond.getAtom(1).getSymbol().equals(productBond.getAtom(0).getSymbol())
                && reactantBond.getAtom(0).getSymbol().equals(productBond.getAtom(1).getSymbol())) {
            if (reactantBond.getAtom(1).getHybridization().equals(productBond.getAtom(0).getHybridization())
                    && reactantBond.getAtom(0).getHybridization().equals(productBond.getAtom(1).getHybridization())) {
                return true;
            }
        }
        return false;
    }

    /**
     * compGraphNodesCZero is used to build up of the edges of the compatibility
     * graph
     *
     * @return
     * @throws IOException
     */
    private Integer compatibilityGraphNodesIfCEdgeIsZero() throws IOException {

        int count_nodes = 1;
        List<String> list = new ArrayList<>();
        compGraphNodesCZero = new ArrayList<>(); //Initialize the compGraphNodesCZero List
        LabelContainer labelContainer = LabelContainer.getInstance();
        compGraphNodes.clear();

        for (int i = 0; i < source.getAtomCount(); i++) {
            for (int j = 0; j < target.getAtomCount(); j++) {
                IAtom atom1 = source.getAtom(i);
                IAtom atom2 = target.getAtom(j);

                //You can also check object equal or charge, hydrogen count etc
                if ((atom1 instanceof IQueryAtom)
                        && ((IQueryAtom) atom1).matches(atom2)
                        && !list.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom2.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    compGraphNodes.add(i);
                    compGraphNodes.add(j);
                    compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                } else if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())
                        && !list.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    compGraphNodes.add(i);
                    compGraphNodes.add(j);
                    compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                }
            }
        }
        list.clear();
        return count_nodes;
    }

    /**
     * compatibilityGraphCEdgeZero is used to build up of the edges of the
     * compatibility graph BIS
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphCEdgeZero() throws IOException {

        int compGraphNodesCZeroListSize = compGraphNodesCZero.size();

        for (int a = 0; a < compGraphNodesCZeroListSize; a += 4) {
            int index_a = compGraphNodesCZero.get(a);
            int index_aPlus1 = compGraphNodesCZero.get(a + 1);
            for (int b = a + 4; b < compGraphNodesCZeroListSize; b += 4) {
                int index_b = compGraphNodesCZero.get(b);
                int index_bPlus1 = compGraphNodesCZero.get(b + 1);

                // if element atomCont !=jIndex and atoms on the adjacent sides of the bonds are not equal
                if ((a != b) && (index_a != index_b)
                        && (index_aPlus1 != index_bPlus1)) {

                    IBond reactantBond;
                    IBond productBond;

                    reactantBond = source.getBond(source.getAtom(index_a), source.getAtom(index_b));
                    productBond = target.getBond(target.getAtom(index_aPlus1), target.getAtom(index_bPlus1));

                    if (reactantBond != null && productBond != null) {
                        addZeroEdges(reactantBond, productBond, a, b);
                    } else if (reactantBond == null && productBond == null
                            && source.getAtomCount() < 50 && target.getAtomCount() < 50) {
                        //50 unique condition to speed up the AAM
                        dEdges.add((a / 4) + 1);
                        dEdges.add((b / 4) + 1);
                    }
                }
            }
        }

        //Size of C and D edges of the compatibility graph
        cEdgesSize = cEdges.size();
        dEdgesSize = dEdges.size();
        return 0;
    }

    private void addZeroEdges(IBond reactantBond, IBond productBond, int indexI, int indexJ) {
        if (isMatchFeasible(reactantBond, productBond, isMatchBond(), isMatchRings(), matchAtomType)) {
            cEdges.add((indexI / 4) + 1);
            cEdges.add((indexJ / 4) + 1);
        } else {
            dEdges.add((indexI / 4) + 1);
            dEdges.add((indexJ / 4) + 1);
        }
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @return
     */
    private boolean isMatchFeasible(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {

        if (bondA1 instanceof IQueryBond) {
            if (((IQueryBond) bondA1).matches(bondA2)) {
                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
                return atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0));
            }
            return false;
        } else {
            /*
             This one also matches atom type, not just symbols
             */
            return DefaultMatcher.matches(bondA1, bondA2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        }
    }

    public synchronized List<Integer> getCEgdes() {
        return Collections.synchronizedList(cEdges);
    }

    public synchronized List<Integer> getDEgdes() {
        return Collections.synchronizedList(dEdges);
    }

    public synchronized List<Integer> getCompGraphNodes() {
        return Collections.synchronizedList(compGraphNodes);
    }

    protected synchronized int getCEdgesSize() {
        return cEdgesSize;
    }

    protected synchronized int getDEdgesSize() {
        return dEdgesSize;
    }

    private List<Integer> getCompGraphNodesCZero() {
        return Collections.unmodifiableList(compGraphNodesCZero);
    }

    private void clearCEgdes() {
        cEdges.clear();
    }

    private void clearDEgdes() {
        dEdges.clear();
    }

    private void clearCompGraphNodes() {
        compGraphNodes.clear();
    }

    private void clearCompGraphNodesCZero() {
        compGraphNodesCZero.clear();
    }

    private void resetCEdgesSize() {
        cEdgesSize = 0;
    }

    private void resetDEdgesSize() {
        dEdgesSize = 0;
    }

    public synchronized void clear() {
        cEdges.clear();
        dEdges.clear();
        compGraphNodes.clear();
        compGraphNodesCZero.clear();
    }

    /**
     * @return the shouldMatchBonds
     */
    private boolean isMatchBond() {
        return shouldMatchBonds;
    }

    /**
     * @return the shouldMatchRings
     */
    private boolean isMatchRings() {
        return shouldMatchRings;
    }
}
