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
package uk.ac.ebi.reactionblast.mapping;

import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.sort;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org._3pq.jgrapht.graph.SimpleGraph;
import org.openscience.cdk.AtomContainer;
import static org.openscience.cdk.CDKConstants.VISITED;
import static org.openscience.cdk.graph.BFSShortestPath.findPathBetween;
import static org.openscience.cdk.graph.MoleculeGraphs.getMoleculeGraph;
import static org.openscience.cdk.graph.PathTools.computeFloydAPSP;
import static org.openscience.cdk.graph.matrix.AdjacencyMatrix.getMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;
import uk.ac.ebi.reactionblast.tools.labelling.SignatureMoleculeLabeller;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CanonicalNumberingGenerator {

    private static synchronized void resetFlags(IAtomContainer atomContainer) {
        for (int f = 0; f < atomContainer.getAtomCount(); f++) {
            atomContainer.getAtom(f).setFlag(VISITED, false);
        }
        for (int f = 0; f < atomContainer.getBondCount(); f++) {
            atomContainer.getBond(f).setFlag(VISITED, false);
        }
    }

    private static synchronized <T extends Comparable<? super T>> List<T> asSortedList(Collection<T> c) {
        List<T> list = new ArrayList<>(c);
        sort(list);
        return list;
    }

    private final IAtomContainer atomContainer;
    private final SimpleGraph simpleGraph;
    private final int costMatrix[][];
    private final int distanceMatrix[][];
    private final List<Integer> canonicalPermutationList;
    private final List<Integer> orbitalCanonicalLabellingList;

    /**
     * Canonical labeling for the given atom container
     *
     * @param atomContainer
     */
    public CanonicalNumberingGenerator(IAtomContainer atomContainer) {
        this.atomContainer = new AtomContainer(atomContainer);
        this.costMatrix = getMatrix(this.atomContainer);
        this.distanceMatrix = computeFloydAPSP(costMatrix);
        this.simpleGraph = getMoleculeGraph(atomContainer);

        resetFlags(this.atomContainer);

        ICanonicalMoleculeLabeller molLabel
                = //Signature based canonical labelling
                new SignatureMoleculeLabeller();
        int[] canonicalPermutation = molLabel.getCanonicalPermutation(atomContainer);
        canonicalPermutationList = asList(canonicalPermutation);

        /*
         * now re-sort the elements of the list iMax.e. C>N>O 
         * etc
         */
//        Collections.reverse(canonicalPermutationList);
        Map<String, List<Label>> subLabels = new TreeMap<>();
        int counter = 0;
        for (int i : canonicalPermutationList) {
            List<Label> labels;
            IAtom atom = this.atomContainer.getAtom(i);
            if (!subLabels.containsKey(atom.getSymbol())) {
                labels = new ArrayList<>();
            } else {
                labels = subLabels.get(atom.getSymbol());
            }
            Label l = new Label();
            l.atom = atom;
            l.rank = i;
            l.postion = counter++;
            labels.add(l);
            subLabels.put(atom.getSymbol(), labels);
        }
        orbitalCanonicalLabellingList = new ArrayList<>(canonicalPermutationList.size());
        subLabels.values().stream().map((List<Label> labels) -> {
            List<Integer> l = new ArrayList<>();
            //            Collections.sort(labels, new Distance());
            labels.stream().forEach((label) -> {
                l.add(label.rank);
            });
            return l;
        }).forEach((l) -> {
            orbitalCanonicalLabellingList.addAll(l);
        });//        Collections.sort(allLabels, new Distance());
//        for (Label label : allLabels) {
//            orbitalCanonicalLabellingList.add(label.rank);
//        }
    }

    private synchronized List<Integer> asList(int[] is) {
        List<Integer> intList = new ArrayList<>();
        for (int index = 0; index < is.length; index++) {
            intList.add(index, is[index]);
        }
        return intList;
    }

    /**
     * @return the canonicalPermutation
     */
    public synchronized int[] getCanonicalPermutation() {
        int[] val = new int[canonicalPermutationList.size()];
        int index = 0;
        for (Integer i : canonicalPermutationList) {
            val[index++] = i;
        }
        return val;
    }

    /**
     * @return the orbitalCanonicalLabelling
     */
    public synchronized int[] getOrbitalCanonicalLabelling() {
        int[] val = new int[orbitalCanonicalLabellingList.size()];
        int index = 0;
        for (Integer i : orbitalCanonicalLabellingList) {
            val[index++] = i;
        }
        return val;
    }

    class Label extends Object {

        int rank;
        IAtom atom;
        int postion;

        @Override
        public synchronized String toString() {
            StringBuilder result = new StringBuilder();
            String NEW_LINE = getProperty("line.separator");
            result.append("Atom: ").append(atom.getSymbol()).append(", Rank: ").append(this.rank);
            result.append(NEW_LINE);
            return result.toString();
        }
    }

    class Distance implements java.util.Comparator<Label> {

        /**
         *
         * @param t2
         * @param t1
         * @return
         */
        @Override
        public synchronized int compare(Label t1, Label t2) {
            List<org._3pq.jgrapht.Edge> sp
                    = findPathBetween(simpleGraph, t1.atom, t2.atom);
            if (t1.atom == t2.atom) {
                return 0;
            } else if (sp.isEmpty()) {
                return 99999;
            } else if (sp.size() == 1) {
                return 1;
            }
            return -1;
        }
    }
}
