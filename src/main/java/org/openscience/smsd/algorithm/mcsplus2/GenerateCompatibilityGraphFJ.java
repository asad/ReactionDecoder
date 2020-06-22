/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.RecursiveTask;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.AtomBondMatcher;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.graph.Edge;
import org.openscience.smsd.helper.LabelContainer;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class GenerateCompatibilityGraphFJ extends RecursiveTask<List<Result>> {

    private final static ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(GenerateCompatibilityGraphFJ.class);
    private final static boolean DEBUG = false;
    private static final int THRESHOLD = 20;
    private static final int COMPLEX_MAX_GRAPH_NODE_COUNT = 100;
    private final int startIndex;
    private final int endIndex;

    private final IAtomContainer source;
    private final IAtomContainer target;
    private final AtomMatcher atomMatcher;
    private final BondMatcher bondMatcher;

    /**
     *
     * @param startIndex
     * @param endIndex
     * @param source
     * @param target
     * @param atomMatcher
     * @param bondMatcher
     */
    public GenerateCompatibilityGraphFJ(int startIndex,
            int endIndex,
            IAtomContainer source,
            IAtomContainer target,
            AtomMatcher atomMatcher,
            BondMatcher bondMatcher) {
        this.endIndex = endIndex;
        this.source = source;
        this.target = target;
        this.startIndex = startIndex;
        this.atomMatcher = atomMatcher;
        this.bondMatcher = bondMatcher;
    }

    @Override

    protected List<Result> compute() {

        if (endIndex - startIndex < THRESHOLD) {
            List<Result> arrayList = new ArrayList<>();
            arrayList.add(processing(startIndex, endIndex));
            return new ArrayList<>(new HashSet<>(arrayList));//remove any duplicates
        } else {

            if (DEBUG) {
                System.out.println("Splitting workLoad startIndex: " + startIndex + ", endIndex: " + endIndex);
            }

            List<GenerateCompatibilityGraphFJ> subtasks
                    = new ArrayList<>();
            subtasks.addAll(createSubtasks());

            //Collection<CustomRecursiveTask> invokeAll = invokeAll(subtasks);
            subtasks.forEach((subtask) -> {
                subtask.fork();
            });

            List<Result> result = new ArrayList<>();
            subtasks.forEach((subtask) -> {
                result.addAll(subtask.join());
            });
            return new ArrayList<>(new HashSet<>(result));//remove any duplicates;
        }
    }

    private List<GenerateCompatibilityGraphFJ> createSubtasks() {
        List<GenerateCompatibilityGraphFJ> dividedTasks = new ArrayList<>();
        int middle = (endIndex + startIndex) / 2;

        GenerateCompatibilityGraphFJ partOne = new GenerateCompatibilityGraphFJ(startIndex, middle, source, target, atomMatcher, bondMatcher);
        GenerateCompatibilityGraphFJ partTwo = new GenerateCompatibilityGraphFJ(middle, endIndex, source, target, atomMatcher, bondMatcher);
        dividedTasks.add(partOne);
        dividedTasks.add(partTwo);

        return dividedTasks;
    }

    private Result processing(int startIndex, int endIndex) {
        Result result;
        if (DEBUG) {
            System.out.println(" GenerateCompatibilityGraphFJ ");
            System.out.println("Splitting workLoad startIndex: " + startIndex + ", endIndex: " + endIndex);
        }

        if (source.getAtomCount() > COMPLEX_MAX_GRAPH_NODE_COUNT
                || target.getAtomCount() > COMPLEX_MAX_GRAPH_NODE_COUNT) {
            result = new Result();
            if (DEBUG) {
                System.out.println("CASE LARGE GRAPH");
            }
            List<Integer> compGraphNodesCZero = new ArrayList<>(); //Initialize the compGraphNodesCZero List
            compatibilityGraphNodesIfCEdgeIsZero(startIndex, endIndex, result, compGraphNodesCZero);
            compatibilityGraphCEdgeZero(result, compGraphNodesCZero);
            compGraphNodesCZero.clear();
        } else {
            result = new Result();
            if (DEBUG) {
                System.out.println("Calling Compatibility Graph Nodes ");
            }
            compatibilityGraphNodes(startIndex, endIndex, result);
            if (DEBUG) {
                System.out.println("Calling Compatibility Graph ");
            }
            compatibilityGraph(result);
            if (DEBUG) {
                System.out.println("c-edges " + result.cEdges.size());
            }
            if (DEBUG) {
                System.out.println("d-edges " + result.dEdges.size());
            }

            if (result.cEdges.isEmpty()) {
                result = new Result();
                List<Integer> compGraphNodesCZero = new ArrayList<>(); //Initialize the compGraphNodesCZero List
                compatibilityGraphNodesIfCEdgeIsZero(startIndex, endIndex, result, compGraphNodesCZero);
                compatibilityGraphCEdgeZero(result, compGraphNodesCZero);
                compGraphNodesCZero.clear();
            }
        }
        return result;
    }

    /**
     * compGraphNodesCZero is used to build up of the edges of the compatibility
     * graph
     *
     * @return
     * @throws IOException
     */
    private Integer compatibilityGraphNodesIfCEdgeIsZero(int startIndex, int endIndex, Result result, List<Integer> compGraphNodesCZero) {

        int count_nodes = 1;
        List<String> list = new ArrayList<>();
        LabelContainer labelContainer = LabelContainer.getInstance();

        for (int i = startIndex; i < endIndex; i++) {
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
                    result.compGraphNodes.add(i);
                    result.compGraphNodes.add(j);
                    result.compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                } else if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())
                        && !list.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    result.compGraphNodes.add(i);
                    result.compGraphNodes.add(j);
                    result.compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                }
            }
        }
        list.clear();
        if (DEBUG) {
            System.out.println("count_nodes " + count_nodes);
        }
        return count_nodes;
    }

    /**
     * compatibilityGraphCEdgeZero is used to build up of the edges of the
     * compatibility graph BIS
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphCEdgeZero(Result result, List<Integer> compGraphNodesCZero) {

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
                        addZeroEdges(result.cEdges, result.dEdges, reactantBond, productBond, a, b);
                    } //                    else if (reactantBond == null && productBond == null
                    //                            && ((source.getAtomCount() < (COMPLEX_MAX_GRAPH_NODE_COUNT)
                    //                            && target.getAtomCount() < (COMPLEX_MAX_GRAPH_NODE_COUNT))
                    //                            || (result.dEdges.size() < result.compGraphNodes.size()))) {
                    //                        //50 unique condition to speed up the AAM
                    //                        Edge edge = new Edge((a / 4) + 1, (b / 4) + 1);
                    //                        if (!result.dEdges.contains(edge)) {
                    //                            result.dEdges.add(edge);
                    //                        }
                    //                    }
                    else if (reactantBond == null && productBond == null) {
                        //50 unique condition to speed up the AAM
                        Edge edge = new Edge(((a / 4) + 1), ((b / 4) + 1));
                        if (!result.dEdges.contains(edge)) {
                            result.dEdges.add(edge);
                        }
                    }
                }
            }
        }
        if (DEBUG) {
            //Size of C and D edges of the compatibility graph
            int cEdgesSize = result.cEdges.size();
            int dEdgesSize = result.dEdges.size();
            System.out.println("cEdgesSize " + cEdgesSize);
            System.out.println("dEdgesSize " + dEdgesSize);
        }
        return 0;
    }

    private void addZeroEdges(List<Edge> cEdges, List<Edge> dEdges,
            IBond reactantBond, IBond productBond,
            int indexI, int indexJ) {
        if (AtomBondMatcher.matchAtomAndBond(reactantBond, productBond, atomMatcher, bondMatcher, true)) {
            Edge edge = new Edge(((indexI / 4) + 1), ((indexJ / 4) + 1));
            if (!cEdges.contains(edge)) {
                cEdges.add(edge);
            }
        } else {
            Edge edge = new Edge(((indexI / 4) + 1), ((indexJ / 4) + 1));
            if (!dEdges.contains(edge)) {
                dEdges.add(edge);
            }
        }
    }

    private Map<IAtom, List<String>> labelAtomsBySymbol(int startIndex, int endIndex, IAtomContainer atomCont) {
        Map<IAtom, List<String>> label_list = new HashMap<>();

        for (int i = startIndex; i < endIndex; i++) {
            List<String> label = new ArrayList<>(7);
            for (int a = 0; a < 7; a++) {
                label.add(a, "Z9");
            }

            IAtom refAtom = atomCont.getAtom(i);
            if (refAtom == null) {
                return label_list;
            }
            /*
             * Important Step: Discriminate between source atom types
             */
            String referenceAtom;
            if (refAtom instanceof IQueryAtom) {
                referenceAtom = ((IQueryAtom) refAtom).getSymbol() == null ? "*" : ((IQueryAtom) refAtom).getSymbol();
                if (DEBUG) {
                    System.out.println("referenceAtom " + referenceAtom);
                }
            } else {
                referenceAtom = refAtom.getSymbol(); //+ refAtom.getAtomicNumber();
            }
            label.set(0, referenceAtom);
            List<IAtom> connAtoms = atomCont.getConnectedAtomsList(refAtom);

            int counter = 1;

            for (IAtom negAtom : connAtoms) {
                String neighbouringAtom;
                if (refAtom instanceof IQueryAtom) {
                    neighbouringAtom = ((IQueryAtom) negAtom).getSymbol() == null ? "*" : ((IQueryAtom) negAtom).getSymbol();
//                    System.out.println("neighbouringAtom " + neighbouringAtom);
                } else {
                    neighbouringAtom = negAtom.getSymbol(); //+ negAtom.getAtomicNumber();
                }
                label.set(counter, neighbouringAtom);
                counter += 1;
            }
            if (DEBUG) {
                System.out.println("label " + label);
            }
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
    private int compatibilityGraphNodes(int startIndex, int endIndex, Result result) {

        Set<Edge> edges = new HashSet<>();

        int nodeCount = 1;
        Map<IAtom, List<String>> labelAtomsBySymbolA = labelAtomsBySymbol(startIndex, endIndex, source);
        Map<IAtom, List<String>> labelAtomsBySymbolB = labelAtomsBySymbol(0, target.getAtomCount(), target);

        for (Map.Entry<IAtom, List<String>> labelA : labelAtomsBySymbolA.entrySet()) {
            if (DEBUG) {
                System.out.println("labelA.getValue() " + labelA.getValue());
            }
            for (Map.Entry<IAtom, List<String>> labelB : labelAtomsBySymbolB.entrySet()) {
                IAtom atom = labelA.getKey();
                if (((atom instanceof IQueryAtom) && ((IQueryAtom) atom).matches(labelB.getKey()))
                        || (!(atom instanceof IQueryAtom) && atom.getSymbol().equals(labelB.getKey().getSymbol()))) {
                    if (DEBUG) {
                        System.out.println("labelB.getValue() " + labelB.getValue());
                    }
                    int atomNumberI = source.indexOf(labelA.getKey());
                    int atomNumberJ = target.indexOf(labelB.getKey());
                    Edge e = new Edge(atomNumberI, atomNumberJ);
                    if (!edges.contains(e)) {
                        edges.add(e);
                        result.compGraphNodes.add(atomNumberI);
                        result.compGraphNodes.add(atomNumberJ);
                        result.compGraphNodes.add(nodeCount);
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
    private int compatibilityGraph(Result result) {
        int comp_graph_nodes_List_size = result.compGraphNodes.size();
        if (DEBUG) {
            System.out.println("Source atom count " + source.getAtomCount());
            System.out.println("target atom count " + target.getAtomCount());
            System.out.println("Expected " + (source.getAtomCount() * target.getAtomCount())
                    + " Found Compatibilty: " + ((result.compGraphNodes.size() / 3) * 2));
            System.out.println("compGraphNodes " + result.compGraphNodes);
        }
        for (int a = 0; a < comp_graph_nodes_List_size; a += 3) {
            for (int b = a; b < comp_graph_nodes_List_size; b += 3) {
                if ((a != b)
                        && (!Objects.equals(result.compGraphNodes.get(a), result.compGraphNodes.get(b)))
                        && (!Objects.equals(result.compGraphNodes.get(a + 1), result.compGraphNodes.get(b + 1)))) {

                    IBond reactantBond;
                    IBond productBond;
                    if (DEBUG) {
                        System.out.println("a " + result.compGraphNodes.get(a) + " b " + result.compGraphNodes.get(b));
                    }//exists a bond in molecule 2, so that molecule 1 pair is connected?
                    reactantBond = source.getBond(source.getAtom(result.compGraphNodes.get(a)), source.getAtom(result.compGraphNodes.get(b)));
                    productBond = target.getBond(target.getAtom(result.compGraphNodes.get(a + 1)), target.getAtom(result.compGraphNodes.get(b + 1)));

                    boolean connectedFlag = false;
                    boolean disConnectedFlag = false;
                    boolean matchBondFlag = false;

                    if (reactantBond != null
                            && productBond != null) {
                        connectedFlag = true;
                    }

                    if (reactantBond == null
                            && productBond == null) {
                        disConnectedFlag = true;
                    }

                    if (connectedFlag
                            && AtomBondMatcher.
                                    matchAtomAndBond(reactantBond, productBond, atomMatcher, bondMatcher, true)) {
                        matchBondFlag = true;
                    }

                    //in case that both molecule pairs are connected a c-edge is generated
                    if (connectedFlag && matchBondFlag) {
                        Edge edge = new Edge(((a / 3) + 1), ((b / 3) + 1));
                        result.cEdges.add(edge);
                    }

                    //in case that both molecule pairs are not connected a d-edge is generated
                    if (disConnectedFlag) {
                        Edge edge = new Edge(((a / 3) + 1), ((b / 3) + 1));
                        result.dEdges.add(edge);
                    }
                }
            }
        }
        if (DEBUG) {
            int cEdgesSize = result.cEdges.size();
            int dEdgesSize = result.dEdges.size();
        }
        return 0;
    }
}
