/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.graph;

import java.util.HashMap;
import java.util.Map;
import org.jgrapht.Graph;
import org.jgrapht.graph.Multigraph;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 *
 * @author Syed Asad Rahman <asad.rahman @ bioinceptionlabs.com>
 */
public class MoleculeAsGraph {

    // make class non-instantiable
    private MoleculeAsGraph() {
    }

    /**
     * Creates a molecule graph for use with jgrapht.Bond - orders can be chosen
     *
     * @param molecule the specified molecule
     * @param respectBondOrder include bond order
     * @param respectRing
     * @param atomType
     * @return a graph representing the molecule
     */
    static public Graph getMoleculeGraph(
            IAtomContainer molecule,
            boolean respectBondOrder,
            boolean respectRing,
            boolean atomType) {
        Map<Integer, StringBuilder> vertices = new HashMap<>();
        Graph graph = new Multigraph<>(StringBuilder.class);
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            IAtom atom = molecule.getAtom(i);
            int id = molecule.indexOf(atom);
            String label;
            if (atomType) {
                label = atom.getSymbol() + ":" + atom.getAtomTypeName();
            } else {
                label = atom.getSymbol();
            }

            StringBuilder node = new StringBuilder(label);
//            StringLabeledObject node = new StringLabeledObject(label + ":" + id);
            graph.addVertex(node);
            vertices.put(id, node);
        }

        for (int i = 0; i < molecule.getBondCount(); i++) {
            IBond bond = molecule.getBond(i);
            int begin = molecule.indexOf(bond.getBegin());
            int end = molecule.indexOf(bond.getEnd());
//            String label = molecule.indexOf(bond) + "";
            String label;
            if (respectBondOrder) {
                label = bond.getOrder().numeric() + "";
            } else {
                label = IBond.Order.UNSET.numeric() + "";
            }

//            System.out.println("Bond Label: " + label);
            StringBuilder node = new StringBuilder(label);
            /*
             * int order = (int) bond.getOrder(); for (int j=0; j<order; j++) {
             * graph.addEdge(bond.getAtoms()[0], bond.getAtoms()[1]); }
             */
            graph.addEdge(vertices.get(begin), vertices.get(end), node);
        }
        return graph;
    }

}
