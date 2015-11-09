/*
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package org.openscience.smsd.ring;

import static java.lang.System.out;
import java.util.ArrayList;
import static java.util.Arrays.asList;
import java.util.List;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * 
 * 
 * @author Richard L. Apodaca <rapodaca at metamolecular.com> 2009-2015,
 *         Syed Asad Rahman <asad @ ebi.ac.uk> 2009-2015
 */
public class PathGraph {
    private static final Logger LOG = getLogger(PathGraph.class.getName());

    private final List<PathEdge> edges;
    private final List<IAtom> atoms;
    private final IAtomContainer mol;

    public PathGraph(IAtomContainer molecule) {
        edges = new ArrayList<>();
        atoms = new ArrayList<>();
        this.mol = molecule;

        loadEdges(molecule);
        loadNodes(molecule);
    }

    public void printPaths() {
        for (PathEdge edge : edges) {
            if (edge.isCycle()) {
                out.print("*");
            }

            for (IAtom atom : edge.getAtoms()) {
                out.print(mol.getAtomNumber(atom) + "-");
            }

            out.println();
        }
    }

    public List<PathEdge> remove(IAtom atom) {
        List<PathEdge> oldEdges = getEdges(atom);
        List<PathEdge> result = new ArrayList<>();

        for (PathEdge edge : oldEdges) {
            if (edge.isCycle()) {
                result.add(edge);
            }
        }

        oldEdges.removeAll(result);
        edges.removeAll(result);

        List<PathEdge> newEdges = spliceEdges(oldEdges);

        edges.removeAll(oldEdges);
        edges.addAll(newEdges);
        atoms.remove(atom);

        return result;
    }

    private List<PathEdge> spliceEdges(List<PathEdge> edges) {
        List<PathEdge> result = new ArrayList<>();

        for (int i = 0; i < edges.size(); i++) {
            for (int j = i + 1; j < edges.size(); j++) {
                PathEdge splice = edges.get(j).splice(edges.get(i));

                if (splice != null) {
                    result.add(splice);
                }
            }
        }

        return result;
    }

    private List<PathEdge> getEdges(IAtom atom) {
        List<PathEdge> result = new ArrayList<>();

        for (PathEdge edge : edges) {
            if (edge.isCycle()) {
                if (edge.getAtoms().contains(atom)) {
                    result.add(edge);
                }
            } else {
                if ((edge.getSource() == atom) || (edge.getTarget() == atom)) {
                    result.add(edge);
                }
            }
        }

        return result;
    }

    private void loadEdges(IAtomContainer molecule) {
        for (int i = 0; i < molecule.getBondCount(); i++) {
            IBond bond = molecule.getBond(i);
            edges.add(new PathEdge(asList(bond.getAtom(0), bond.getAtom(1))));
        }
    }

    private void loadNodes(IAtomContainer molecule) {
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            atoms.add(molecule.getAtom(i));
        }
    }
}
