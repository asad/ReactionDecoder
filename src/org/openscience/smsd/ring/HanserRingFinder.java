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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

/**
 * Finds the Set of all Rings. This is an implementation of the algorithm
 * published in {@cdk.cite HAN96}. Some of the comments refer to pseudo code
 * fragments listed in this article. The concept is that a regular molecular
 * graph is first converted into a path graph (refer PathGraph.java),
 * i.e. a graph where the edges are actually paths. This can list several
 * nodes that are implicitly connecting the two nodes between the path
 * is formed (refer PathEdge.java).
 *
 * The paths that join source and sink node are step by step fused and the joined
 * nodes are deleted from the path graph (collapsed path). What remains is a graph
 * of paths that have the same start and endpoint and are thus rings (source=sink=ring).
 * 
 * 
 * 
 * @author Syed Asad Rahman <asad @ ebi.ac.uk> 2009-2015
 * 
 */
final public class HanserRingFinder {

    /**
     * Returns Collection of atoms in Rings based on Hanser Ring Finding method
     * @param molecule
     * @return report collected the rings
     */
    public static synchronized Collection<List<IAtom>> findRings(IAtomContainer molecule) {
        List<List<IAtom>> rings = new ArrayList<List<IAtom>>();
        PathGraph graph = new PathGraph(molecule);

        for (int i = 0; i < molecule.getAtomCount(); i++) {
            List<PathEdge> edges = graph.remove(molecule.getAtom(i));

            for (PathEdge edge : edges) {
                List<IAtom> ring = edge.getAtoms();
                rings.add(ring);
            }
        }
        return Collections.synchronizedList(new ArrayList<List<IAtom>>(rings));
    }

    /**
     * Returns CDK object Ring set based on Hanser Ring Finding method
     * @param molecule
     * @return report collected the rings
     * @throws CDKException 
     */
    public static IRingSet getRingSet(IAtomContainer molecule) throws CDKException {
        IRingSet ringSet = DefaultChemObjectBuilder.getInstance().newInstance(IRingSet.class);
        Collection<List<IAtom>> cycles = findRings(molecule);
        for (List<IAtom> ringAtoms : cycles) {
            IRing ring = molecule.getBuilder().newInstance(IRing.class);
            for (IAtom atom : ringAtoms) {
                atom.setFlag(CDKConstants.ISINRING, true);
                ring.addAtom(atom);
                for (IAtom atomNext : ringAtoms) {
                    if (!atom.equals(atomNext)) {
                        IBond bond = molecule.getBond(atom, atomNext);
                        if (bond != null) {
                            bond.setFlag(CDKConstants.ISINRING, true);
                            ring.addElectronContainer(bond);
                        }
                    }
                }
            }
            ringSet.addAtomContainer(ring);
        }
        return ringSet;
    }
    private static final Logger LOG = Logger.getLogger(HanserRingFinder.class.getName());
}
