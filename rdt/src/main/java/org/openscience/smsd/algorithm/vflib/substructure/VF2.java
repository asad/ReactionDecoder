/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * 
 ** Copyright (C) 2009-2015 Kyle Lutz <kyle.r.lutz@gmail.com>
 **
 ** This file is part of chemkit. For more information see
 ** <http://www.chemkit.org>.
 **
 ** chemkit is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU Lesser General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** chemkit is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU Lesser General Public License for more details.
 **
 ** You should have received a copy of the GNU Lesser General Public License
 ** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
 **
 ******************************************************************************/
package org.openscience.smsd.algorithm.vflib.substructure;

import java.util.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class finds mapping states between query and target molecules.
 *
 *  
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class VF2 implements IResults {

    private List<AtomAtomMapping> allAtomMCS = null;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchRings;
    private final boolean shouldMatchBonds;
    private final boolean matchAtomType;

    private boolean isSubgraph = false;
    private final ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(VF2.class);

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public VF2(IAtomContainer source, IAtomContainer target, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.source = source;
        this.target = target;
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.allAtomMCS = new ArrayList<>();
        this.isSubgraph = findSubgraph();
        this.matchAtomType = matchAtomType;
    }

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     *
     * @param source
     * @param target
     */
    public VF2(IQueryAtomContainer source, IAtomContainer target) {
        this.source = source;
        this.target = target;
        this.shouldMatchRings = true;
        this.shouldMatchBonds = true;
        this.matchAtomType = true;
        allAtomMCS = new ArrayList<>();
        this.isSubgraph = findSubgraph();
    }

    /**
     * The isomorphism method returns an isomorphism between two molecular graphs using the VF2Automorphism algorithm.
     * This can be used for finding both graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter case
     * graph 'a' is the subgraph, implying a.size() < b.size(). In the case that no isomorphism is found an empty
     * mapping is returned.
     *
     *
     *
     *
     *
     *
     *
     *
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @return
     */
    private synchronized void isomorphism() {

        if (!isDead(source, target) && MoleculeInitializer.testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings, matchAtomType);
            if (!state.isDead()) {
                state.matchFirst(state, allAtomMCS);
            }
        }
    }

    /**
     * The isomorphism method returns an isomorphism between two molecular graphs using the VF2Automorphism algorithm.
     * This can be used for finding both graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter case
     * graph 'a' is the subgraph, implying a.size() < b.size(). In the case that no isomorphism is found an empty
     * mapping is returned.
     *
     *
     */
    private synchronized void isomorphisms() {

        if (!isDead(source, target) && MoleculeInitializer.testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings, matchAtomType);
            if (!state.isDead()) {
                state.matchAll(state, allAtomMCS);
            }
        }
    }

    // Returns true substructure is bigger than the target
    private synchronized boolean isDead(IAtomContainer a, IAtomContainer b) {
        return a.getAtomCount() > b.getAtomCount();
    }

    private boolean findSubgraph() {
        isomorphism();
        return !allAtomMCS.isEmpty();
    }

    private boolean findSubgraphs() {
        isomorphisms();
        return !allAtomMCS.isEmpty();
    }

    @Override
    public List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    @Override
    public AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }

    /**
     * @return the findSubgraph
     */
    public boolean isSubgraph() {
        return isSubgraph;
    }
}
