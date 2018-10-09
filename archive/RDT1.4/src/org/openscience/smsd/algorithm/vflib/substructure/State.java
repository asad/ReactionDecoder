/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *                          Gilleain Torrance <gilleain.torrance@gmail.com>
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

import java.util.List;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomTypeMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultBondMatcher;

/**
 * This class finds mapping states between query and target molecules.
 *
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
// The State class represents a single state in the isomorphism detection
// algorithm. Every state uses and modifies the same SharedState object.
final class State {

    private final boolean shouldMatchBonds;
    private final boolean shouldMatchRings;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchAtomType;

    // Returns true if the state contains an isomorphism.
    public boolean isGoal() {
        return size == source.getAtomCount();
    }

    public boolean isDead() {
        return (!isMatchPossible || source.getAtomCount() > target.getAtomCount());
    }

    public boolean hasNextCandidate(Pair<Integer, Integer> candidate) {
        return candidate.getSourceAtom() != -1;
    }

    int getSize() {
        return size;
    }

    IAtomContainer getSource() {
        return source;
    }

    IAtomContainer getTarget() {
        return target;
    }

    IAtom sourceAtom(int index) {
        return source.getAtom(index);
    }

    IAtom targetAtom(int index) {
        return target.getAtom(index);
    }
    private int size;
    private int sourceTerminalSize;
    private int targetTerminalSize;
    private Pair<Integer, Integer> lastAddition;
    private SharedState sharedState;
    private final boolean ownSharedState;
    private boolean[][] matches;
    private boolean isMatchPossible = false;

    State(IAtomContainer source, IAtomContainer target,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.size = 0;
        this.sourceTerminalSize = 0;
        this.targetTerminalSize = 0;
        this.source = source;
        this.target = target;
        this.ownSharedState = true;
        this.matches = new boolean[this.source.getAtomCount()][this.target.getAtomCount()];
        this.isMatchPossible = isFeasible();

        this.lastAddition = new Pair<>(-1, -1);
        this.sharedState = new SharedState(source.getAtomCount(),
                target.getAtomCount());
        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchAtomType = matchAtomType;
    }

    State(IQueryAtomContainer source, IAtomContainer target) {
        this.size = 0;
        this.sourceTerminalSize = 0;
        this.targetTerminalSize = 0;
        this.source = source;
        this.target = target;
        this.ownSharedState = true;
        this.matches = new boolean[this.source.getAtomCount()][this.target.getAtomCount()];
        this.isMatchPossible = isFeasible();

        this.lastAddition = new Pair<>(-1, -1);
        this.sharedState = new SharedState(source.getAtomCount(),
                target.getAtomCount());
        this.shouldMatchBonds = true;
        this.shouldMatchRings = true;
        this.shouldMatchAtomType = true;
    }

    State(State state) {
        this.size = state.size;
        this.sourceTerminalSize = state.sourceTerminalSize;
        this.targetTerminalSize = state.targetTerminalSize;
        this.source = state.source;
        this.target = state.target;
        this.ownSharedState = false;
        this.matches = state.matches;
        this.lastAddition = new Pair<>(-1, -1);
        this.sharedState = state.sharedState;
        this.shouldMatchBonds = state.shouldMatchBonds;
        this.shouldMatchRings = state.shouldMatchRings;
        this.shouldMatchAtomType = state.shouldMatchAtomType;
    }

    private boolean isFeasible() {
        for (int i = 0; i < source.getAtomCount(); i++) {
            boolean flag = false;
            for (int j = 0; j < target.getAtomCount(); j++) {
                if (matcher(i, j)) {
                    this.matches[i][j] = true;
                    flag = true;
                } else {
                    this.matches[i][j] = false;
                }
            }
            if (!flag) {
                this.matches = null;
                return false;
            }
        }
        return true;
//        System.out.println("Compatibility graph " + candidates.size());
    }

    public void dispose() {
        if (this.ownSharedState) {
            if (this.sharedState != null) {
                this.sharedState = null;
            }
        }
    }

    // Returns the current isomorphism for the state in an AtomMapping
    // object.
    AtomAtomMapping getMapping() {
        AtomAtomMapping mapping = new AtomAtomMapping(source, target);

        for (int i = 0; i < size; i++) {
            mapping.put(source.getAtom(i),
                    target.getAtom(sharedState.sourceMapping[i]));
        }
        return mapping;
    }

    // Returns the next candidate pair (sourceAtom, targetAtom) to be added
    // to the state. The candidate should be checked for feasibility and then added
    // using the addPair() method.
    Pair<Integer, Integer> nextCandidate(
            Pair<Integer, Integer> lastCandidate) {
        int lastSourceAtom = lastCandidate.getSourceAtom();
        int lastTargetAtom = lastCandidate.getTargetAtom();

        int sourceSize = source.getAtomCount();
        int targetSize = target.getAtomCount();

        if (lastSourceAtom == -1) {
            lastSourceAtom = 0;
        }

        if (lastTargetAtom == -1) {
            lastTargetAtom = 0;
        } else {
            lastTargetAtom++;
        }

        if (sourceTerminalSize > size && targetTerminalSize > size) {
            while (lastSourceAtom < sourceSize
                    && (sharedState.sourceMapping[lastSourceAtom] != -1
                    || sharedState.sourceTerminalSet[lastSourceAtom] == 0)) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        } else {
            while (lastSourceAtom < sourceSize
                    && sharedState.sourceMapping[lastSourceAtom] != -1) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        }

        if (sourceTerminalSize > size && targetTerminalSize > size) {
            while (lastTargetAtom < targetSize
                    && (sharedState.targetMapping[lastTargetAtom] != -1
                    || sharedState.targetTerminalSet[lastTargetAtom] == 0)) {
                lastTargetAtom++;
            }
        } else {
            while (lastTargetAtom < targetSize
                    && sharedState.targetMapping[lastTargetAtom] != -1) {
                lastTargetAtom++;
            }
        }

        if (lastSourceAtom < sourceSize && lastTargetAtom < targetSize) {
            return new Pair<>(lastSourceAtom, lastTargetAtom);
        }

        return new Pair<>(-1, -1);
    }

    // Adds the candidate pair (sourceAtom, targetAtom) to the state. The
    // candidate pair must be feasible to add it to the state.
    void nextState(Pair<Integer, Integer> candidate) {
        size++;
        lastAddition = candidate;

        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

        if (sharedState.sourceTerminalSet[sourceAtom] < 1) {
            sharedState.sourceTerminalSet[sourceAtom] = size;
//                sourceTerminalSize++;
        }

        if (sharedState.targetTerminalSet[targetAtom] < 1) {
            sharedState.targetTerminalSet[targetAtom] = size;
//                targetTerminalSize++;
        }

        sharedState.sourceMapping[sourceAtom] = targetAtom;
        sharedState.targetMapping[targetAtom] = sourceAtom;

        List<IAtom> sourceNeighbours
                = source.getConnectedAtomsList(source.getAtom(sourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.indexOf(neighbor);
            if (sharedState.sourceTerminalSet[neighbourIndex] < 1) {
                sharedState.sourceTerminalSet[neighbourIndex] = size;
                sourceTerminalSize++;
            }
        }

        List<IAtom> targetNeighbours = target.getConnectedAtomsList(target.getAtom(targetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.indexOf(neighbor);
            if (sharedState.targetTerminalSet[neighbourIndex] < 1) {
                sharedState.targetTerminalSet[neighbourIndex] = size;
                targetTerminalSize++;
            }
        }
    }

    // Restores the shared state to how it was before adding the last
    // candidate pair. Assumes addPair() has been called on the state only once.
    void backTrack() {
        if (isGoal()) {
            lastAddition = new Pair<>(-1, -1);
            return;
        }
        int addedSourceAtom = lastAddition.getSourceAtom();

        if (sharedState.sourceTerminalSet[addedSourceAtom] == size) {
            sharedState.sourceTerminalSet[addedSourceAtom] = 0;
        }

        List<IAtom> sourceNeighbours
                = source.getConnectedAtomsList(source.getAtom(addedSourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.indexOf(neighbor);
            if (sharedState.sourceTerminalSet[neighbourIndex] == size) {
                sharedState.sourceTerminalSet[neighbourIndex] = 0;
            }
        }

        int addedTargetAtom = lastAddition.getTargetAtom();

        if (sharedState.targetTerminalSet[addedTargetAtom] == size) {
            sharedState.targetTerminalSet[addedTargetAtom] = 0;
        }

        List<IAtom> targetNeighbours
                = target.getConnectedAtomsList(target.getAtom(addedTargetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.indexOf(neighbor);
            if (sharedState.targetTerminalSet[neighbourIndex] == size) {
                sharedState.targetTerminalSet[neighbourIndex] = 0;
            }
        }

        sharedState.sourceMapping[addedSourceAtom] = -1;
        sharedState.targetMapping[addedTargetAtom] = -1;
        size--;
        lastAddition = new Pair<>(-1, -1);
    }

    boolean isMatchFeasible(Pair<Integer, Integer> candidate) {
        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

//        if (sourceNeighbours.size() > targetNeighbours.size()) {
//            return false;
//        }
//        if (!matchAtoms(source.getAtom(sourceAtom), target.getAtom(targetAtom))) {
//            return false;
//        }
        if (!this.matches[sourceAtom][targetAtom]) {
            return false;
        }

        int sourceTerminalNeighborCount = 0;
        int targetTerminalNeighborCount = 0;
        int sourceNewNeighborCount = 0;
        int targetNewNeighborCount = 0;

        List<IAtom> sourceNeighbours
                = source.getConnectedAtomsList(source.getAtom(sourceAtom));

        for (IAtom neighbour : sourceNeighbours) {
            int neighbourIndex = source.indexOf(neighbour);

            IAtom sourceAtomAtom = source.getAtom(sourceAtom);
            IBond sourceBond = source.getBond(sourceAtomAtom, neighbour);

            if (sharedState.sourceMapping[neighbourIndex] != -1) {
                int targetNeighbor = sharedState.sourceMapping[neighbourIndex];
                IAtom targetNeighbourAtom = target.getAtom(targetNeighbor);
                IAtom targetAtomAtom = target.getAtom(targetAtom);

                if (target.getBond(targetAtomAtom, targetNeighbourAtom) == null) {
                    return false;
                }

                IBond targetBond = target.getBond(targetAtomAtom, targetNeighbourAtom);
                if (!matchBonds(sourceBond, targetBond)) {
                    return false;
                }

            } else {
                if (sharedState.sourceTerminalSet[neighbourIndex] > 0) {
                    sourceTerminalNeighborCount++;
                } else {
                    sourceNewNeighborCount++;
                }
            }
        }

        List<IAtom> targetNeighbours
                = target.getConnectedAtomsList(target.getAtom(targetAtom));
        for (IAtom neighbour : targetNeighbours) {
            int neighbourIndex = target.indexOf(neighbour);
            if (sharedState.targetMapping[neighbourIndex] != -1) {
//                    int sourceNeighbor = sharedState.targetMapping[neighbourIndex];
//                    IAtom sourceNeighbourAtom = source.getAtom(sourceNeighbor);
//                    IAtom sourceAtomAtom = source.getAtom(targetAtom);
//
//                    if (source.getBond(sourceAtomAtom, sourceNeighbourAtom) == null) {
//                        return false;
//                    }
            } else {
                if (sharedState.targetTerminalSet[neighbourIndex] > 0) {
                    targetTerminalNeighborCount++;
                } else {
                    targetNewNeighborCount++;
                }
            }
        }
        return (sourceTerminalNeighborCount <= targetTerminalNeighborCount)
                && (sourceNewNeighborCount <= targetNewNeighborCount);
    }

    boolean matchFirst(State state, List<AtomAtomMapping> mappings) {
//            System.out.println("Matched " + state.size + " out of " + state.source.getAtomCount());
        if (state.isGoal()) {
            mappings.add(state.getMapping());
            return true;
        }

        Pair<Integer, Integer> lastCandidate = new Pair<>(-1, -1);

        boolean found = false;
        while (!found) {
            Pair<Integer, Integer> candidate = state.nextCandidate(lastCandidate);

            if (!state.hasNextCandidate(candidate)) {
                return false;
            }

            lastCandidate = candidate;

            if (state.isMatchFeasible(candidate)) {
                State nextState = new State(state);
                nextState.nextState(candidate);
                found = matchFirst(nextState, mappings);
                if (found) {
                    return true;
                }
                nextState.backTrack();
            }
        }

        return found;
    }

    /* TO DO: Fix the match all results*/
    void matchAll(State state, List<AtomAtomMapping> mappings) {
//        System.out.println("Matched " + state.size + " out of " + state.source.getAtomCount());

        if (state.isGoal()) {
            AtomAtomMapping map = state.getMapping();
            if (!hasMap(map, mappings)) {
                mappings.add(state.getMapping());
            }
            return;
        }

        Pair<Integer, Integer> lastCandidate = new Pair<>(-1, -1);
        Pair<Integer, Integer> candidate = state.nextCandidate(lastCandidate);

        while (state.hasNextCandidate(candidate)) {
            lastCandidate = candidate;
            if (state.isMatchFeasible(lastCandidate)) {
                State nextState = new State(state);
                nextState.nextState(candidate);
                matchAll(nextState, mappings);
                nextState.backTrack();
            }
        }
    }

    private boolean matcher(int queryAtom, int targetAtom) {
        List<IAtom> sourceNeighbours
                = source.getConnectedAtomsList(source.getAtom(queryAtom));
        List<IAtom> targetNeighbours
                = target.getConnectedAtomsList(target.getAtom(targetAtom));
        if (!matchAtoms(source.getAtom(queryAtom), target.getAtom(targetAtom))) {
            return false;
        }
        return sourceNeighbours.size() <= targetNeighbours.size();
    }

    boolean matchBonds(IBond queryBond, IBond targetBond) {
        BondMatcher defaultVFBondMatcher
                = new DefaultBondMatcher(
                        queryBond, shouldMatchBonds, shouldMatchRings, shouldMatchAtomType);
        return defaultVFBondMatcher.matches(targetBond);
    }

    boolean matchAtoms(IAtom sourceAtom, IAtom targetAtom) {
        AtomMatcher defaultVFAtomMatcher;
        if (shouldMatchAtomType) {
            defaultVFAtomMatcher
                    = new DefaultAtomTypeMatcher(sourceAtom, shouldMatchRings);
        } else {
            defaultVFAtomMatcher = new DefaultAtomMatcher(sourceAtom, shouldMatchRings);
        }
        return defaultVFAtomMatcher.matches(targetAtom);
    }

    private boolean hasMap(AtomAtomMapping map, List<AtomAtomMapping> mappings) {
        for (AtomAtomMapping test : mappings) {
            if (test.equals(map)) {
                return true;
            }
        }
        return false;
    }
}
