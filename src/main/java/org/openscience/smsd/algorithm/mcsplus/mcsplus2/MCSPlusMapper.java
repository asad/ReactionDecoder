/* Copyright (C) 2009-2018  Syed Asad Rahman <asad at ebi.ac.uk>
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
 * MERCHANTABILITY or FITNESS FOR sourceAtom PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus.mcsplus2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.filters.PostFilter;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class acts as a handler class for MCSPlus algorithm.
 * {@link org.openscience.smsd.algorithm.mcsplus.mcsplus2.MCSPlus}
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public final class MCSPlusMapper implements IResults {

    private final List<AtomAtomMapping> allAtomMCS;
    private final List<Map<Integer, Integer>> allMCS;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private boolean flagExchange = false;
    private final boolean shouldMatchRings;
    private final boolean shouldMatchBonds;
    private final boolean matchAtomType;
    private final boolean timeout;

    /**
     * Constructor for the MCS Plus algorithm class
     *
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public MCSPlusMapper(IAtomContainer source, IAtomContainer target,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.source = source;
        this.target = target;
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.matchAtomType = matchAtomType;
        allAtomMCS = Collections.synchronizedList(new ArrayList<>());
        allMCS = Collections.synchronizedList(new ArrayList<>());
        this.timeout = searchMCS();
    }

    /**
     * Constructor for the MCS Plus algorithm class
     *
     * @param source
     * @param target
     */
    public MCSPlusMapper(IQueryAtomContainer source, IAtomContainer target) {
        this.source = source;
        this.target = target;
        this.shouldMatchRings = true;
        this.shouldMatchBonds = true;
        this.matchAtomType = true;
        allAtomMCS = Collections.synchronizedList(new ArrayList<>());
        allMCS = Collections.synchronizedList(new ArrayList<>());
        this.timeout = searchMCS();
    }

    /**
     * {@inheritDoc} Function is called by the main program and serves as a
     * starting point for the comparison procedure.
     *
     */
    private synchronized boolean searchMCS() {
        List<List<Integer>> mappings;
        MCSPlus mcsplus;

        if (source instanceof IQueryAtomContainer) {
            mcsplus = new MCSPlus((IQueryAtomContainer) source, target);
            List<List<Integer>> overlaps = mcsplus.getOverlaps();
            mappings = Collections.synchronizedList(overlaps);

        } else if (!(source instanceof IQueryAtomContainer) && source.getAtomCount() < target.getAtomCount()) {
            mcsplus = new MCSPlus(source, target, shouldMatchBonds, shouldMatchRings, matchAtomType);
            List<List<Integer>> overlaps = mcsplus.getOverlaps();
            mappings = Collections.synchronizedList(overlaps);

        } else {
            flagExchange = true;
            mcsplus = new MCSPlus(target, source, shouldMatchBonds, shouldMatchRings, matchAtomType);
            List<List<Integer>> overlaps = mcsplus.getOverlaps();
            mappings = Collections.synchronizedList(overlaps);
        }
        if (flagExchange) {
            mappings = reverseMappings(mappings);
        }
//        System.out.println("PreFilter.filter " + mappings);
        List<Map<Integer, Integer>> solutions = PostFilter.filter(mappings);
//        System.out.println("PostFilter.filter " + solutions);
        setAllMapping(solutions);
        setAllAtomMapping();
        return mcsplus.isTimeout();
    }

    private synchronized void setAllMapping(List<Map<Integer, Integer>> solutions) {
        try {
            int counter = 0;
            int bestSolSize = 0;
            for (Map<Integer, Integer> solution : solutions) {
//                System.out.println("Number of MCS solution: " + solution);
                Map<Integer, Integer> validSolution = Collections.synchronizedSortedMap(new TreeMap<>());

                solution.entrySet().stream().forEach((map) -> {
                    validSolution.put(map.getKey(), map.getValue());
                });

                if (validSolution.size() > bestSolSize) {
                    bestSolSize = validSolution.size();
                    counter = 0;
                    allMCS.clear();
                }
                if (validSolution.size() == bestSolSize) {
                    allMCS.add(counter++, validSolution);
                }
            }

        } catch (Exception ex) {
        }
    }

    private synchronized void setAllAtomMapping() {
        try {

            int counter = 0;
            for (Map<Integer, Integer> solution : allMCS) {
                AtomAtomMapping atomMappings = new AtomAtomMapping(source, target);
                solution.entrySet().forEach((map) -> {
                    int IIndex = map.getKey();
                    int JIndex = map.getValue();

                    IAtom sourceAtom;
                    IAtom targetAtom;

                    sourceAtom = source.getAtom(IIndex);
                    targetAtom = target.getAtom(JIndex);
                    atomMappings.put(sourceAtom, targetAtom);
                });
                allAtomMCS.add(counter++, atomMappings);
            }
        } catch (Exception I) {
            I.getCause();
        }

    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }

    /**
     * @return the timeout
     */
    public synchronized boolean isTimeout() {
        return timeout;
    }

    private List<List<Integer>> reverseMappings(List<List<Integer>> mappings) {
//        System.out.println("Before reverse " + mappings);
        List<List<Integer>> reverse = new ArrayList<>();
        mappings.stream().map((mapping) -> {
            Collections.reverse(mapping);
            return mapping;
        }).forEach((mapping) -> {
            reverse.add(mapping);
        });

//        System.out.println("reverse " + reverse);
        return reverse;
    }
}
