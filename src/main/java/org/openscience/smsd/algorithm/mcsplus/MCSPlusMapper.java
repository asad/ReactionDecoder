/* Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.exception.CDKException;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.filters.PostFilter;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class acts as a handler class for MCSPlus algorithm.
 * {@link org.openscience.smsd.algorithm.mcsplus.MCSPlus}
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public final class MCSPlusMapper implements IResults {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MCSPlusMapper.class);
    private final List<AtomAtomMapping> allAtomMCS;
    private final List<Map<Integer, Integer>> allMCS;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private boolean flagExchange = false;
    private final boolean timeout;
    private final boolean DEBUG = false;
    private final AtomMatcher atomMatcher;
    private final BondMatcher bondMatcher;

    /**
     * Constructor for the MCSPlus Plus algorithm class
     *
     * @param source
     * @param target
     * @throws org.openscience.cdk.exception.CDKException
     */
    public MCSPlusMapper(IAtomContainer source, IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm) throws CDKException {
        this.source = source;
        this.target = target;
        this.atomMatcher = am;
        this.bondMatcher = bm;
        allAtomMCS = Collections.synchronizedList(new ArrayList<>());
        allMCS = Collections.synchronizedList(new ArrayList<>());
        this.timeout = searchMCS();
    }

    /**
     * Constructor for the MCSPlus Plus algorithm class
     *
     * @param source
     * @param target
     * @throws org.openscience.cdk.exception.CDKException
     */
    public MCSPlusMapper(IQueryAtomContainer source, IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm) throws CDKException {
        this.source = source;
        this.target = target;
        this.atomMatcher = am;
        this.bondMatcher = bm;
        this.allAtomMCS = Collections.synchronizedList(new ArrayList<>());
        this.allMCS = Collections.synchronizedList(new ArrayList<>());
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
            mcsplus = new MCSPlus((IQueryAtomContainer) source, target, atomMatcher, bondMatcher);
            List<List<Integer>> overlaps = mcsplus.getOverlaps();
            mappings = Collections.synchronizedList(overlaps);

        } else if (!(source instanceof IQueryAtomContainer) && source.getAtomCount() <= target.getAtomCount()) {
            mcsplus
                    = new MCSPlusBuilder().setQuery(source).setTarget(target).setAtomMatcher(atomMatcher).setBondMatcher(bondMatcher).createMCSPlus();
            List<List<Integer>> overlaps = mcsplus.getOverlaps();
            mappings = Collections.synchronizedList(overlaps);

        } else {
            flagExchange = true;
            mcsplus = new MCSPlusBuilder().setQuery(target).setTarget(source).setAtomMatcher(atomMatcher).setBondMatcher(bondMatcher).createMCSPlus();
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
        return mappings.isEmpty();
    }

    private synchronized void setAllMapping(List<Map<Integer, Integer>> solutions) {
        try {
            int bestSolSize = 0;
            for (Map<Integer, Integer> solution : solutions) {
//                System.out.println("Number of MCSPlus solution: " + solution.size());
                Map<Integer, Integer> validSolution = Collections.synchronizedSortedMap(new TreeMap<>());
//                System.out.println("validSolution " + solutions);
                validSolution.putAll(solution);

                if (validSolution.size() > bestSolSize
                        && (validSolution.size() <= source.getAtomCount()
                        && validSolution.size() <= target.getAtomCount())) {
                    bestSolSize = validSolution.size();
                    allMCS.clear();
                }
                if (validSolution.size() == bestSolSize) {
                    allMCS.add(validSolution);
                }
            }

        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        if (DEBUG) {
            System.out.println("Number of MCSPlus solution - : allMCS " + allMCS.size());
        }
    }

    private synchronized void setAllAtomMapping() {
        if (DEBUG) {
            System.out.println("setAllAtomMapping");
            System.out.println("source size " + source.getAtomCount());
            System.out.println("target size " + target.getAtomCount());
        }
        try {
            allMCS.stream().map((solution) -> {
                AtomAtomMapping atomMapping = new AtomAtomMapping(source, target);
//                System.out.println("solution " + solution);
                solution.entrySet().stream().forEach((m) -> {
                    int indexI = m.getKey();
                    int indexJ = m.getValue();
//                    System.out.println("indexI " + indexI + ", " + "indexJ " + indexJ);
                    IAtom sourceAtom = this.source.getAtom(indexI);
                    IAtom targetAtom = this.target.getAtom(indexJ);
                    atomMapping.put(sourceAtom, targetAtom);
                });
                return atomMapping;
            }).forEach((atomMapping) -> {
                allAtomMCS.add(atomMapping);
            });
        } catch (Exception e) {
            LOGGER.error(SEVERE, null, e);
        }
        if (DEBUG) {
            System.out.println("Number of MCSPlus solution - : allAtomMCS " + allAtomMCS.size());
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
