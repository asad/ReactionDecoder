/* Copyright (C) 2009-2017  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus1;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.filters.PostFilter;
import org.openscience.smsd.interfaces.IResults;

/**
 * This class acts as a handler class for MCSPlus algorithm.
 * {@link org.openscience.smsd.algorithm.mcsplus.MCSPlus}
 *
 *
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class MCSPlusMapper implements IResults {

    private final List<AtomAtomMapping> allAtomMCS;
    private final List<Map<Integer, Integer>> allMCS;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private boolean flagExchange = false;
    private final boolean timeout;
    private boolean shouldMatchBonds;
    private boolean shouldMatchRings;
    private boolean matchAtomType;

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

        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomType = matchAtomType;

        allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        allMCS = Collections.synchronizedList(new ArrayList<Map<Integer, Integer>>());
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
        allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        allMCS = Collections.synchronizedList(new ArrayList<Map<Integer, Integer>>());
        this.timeout = searchMCS();
    }

    /**
     * {@inheritDoc} Function is called by the main program and serves as a
     * starting point for the comparison procedure.
     *
     */
    private synchronized boolean searchMCS() {
        List<List<Integer>> mappings = new ArrayList<>();

        if (source instanceof IQueryAtomContainer || target instanceof IQueryAtomContainer) {
            new CDKException("Not supported");

        } else if (source.getAtomCount() >= target.getAtomCount()) {
            this.flagExchange = false;
            MoleculeHandler file1 = new MoleculeHandler(source);
            MoleculeHandler file2 = new MoleculeHandler(target);

            MCS mcs = new MCS(
                    file1.getAtomNumber(),
                    file2.getAtomNumber(),
                    file1.getStartHatom_num(),
                    file2.getStartHatom_num(),
                    file1.getBondNumber(),
                    file2.getBondNumber(),
                    file1.getAtomString(),
                    file2.getAtomString(),
                    file1.intTable,
                    file2.intTable,
                    file1.charTable,
                    file2.charTable,
                    file1.getAtomContainer(),
                    file2.getAtomContainer());
            mcs.search_cliques();
            System.out.println("mcs.final_MAPPINGS " + mcs.final_MAPPINGS.size());
            mappings = Collections.synchronizedList(mcs.final_MAPPINGS);

        } else {
            this.flagExchange = true;
            MoleculeHandler file2 = new MoleculeHandler(source);
            MoleculeHandler file1 = new MoleculeHandler(target);

            MCS mcs = new MCS(
                    file1.getAtomNumber(),
                    file2.getAtomNumber(),
                    file1.getStartHatom_num(),
                    file2.getStartHatom_num(),
                    file1.getBondNumber(),
                    file2.getBondNumber(),
                    file1.getAtomString(),
                    file2.getAtomString(),
                    file1.intTable,
                    file2.intTable,
                    file1.charTable,
                    file2.charTable,
                    file1.getAtomContainer(),
                    file2.getAtomContainer());
            mcs.search_cliques();
//            System.out.println("mcs.final_MAPPINGS " + mcs.final_MAPPINGS.size());
            mappings = Collections.synchronizedList(mcs.final_MAPPINGS);
        }
        List<Map<Integer, Integer>> solutions = PostFilter.filter(mappings);
//        System.out.println("PostFilter.filter " + solutions.size());
        setAllMapping(solutions);
        setAllAtomMapping();

        return !mappings.isEmpty();
    }

    private synchronized void setAllMapping(List<Map<Integer, Integer>> solutions) {
        try {
            int bestSolSize = 0;
            for (Map<Integer, Integer> solution : solutions) {
//                System.out.println("Number of MCS solution: " + solution.size());
                Map<Integer, Integer> validSolution = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
                if (!flagExchange) {
                    solution.entrySet().stream().forEach((map) -> {
                        validSolution.put(map.getKey(), map.getValue());
                    });
                } else {
                    solution.entrySet().stream().forEach((map) -> {
                        validSolution.put(map.getValue(), map.getKey());
                    });
                }
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
        }

//        System.out.println("Number of MCS solution - : allMCS " + allMCS.size());
    }

    private synchronized void setAllAtomMapping() {
//        System.out.println("setAllAtomMapping");
//        System.out.println("source size " + source.getAtomCount());
//        System.out.println("target size " + target.getAtomCount());
        try {
            allMCS.stream().map((solution) -> {
                AtomAtomMapping atomMapping = new AtomAtomMapping(source, target);
                //                System.out.println("solution " + solution);
                solution.entrySet().stream().forEach((m) -> {
                    int indexI = m.getKey() - 1;
                    int indexJ = m.getValue() - 1;
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
            e.printStackTrace();
        }

//        System.out.println("Number of MCS solution - : allAtomMCS " + allAtomMCS.size());
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
}
