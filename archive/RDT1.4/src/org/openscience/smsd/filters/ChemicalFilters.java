/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
 */
package org.openscience.smsd.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;

/**
 *
 * A set of filters applied to the results.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * 
 *
 */
public class ChemicalFilters extends BaseFilter {

    private final List<AtomAtomMapping> allAtomMCS;
    private final IChemicalFilter<Double> energyFilter;
    private final IChemicalFilter<Integer> fragmentFilter;
    private final IChemicalFilter<Double> stereoFilter;

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public ChemicalFilters(
            IAtomContainer sourceMol,
            IAtomContainer targetMol) {
        super(sourceMol, targetMol);
        this.allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        this.stereoFilter = new StereoFilter(this);
        this.fragmentFilter = new FragmentFilter(this);
        this.energyFilter = new EnergyFilter(this);
    }

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public ChemicalFilters(
            IQueryAtomContainer sourceMol,
            IAtomContainer targetMol) {
        super(sourceMol, targetMol);
        this.allAtomMCS = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        this.stereoFilter = new StereoFilter(this);
        this.fragmentFilter = new FragmentFilter(this);
        this.energyFilter = new EnergyFilter(this);
    }

    private synchronized void clear(
            Map<Integer, AtomAtomMapping> sortedAllAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Integer> fragmentScoreMap,
            Map<Integer, Double> energySelectionMap) {
        sortedAllAtomMCS.clear();
        stereoScoreMap.clear();
        fragmentScoreMap.clear();
        energySelectionMap.clear();
    }

    /**
     * Sort MCS solution by bond breaking energy.
     *
     * @throws CDKException
     */
    public synchronized void sortResultsByEnergies() throws CDKException {
        Map<Integer, AtomAtomMapping> allEnergyAtomMCS = Collections.synchronizedSortedMap(new TreeMap<Integer, AtomAtomMapping>());
        Map<Integer, Double> stereoScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Double>());
        Map<Integer, Integer> fragmentScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
        Map<Integer, Double> energySelectionMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Double>());

        initializeMaps(allEnergyAtomMCS, stereoScoreMap, fragmentScoreMap, energySelectionMap);
        double lowestEnergyScore = energyFilter.sortResults(allEnergyAtomMCS, energySelectionMap);
        clear();

        int counter = 0;
        for (Map.Entry<Integer, Double> map : energySelectionMap.entrySet()) {
            if (lowestEnergyScore == map.getValue()) {
                addSolution(counter, map.getKey(),
                        allEnergyAtomMCS,
                        stereoScoreMap,
                        energySelectionMap,
                        fragmentScoreMap);
                counter++;
            }
        }

        if (lowestEnergyScore != EnergyFilter.MAX_ENERGY) {
            clear(allEnergyAtomMCS, stereoScoreMap, fragmentScoreMap, energySelectionMap);
        }
    }

    /**
     * Sort solution by ascending order of the fragment count.
     */
    public synchronized void sortResultsByFragments() {
        Map<Integer, AtomAtomMapping> allFragmentAtomMCS = Collections.synchronizedSortedMap(new TreeMap<Integer, AtomAtomMapping>());
        Map<Integer, Double> stereoScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Double>());
        Map<Integer, Double> energyScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Double>());
        Map<Integer, Integer> fragmentScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());

        initializeMaps(
                allFragmentAtomMCS,
                stereoScoreMap,
                fragmentScoreMap,
                energyScoreMap);

        try {
            int minFragmentScore = fragmentFilter.sortResults(allFragmentAtomMCS, fragmentScoreMap);

            boolean flag = false;
            if (minFragmentScore < 9999) {
                flag = true;
                clear();
            }
            int counter = 0;
            for (Map.Entry<Integer, Integer> map : fragmentScoreMap.entrySet()) {
                if (minFragmentScore == map.getValue()) {
                    addSolution(counter, map.getKey(),
                            allFragmentAtomMCS,
                            stereoScoreMap,
                            energyScoreMap,
                            fragmentScoreMap);
                    counter++;
                }
            }

            if (flag) {
                clear(allFragmentAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);
            }
        } catch (CDKException c) {
            // actually, never thrown, but in the interface
        }

    }

    /**
     * Sort MCS solution by stereo and bond type matches.
     *
     * @throws CDKException
     */
    public synchronized void sortResultsByStereoAndBondMatch() throws CDKException {
        Map<Integer, AtomAtomMapping> allStereoAtomMCS = Collections.synchronizedMap(new HashMap<Integer, AtomAtomMapping>());
        Map<Integer, Integer> fragmentScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
        Map<Integer, Double> energyScoreMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Double>());
        Map<Integer, Double> stereoScoreMap = Collections.synchronizedMap(new HashMap<Integer, Double>());

        initializeMaps(
                allStereoAtomMCS,
                stereoScoreMap,
                fragmentScoreMap,
                energyScoreMap);
        double highestStereoScore = stereoFilter.sortResults(allStereoAtomMCS, stereoScoreMap);

        if (highestStereoScore != 0) {
            boolean flag = false;

            //Higher Score is mapped preferred over lower
            double secondhigestStereoScore = highestStereoScore;
            for (Integer key : stereoScoreMap.keySet()) {
                if (secondhigestStereoScore < highestStereoScore
                        && stereoScoreMap.get(key) > secondhigestStereoScore) {
                    secondhigestStereoScore = stereoScoreMap.get(key);
                } else if (secondhigestStereoScore == highestStereoScore
                        && stereoScoreMap.get(key) < secondhigestStereoScore) {
                    secondhigestStereoScore = stereoScoreMap.get(key);
                }
            }

            if (!stereoScoreMap.isEmpty()) {
                flag = true;
                clear();
            }

            /*Put back the sorted solutions*/
            int counter = 0;
            for (Integer I : stereoScoreMap.keySet()) {
                if (highestStereoScore == stereoScoreMap.get(I)) {
                    addSolution(counter, I,
                            allStereoAtomMCS,
                            stereoScoreMap,
                            energyScoreMap,
                            fragmentScoreMap);
                    counter++;

                }
            }
            if (flag) {
                clear(allStereoAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);
            }
        }
    }

    /**
     * Return sorted energy in ascending order.
     *
     * @return sorted bond breaking energy
     */
    public synchronized List<Double> getSortedEnergy() {
        return Collections.unmodifiableList(energyFilter.getScores());
    }

    /**
     * Return sorted fragment in ascending order of the size.
     *
     * @return sorted fragment count
     */
    public synchronized List<Integer> getSortedFragment() {
        return Collections.unmodifiableList(fragmentFilter.getScores());
    }

    /**
     * Return Stereo matches in descending order.
     *
     * @return sorted stereo matches
     */
    public synchronized List<Double> getStereoMatches() {
        return Collections.unmodifiableList(stereoFilter.getScores());
    }

    private synchronized void initializeMaps(
            Map<Integer, AtomAtomMapping> sortedAllAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Integer> fragmentScoreMap,
            Map<Integer, Double> energySelectionMap) {

        Integer Index = 0;
        for (AtomAtomMapping atomsMCS : allAtomMCS) {
            sortedAllAtomMCS.put(Index, atomsMCS);
            fragmentScoreMap.put(Index, 0);
            energySelectionMap.put(Index, 0.0);
            stereoScoreMap.put(Index, 0.0);
            Index++;
        }

        energyFilter.fillMap(energySelectionMap);
        fragmentFilter.fillMap(fragmentScoreMap);
        stereoFilter.fillMap(stereoScoreMap);

    }

    private synchronized void addSolution(int counter, int key,
            Map<Integer, AtomAtomMapping> allFragmentAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Double> energyScoreMap,
            Map<Integer, Integer> fragmentScoreMap) {

        allAtomMCS.add(counter, allFragmentAtomMCS.get(key));
        stereoFilter.addScore(counter, stereoScoreMap.get(key));
        fragmentFilter.addScore(counter, fragmentScoreMap.get(key));
        energyFilter.addScore(counter, energyScoreMap.get(key));

    }

    private synchronized void clear() {
        allAtomMCS.clear();
        energyFilter.clearScores();
        fragmentFilter.clearScores();
        stereoFilter.clearScores();
    }

    /**
     * @return the mcsList
     */
    protected synchronized List<AtomAtomMapping> getMCSList() {
        return Collections.synchronizedList(allAtomMCS);
    }
}
