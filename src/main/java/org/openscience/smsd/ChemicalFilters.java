/* Copyright (C) 2009-2020  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
package org.openscience.smsd;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.logging.Level;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * A set of filters applied to the results.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ChemicalFilters {

    private static final ILoggingTool BASE_FILTER_LOGGER
            = LoggingToolFactory.createLoggingTool(ChemicalFilters.class);

    // Fields formerly in BaseFilter
    private final IAtomContainer mol1;
    private final IAtomContainer mol2;

    private final List<AtomAtomMapping> allAtomMCS;
    private final IChemicalFilter<Double> energyFilter;
    private final IChemicalFilter<Integer> fragmentFilter;
    private final IChemicalFilter<Double> stereoFilter;

    public ChemicalFilters(IAtomContainer sourceMol, IAtomContainer targetMol) {
        this.mol1 = sourceMol;
        this.mol2 = targetMol;
        this.allAtomMCS = new ArrayList<>();
        this.stereoFilter = new StereoFilter(this);
        this.fragmentFilter = new FragmentFilter(this);
        this.energyFilter = new EnergyFilter(this);
    }

    public ChemicalFilters(IQueryAtomContainer sourceMol, IAtomContainer targetMol) {
        this.mol1 = sourceMol;
        this.mol2 = targetMol;
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
        } catch (CDKException ex) {
            BASE_FILTER_LOGGER.error(Level.SEVERE, null, ex);
        }
        this.allAtomMCS = new ArrayList<>();
        this.stereoFilter = new StereoFilter(this);
        this.fragmentFilter = new FragmentFilter(this);
        this.energyFilter = new EnergyFilter(this);
    }

    /** @return the source molecule */
    public IAtomContainer getQuery() {
        return mol1;
    }

    /** @return the target molecule */
    public IAtomContainer getTarget() {
        return mol2;
    }

    private void clear(
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
    public void sortResultsByEnergies() throws CDKException {
        Map<Integer, AtomAtomMapping> allEnergyAtomMCS = new TreeMap<>();
        Map<Integer, Double> stereoScoreMap = new TreeMap<>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<>();
        Map<Integer, Double> energySelectionMap = new TreeMap<>();

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
    public void sortResultsByFragments() {
        Map<Integer, AtomAtomMapping> allFragmentAtomMCS = new TreeMap<>();
        Map<Integer, Double> stereoScoreMap = new TreeMap<>();
        Map<Integer, Double> energyScoreMap = new TreeMap<>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<>();

        initializeMaps(allFragmentAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);

        try {
            int minFragmentScore = fragmentFilter.sortResults(allFragmentAtomMCS, fragmentScoreMap);

            boolean flag = false;
            if (minFragmentScore < FragmentFilter.MAX_FRAGMENT_SCORE) {
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
    public void sortResultsByStereoAndBondMatch() throws CDKException {
        Map<Integer, AtomAtomMapping> allStereoAtomMCS = new HashMap<>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<>();
        Map<Integer, Double> energyScoreMap = new TreeMap<>();
        Map<Integer, Double> stereoScoreMap = new HashMap<>();

        initializeMaps(allStereoAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);
        double highestStereoScore = stereoFilter.sortResults(allStereoAtomMCS, stereoScoreMap);

        if (highestStereoScore != 0) {
            boolean flag = false;

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

    /** @return sorted bond breaking energy */
    public List<Double> getSortedEnergy() {
        return Collections.unmodifiableList(energyFilter.getScores());
    }

    /** @return sorted fragment count */
    public List<Integer> getSortedFragment() {
        return Collections.unmodifiableList(fragmentFilter.getScores());
    }

    /** @return sorted stereo matches */
    public List<Double> getStereoMatches() {
        return Collections.unmodifiableList(stereoFilter.getScores());
    }

    private void initializeMaps(
            Map<Integer, AtomAtomMapping> sortedAllAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Integer> fragmentScoreMap,
            Map<Integer, Double> energySelectionMap) {

        int index = 0;
        for (AtomAtomMapping atomsMCS : allAtomMCS) {
            sortedAllAtomMCS.put(index, atomsMCS);
            fragmentScoreMap.put(index, 0);
            energySelectionMap.put(index, 0.0);
            stereoScoreMap.put(index, 0.0);
            index++;
        }

        energyFilter.fillMap(energySelectionMap);
        fragmentFilter.fillMap(fragmentScoreMap);
        stereoFilter.fillMap(stereoScoreMap);
    }

    private void addSolution(int counter, int key,
            Map<Integer, AtomAtomMapping> allFragmentAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Double> energyScoreMap,
            Map<Integer, Integer> fragmentScoreMap) {

        allAtomMCS.add(counter, allFragmentAtomMCS.get(key));
        stereoFilter.addScore(counter, stereoScoreMap.get(key));
        fragmentFilter.addScore(counter, fragmentScoreMap.get(key));
        energyFilter.addScore(counter, energyScoreMap.get(key));
    }

    private void clear() {
        allAtomMCS.clear();
        energyFilter.clearScores();
        fragmentFilter.clearScores();
        stereoFilter.clearScores();
    }

    /** @return the mcsList */
    protected List<AtomAtomMapping> getMCSList() {
        return allAtomMCS;
    }

    // ==================== Inner interface ====================

    /**
     * A filter on SMSD results.
     *
     * @param <T>
     * @author Syed Asad Rahman
     * @author maclean
     */
    public static interface IChemicalFilter<T> {

        public T sortResults(
                Map<Integer, AtomAtomMapping> allAtomMCS,
                Map<Integer, T> selectionMap) throws CDKException;

        public List<T> getScores();

        public void clearScores();

        public void addScore(int counter, T value);

        public void fillMap(Map<Integer, T> map);
    }

    // ==================== Inner interface IAtomMapping ====================

    /**
     * Interface for all MCS/Substructure algorithms.
     *
     * @author Syed Asad Rahman
     */
    public static interface IAtomMapping {

        public abstract void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter);

        public abstract Double getEnergyScore(int Key);

        public abstract Integer getFragmentSize(int Key);

        public abstract Integer getStereoScore(int Key);

        public abstract List<AtomAtomMapping> getAllAtomMapping();

        public abstract AtomAtomMapping getFirstAtomMapping();

        public abstract double getTanimotoSimilarity();

        public abstract double getEuclideanDistance();

        public abstract boolean isStereoMisMatch();

        public abstract int getMappingCount();

        @Override
        public abstract String toString();
    }

    // ==================== Inner abstract class BaseFilter ====================

    /**
     * @author Syed Asad Rahman
     * @author maclean
     */
    public static abstract class BaseFilter {

        private final IAtomContainer mol1;
        private final IAtomContainer mol2;
        private final static ILoggingTool BASE_LOGGER
                = LoggingToolFactory.createLoggingTool(BaseFilter.class);

        public BaseFilter(IAtomContainer sourceMol, IAtomContainer targetMol) {
            this.mol1 = sourceMol;
            this.mol2 = targetMol;
        }

        public BaseFilter(IQueryAtomContainer sourceMol, IAtomContainer targetMol) {
            this.mol1 = sourceMol;
            this.mol2 = targetMol;
            try {
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
            } catch (CDKException ex) {
                BASE_LOGGER.error(Level.SEVERE, null, ex);
            }
        }

        /** @return the mol1 */
        public IAtomContainer getQuery() {
            return mol1;
        }

        /** @return the mol2 */
        public IAtomContainer getTarget() {
            return mol2;
        }
    }

    // ==================== Inner class Sotter ====================

    /**
     * @author Syed Asad Rahman
     * @author maclean
     */
    public static class Sotter {

        public static Map<Integer, Double> sortMapByValueInAscendingOrder(Map<Integer, Double> map) {
            List<Map.Entry<Integer, Double>> list = new LinkedList<>(map.entrySet());
            Collections.sort(list, (Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1)));
            Map<Integer, Double> result = new LinkedHashMap<>();
            list.stream().forEach((entry) -> {
                result.put(entry.getKey(), entry.getValue());
            });
            return result;
        }

        public static Map<Integer, Double> sortMapByValueInDescendingOrder(Map<Integer, Double> map) {
            List<Map.Entry<Integer, Double>> list = new LinkedList<>(map.entrySet());
            Collections.sort(list, (Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0
                    : (entry.getValue() < entry1.getValue() ? 1 : -1)));
            Map<Integer, Double> result = new LinkedHashMap<>();
            list.stream().forEach((entry) -> {
                result.put(entry.getKey(), entry.getValue());
            });
            return result;
        }
    }

    // ==================== Inner class EnergyFilter ====================

    /**
     * Filter based on energies.
     * @author Syed Asad Rahman
     */
    static final class EnergyFilter extends Sotter implements IChemicalFilter<Double> {

        public static final Double MAX_ENERGY = Double.MAX_VALUE;
        private final List<Double> bEnergies;
        private final ChemicalFilters chemfilter;

        EnergyFilter(ChemicalFilters chemfilter) {
            this.chemfilter = chemfilter;
            bEnergies = new ArrayList<>();
        }

        @Override
        public Double sortResults(
                Map<Integer, AtomAtomMapping> allAtomEnergyMCS,
                Map<Integer, Double> energySelectionMap) throws CDKException {
            for (Integer Key : allAtomEnergyMCS.keySet()) {
                AtomAtomMapping mcsAtom = allAtomEnergyMCS.get(Key);
                Double energies = getMappedMoleculeEnergies(mcsAtom);
                energySelectionMap.put(Key, energies);
            }

            energySelectionMap = sortMapByValueInAscendingOrder(energySelectionMap);

            double lowestEnergyScore = MAX_ENERGY;
            for (Integer key : energySelectionMap.keySet()) {
                lowestEnergyScore = energySelectionMap.get(key);
                break;
            }
            return lowestEnergyScore;
        }

        @Override
        public List<Double> getScores() {
            return Collections.unmodifiableList(bEnergies);
        }

        @Override
        public void clearScores() {
            bEnergies.clear();
        }

        @Override
        public void addScore(int counter, Double value) {
            bEnergies.add(counter, value);
        }

        @Override
        public void fillMap(Map<Integer, Double> energySelectionMap) {
            int Index = 0;
            for (Double score : bEnergies) {
                energySelectionMap.put(Index, score);
                Index++;
            }
        }

        private Double getMappedMoleculeEnergies(AtomAtomMapping mcsAtomSolution) throws CDKException {
            double totalBondEnergy = -9999.0;

            IAtomContainer educt = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getQuery());
            IAtomContainer product = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getTarget());

            for (int i = 0; i < educt.getAtomCount(); i++) {
                educt.getAtom(i).setProperty("Energy", false);
            }

            for (int i = 0; i < product.getAtomCount(); i++) {
                product.getAtom(i).setProperty("Energy", false);
            }

            if (mcsAtomSolution != null) {
                Map<IAtom, IAtom> mappingsByAtoms = mcsAtomSolution.getMappingsByAtoms();
                mappingsByAtoms.entrySet().stream().map((mapping) -> {
                    mapping.getKey().setProperty("Energy", true);
                    return mapping;
                }).forEach((mapping) -> {
                    mapping.getValue().setProperty("Energy", true);
                });
                totalBondEnergy = getEnergy(educt, product);
            }

            for (int i = 0; i < educt.getAtomCount(); i++) {
                educt.getAtom(i).setProperty("Energy", false);
            }

            for (int i = 0; i < product.getAtomCount(); i++) {
                product.getAtom(i).setProperty("Energy", false);
            }

            return totalBondEnergy;
        }

        private static double getEnergy(IAtomContainer educt, IAtomContainer product) throws CDKException {
            Double eEnergy = 0.0;
            BondEnergies bondEnergy = BondEnergies.getInstance();
            for (int i = 0; i < educt.getBondCount(); i++) {
                IBond bond = educt.getBond(i);
                eEnergy += getBondEnergy(bond, bondEnergy);
            }
            Double pEnergy = 0.0;
            for (int j = 0; j < product.getBondCount(); j++) {
                IBond bond = product.getBond(j);
                pEnergy += getBondEnergy(bond, bondEnergy);
            }
            return (eEnergy + pEnergy);
        }

        private static double getBondEnergy(IBond bond, BondEnergies bondEnergy) {
            double energy = 0.0;
            if ((bond.getAtom(0).getProperty("Energy").equals(true) && bond.getAtom(1).getProperty("Energy").equals(false))
                    || (bond.getAtom(0).getProperty("Energy").equals(false) && bond.getAtom(1).getProperty("Energy").equals(true))) {
                int val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
                energy = val;
            }
            return energy;
        }
    }

    // ==================== Inner class FragmentFilter ====================

    /**
     * Filter the results based on fragment size.
     * @author Syed Asad Rahman
     */
    static final class FragmentFilter extends Sotter implements IChemicalFilter<Integer> {

        static final int MAX_FRAGMENT_SCORE = 9999;

        private final List<Integer> fragmentSize;
        private final ChemicalFilters chemfilter;

        FragmentFilter(ChemicalFilters chemfilter) {
            this.chemfilter = chemfilter;
            fragmentSize = new ArrayList<>();
        }

        @Override
        public Integer sortResults(
                Map<Integer, AtomAtomMapping> allFragmentAtomMCS,
                Map<Integer, Integer> fragmentScoreMap) throws CDKException {

            int _minFragmentScore = MAX_FRAGMENT_SCORE;
            for (Integer key : allFragmentAtomMCS.keySet()) {
                AtomAtomMapping mcsAtom = allFragmentAtomMCS.get(key);
                int fragmentCount = getMappedMoleculeFragmentSize(mcsAtom);
                fragmentScoreMap.put(key, fragmentCount);
                if (_minFragmentScore > fragmentCount) {
                    _minFragmentScore = fragmentCount;
                }
            }

            return _minFragmentScore;
        }

        @Override
        public List<Integer> getScores() {
            return Collections.unmodifiableList(fragmentSize);
        }

        @Override
        public void clearScores() {
            fragmentSize.clear();
        }

        @Override
        public void addScore(int counter, Integer value) {
            fragmentSize.add(counter, value);
        }

        @Override
        public void fillMap(Map<Integer, Integer> fragmentScoreMap) {
            int Index = 0;
            for (Integer score : fragmentSize) {
                fragmentScoreMap.put(Index, score);
                Index++;
            }
        }

        private int getMappedMoleculeFragmentSize(AtomAtomMapping mcsAtomSolution) {
            IAtomContainer Educt = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getQuery());
            IAtomContainer product = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getTarget());

            if (mcsAtomSolution != null) {
                mcsAtomSolution.getMappingsByAtoms().entrySet().stream().forEach((map) -> {
                    IAtom atomE = map.getKey();
                    IAtom atomP = map.getValue();
                    Educt.removeAtom(atomE);
                    product.removeAtom(atomP);
                });
            }
            return getFragmentCount(Educt) + getFragmentCount(product);
        }

        private int getFragmentCount(IAtomContainer molecule) {
            boolean fragmentFlag = true;
            IAtomContainerSet fragmentMolSet = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
            int countFrag = 0;
            if (molecule.getAtomCount() > 0) {
                fragmentFlag = ConnectivityChecker.isConnected(molecule);
                if (!fragmentFlag) {
                    fragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(molecule));
                } else {
                    fragmentMolSet.addAtomContainer(molecule);
                }
                countFrag = fragmentMolSet.getAtomContainerCount();
            }
            return countFrag;
        }
    }

    // ==================== Inner class StereoFilter ====================

    /**
     * Filter on stereo and bond matches.
     * @author Syed Asad Rahman
     */
    static final class StereoFilter extends Sotter implements IChemicalFilter<Double> {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(StereoFilter.class);

        private final List<Double> stereoScore;
        private final ChemicalFilters chemfilter;

        StereoFilter(ChemicalFilters chemfilter) {
            this.chemfilter = chemfilter;
            stereoScore = new ArrayList<>();
        }

        @Override
        public Double sortResults(
                Map<Integer, AtomAtomMapping> allStereoAtomMCS,
                Map<Integer, Double> stereoScoreMap) throws CDKException {

            getStereoBondChargeMatch(stereoScoreMap, allStereoAtomMCS);

            Map<Integer, Double> sortedStereoScoreMap = sortMapByValueInDescendingOrder(stereoScoreMap);
            double highestStereoScore;
            highestStereoScore = sortedStereoScoreMap.isEmpty() ? 0
                    : sortedStereoScoreMap.values().iterator().next();
            return highestStereoScore;
        }

        @Override
        public List<Double> getScores() {
            return Collections.unmodifiableList(stereoScore);
        }

        @Override
        public void clearScores() {
            stereoScore.clear();
        }

        @Override
        public void addScore(int counter, Double score) {
            stereoScore.add(counter, score);
        }

        @Override
        public void fillMap(Map<Integer, Double> stereoScoreMap) {
            int Index = 0;
            for (Double score : stereoScore) {
                stereoScoreMap.put(Index, score);
                Index++;
            }
        }

        private boolean getStereoBondChargeMatch(Map<Integer, Double> stereoScoreMap,
                Map<Integer, AtomAtomMapping> allStereoAtomMCS) throws CDKException {

            boolean stereoMatchFlag = false;
            for (Integer Key : allStereoAtomMCS.keySet()) {
                try {
                    double score = 0.0;

                    AtomAtomMapping atomMapMCS = allStereoAtomMCS.get(Key);
                    double atomScore = getAtomScore(score, atomMapMCS, chemfilter.getQuery(), chemfilter.getTarget());
                    Map<IBond, IBond> bondMaps = makeBondMapsOfAtomMaps(chemfilter.getQuery(), chemfilter.getTarget(), atomMapMCS);
                    double ringScore = 0.0;
                    if (chemfilter.getQuery().getBondCount() > 1
                            && chemfilter.getTarget().getBondCount() > 1
                            && !(chemfilter.getQuery() instanceof IQueryAtomContainer
                            || chemfilter.getTarget() instanceof IQueryAtomContainer)) {
                        List<IAtomContainer> subgraphRList = getMappedFragment(chemfilter.getQuery(), atomMapMCS.getMappingsByAtoms().keySet());
                        double rscore = getRingMatchScore(subgraphRList);
                        List<IAtomContainer> subgraphPList = getMappedFragment(chemfilter.getTarget(), atomMapMCS.getMappingsByAtoms().values());
                        double pscore = getRingMatchScore(subgraphPList);
                        ringScore = rscore + pscore;
                    }
                    double bondScore = getBondScore(score, bondMaps);

                    score = atomScore + ringScore + bondScore;
                    if (!stereoMatchFlag) {
                        stereoMatchFlag = true;
                    }
                    stereoScoreMap.put(Key, score);
                } catch (CloneNotSupportedException ex) {
                    LOGGER.error(Level.SEVERE, null, ex);
                }
            }
            return stereoMatchFlag;
        }

        private Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2,
                AtomAtomMapping mappings) {

            Map<IBond, IBond> bondbondMappingMap = new HashMap<>();

            mappings.getMappingsByAtoms().entrySet().stream().forEach((Map.Entry<IAtom, IAtom> map1) -> {
                mappings.getMappingsByAtoms().entrySet().stream().filter((map2) -> (map1.getKey() != map2.getKey())).forEach((map2) -> {
                    IBond bond1 = ac1.getBond(map1.getKey(), map2.getKey());
                    IBond bond2 = ac2.getBond(map1.getValue(), map2.getValue());
                    if (bond1 != null && bond2 != null && !bondbondMappingMap.containsKey(bond1)) {
                        bondbondMappingMap.put(bond1, bond2);
                    }
                });
            });

            return bondbondMappingMap;
        }

        private double getAtomScore(double scoreGlobal, AtomAtomMapping atomMapMCS, IAtomContainer reactant,
                IAtomContainer product) {
            double score = scoreGlobal;
            for (Map.Entry<IAtom, IAtom> mappings : atomMapMCS.getMappingsByAtoms().entrySet()) {
                IAtom rAtom = mappings.getKey();
                IAtom pAtom = mappings.getValue();

                int rHCount = 0;
                int pHCount = 0;
                double rBO = reactant.getBondOrderSum(rAtom);
                double pBO = product.getBondOrderSum(pAtom);

                if (rAtom.getImplicitHydrogenCount() != null) {
                    rHCount = rAtom.getImplicitHydrogenCount();
                }
                if (pAtom.getImplicitHydrogenCount() != null) {
                    pHCount = pAtom.getImplicitHydrogenCount();
                }

                int HScore = Math.abs(rHCount - pHCount);
                double BOScore = Math.abs(rBO - pBO);

                if (rHCount != pHCount) {
                    score -= HScore;
                } else {
                    score += HScore;
                }

                if (rBO != pBO) {
                    score -= BOScore;
                } else {
                    score += BOScore;
                }

                if (Objects.equals(rAtom.getFormalCharge(), pAtom.getFormalCharge())) {
                    score += 5.0;
                }
            }
            return score;
        }

        private double getBondScore(double scoreGlobal, Map<IBond, IBond> bondMaps) {
            double score = scoreGlobal;
            for (Map.Entry<IBond, IBond> matchedBonds : bondMaps.entrySet()) {
                IBond RBond = matchedBonds.getKey();
                IBond PBond = matchedBonds.getValue();
                score += getBondTypeMatches(RBond, PBond);
            }
            return score;
        }

        private double getBondTypeMatches(IBond queryBond, IBond targetBond) {
            double score = 0;

            if (targetBond instanceof IQueryBond && queryBond instanceof IBond) {
                IQueryBond bond = (IQueryBond) targetBond;
                IQueryAtom atom1 = (IQueryAtom) (targetBond.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (targetBond.getAtom(1));
                if (bond.matches(queryBond)) {
                    if (atom1.matches(queryBond.getAtom(0)) && atom2.matches(queryBond.getAtom(1))
                            || atom1.matches(queryBond.getAtom(1)) && atom2.matches(queryBond.getAtom(0))) {
                        score += 4;
                    }
                } else {
                    score -= 4;
                }
            } else if (queryBond instanceof IQueryBond && targetBond instanceof IBond) {
                IQueryBond bond = (IQueryBond) queryBond;
                IQueryAtom atom1 = (IQueryAtom) (queryBond.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (queryBond.getAtom(1));
                if (bond.matches(targetBond)) {
                    if (atom1.matches(targetBond.getAtom(0)) && atom2.matches(targetBond.getAtom(1))
                            || atom1.matches(targetBond.getAtom(1)) && atom2.matches(targetBond.getAtom(0))) {
                        score += 4;
                    }
                } else {
                    score -= 4;
                }
            } else {
                int reactantBondType = convertBondOrder(queryBond);
                int productBondType = convertBondOrder(targetBond);
                int rStereo = convertBondStereo(queryBond);
                int pStereo = convertBondStereo(targetBond);
                if ((queryBond.isAromatic() == targetBond.isAromatic())
                        && (reactantBondType == productBondType)) {
                    score += 8;
                } else if (queryBond.isAromatic() && targetBond.isAromatic()) {
                    score += 4;
                }

                if (reactantBondType == productBondType) {
                    score += productBondType;
                } else {
                    score -= 4 * Math.abs(reactantBondType - productBondType);
                }

                if (rStereo != 4 || pStereo != 4 || rStereo != 3 || pStereo != 3) {
                    if (rStereo == pStereo) {
                        score += 1;
                    } else {
                        score -= 1;
                    }
                }
            }
            return score;
        }

        /**
         * Get stereo value as integer
         */
        @SuppressWarnings("deprecation")
        public static int convertBondStereo(IBond bond) {
            int value;
            switch (bond.getStereo()) {
                case UP:
                    value = 1;
                    break;
                case UP_INVERTED:
                    value = 1;
                    break;
                case DOWN:
                    value = 6;
                    break;
                case DOWN_INVERTED:
                    value = 6;
                    break;
                case UP_OR_DOWN:
                    value = 4;
                    break;
                case UP_OR_DOWN_INVERTED:
                    value = 4;
                    break;
                case E_OR_Z:
                    value = 3;
                    break;
                default:
                    value = 0;
            }
            return value;
        }

        /**
         * Get bond order value as integer
         */
        public static int convertBondOrder(IBond bond) {
            int value;
            switch (bond.getOrder()) {
                case QUADRUPLE:
                    value = 4;
                    break;
                case TRIPLE:
                    value = 3;
                    break;
                case DOUBLE:
                    value = 2;
                    break;
                case SINGLE:
                    value = 1;
                    break;
                default:
                    value = 1;
            }
            return value;
        }

        private double getRingMatchScore(List<IAtomContainer> list) throws CloneNotSupportedException {
            double lScore = 0;
            IAtomContainer listMap = list.get(0).clone();
            IAtomContainer subGraph = list.get(1).clone();
            try {
                Cycles cycles = Cycles.all(subGraph);
                lScore = getRingMatch(cycles.toRingSet(), listMap);
            } catch (Intractable ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }
            return lScore;
        }

        private double getRingMatch(IRingSet rings, IAtomContainer atoms) {
            double score = 0.0;
            for (IAtom a : atoms.atoms()) {
                for (IAtomContainer ring : rings.atomContainers()) {
                    if (ring.contains(a)) {
                        score += 10;
                    } else {
                        score -= 10;
                    }
                }
            }
            return score;
        }

        private List<IAtomContainer> getMappedFragment(IAtomContainer molecule, Collection<IAtom> atomsMCS) throws CloneNotSupportedException {
            IAtomContainer subgraphContainer;

            if (molecule instanceof IAtomContainer) {
                subgraphContainer = molecule.getBuilder().newInstance(IAtomContainer.class, molecule);
            } else {
                return new ArrayList<>(2);
            }
            List<IAtom> list = new ArrayList<>(atomsMCS.size());
            atomsMCS.stream().map((atom) -> molecule.indexOf(atom)).forEach((post) -> {
                list.add(subgraphContainer.getAtom(post));
            });

            IAtomContainer rlist = new AtomContainer();
            for (IAtom atoms : subgraphContainer.atoms()) {
                if (!list.contains(atoms)) {
                    rlist.addAtom(atoms);
                }
            }

            for (IAtom atoms : rlist.atoms()) {
                subgraphContainer.removeAtom(atoms);
            }
            List<IAtomContainer> l = new ArrayList<>(2);
            l.add(rlist);
            l.add(subgraphContainer);
            return l;
        }
    }
}
