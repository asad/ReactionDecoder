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
package org.openscience.smsd;

import java.math.BigDecimal;
import java.util.*;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.filters.ChemicalFilters;
import org.openscience.smsd.interfaces.IAtomMapping;

/**
 *
 *  java1.8+
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
public class BaseMapping extends ChemicalFilters implements IAtomMapping {

    private final boolean matchBonds;
    private final boolean matchRings;
    private final boolean matchAtomType;
    private boolean subgraph;
    private List<Double> stereoScoreList;
    private List<Integer> fragmentSizeList;
    private List<Double> bondEnergiesList;
    private final static ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(BaseMapping.class);

    /**
     *
     * @param matchBonds
     * @param matchRings
     * @param matchAtomType
     * @param mol1
     * @param mol2
     */
    public BaseMapping(IAtomContainer mol1, IAtomContainer mol2, boolean matchBonds,
            boolean matchRings, boolean matchAtomType) {
        super(mol1, mol2);
        this.matchBonds = matchBonds;
        this.matchRings = matchRings;
        this.matchAtomType = matchAtomType;
    }

    /**
     * @param mol1
     * @param mol2
     */
    public BaseMapping(IQueryAtomContainer mol1, IAtomContainer mol2) {
        super(mol1, mol2);
        this.matchBonds = true;
        this.matchRings = true;
        this.matchAtomType = true;

    }

    @Override
    public synchronized void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (getMappingCount() > 0) {

            if (energyFilter) {
                try {
                    sortResultsByEnergies();
                    this.bondEnergiesList = getSortedEnergy();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }

            if (fragmentFilter) {
                sortResultsByFragments();
                this.fragmentSizeList = getSortedFragment();
            }

            if (stereoFilter) {
                try {
                    sortResultsByStereoAndBondMatch();
                    this.stereoScoreList = getStereoMatches();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }
        }
    }

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        return (fragmentSizeList != null && !fragmentSizeList.isEmpty())
                ? fragmentSizeList.get(Key) : null;
    }

    @Override
    public synchronized Integer getStereoScore(int Key) {
        return (stereoScoreList != null && !stereoScoreList.isEmpty()) ? stereoScoreList.get(Key).intValue() : null;
    }

    @Override
    public synchronized Double getEnergyScore(int Key) {
        return (bondEnergiesList != null && !bondEnergiesList.isEmpty()) ? bondEnergiesList.get(Key) : null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized double getTanimotoSimilarity() {
        int decimalPlaces = 4;
        double rAtomCount;
        double pAtomCount;
        double tanimotoAtom = 0.0;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                rAtomCount = (double) this.getMCSList().iterator().next().getQuery().getAtomCount();
                pAtomCount = (double) this.getMCSList().iterator().next().getTarget().getAtomCount();

                double matchCount = (double) firstAtomMCS.getCount();
                tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
                BigDecimal tan = new BigDecimal(tanimotoAtom);
                tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                tanimotoAtom = tan.doubleValue();
            }
        }
        return tanimotoAtom;
    }

    /**
     * {@inheritDoc}
     *
     */
    @Override
    public synchronized boolean isStereoMisMatch() {
        boolean flag = false;
        IAtomContainer reactant = getQuery();
        IAtomContainer product = getTarget();
        int stereoMisMatchScore = 0;
        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();
            for (IAtom indexI : firstAtomMCS.getMappingsByAtoms().keySet()) {
                IAtom indexJ = firstAtomMCS.getMappingsByAtoms().get(indexI);
                for (IAtom indexIPlus : firstAtomMCS.getMappingsByAtoms().keySet()) {
                    IAtom indexJPlus = firstAtomMCS.getMappingsByAtoms().get(indexIPlus);
                    if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                        IAtom sourceAtom1 = indexI;
                        IAtom sourceAtom2 = indexIPlus;

                        IBond rBond = reactant.getBond(sourceAtom1, sourceAtom2);

                        IAtom targetAtom1 = indexJ;
                        IAtom targetAtom2 = indexJPlus;
                        IBond pBond = product.getBond(targetAtom1, targetAtom2);

                        if ((rBond != null && pBond != null)
                                && (rBond.getStereo() != pBond.getStereo())) {
                            stereoMisMatchScore++;
                        }
                    }
                }
            }
        }
        if (stereoMisMatchScore > 0) {
            flag = true;
        }
        return flag;
    }

    @Override
    public synchronized int getMappingCount() {
        return this.getMCSList().isEmpty() ? 0 : this.getMCSList().size();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized double getEuclideanDistance() {
        int decimalPlaces = 4;
        double sourceAtomCount;
        double targetAtomCount;
        double euclidean = -1.;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = getMCSList().iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                sourceAtomCount = (double) this.getMCSList().iterator()
                        .next().getQuery().getAtomCount();
                targetAtomCount = (double) this.getMCSList().iterator()
                        .next().getTarget().getAtomCount();

                double common = (double) firstAtomMCS.getCount();
                euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
                BigDecimal dist = new BigDecimal(euclidean);
                dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                euclidean = dist.doubleValue();
            }
        }
        return euclidean;
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(new ArrayList<>(getMCSList()));
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        return getMCSList().isEmpty() ? new AtomAtomMapping(getQuery(), getTarget())
                : getMCSList().iterator().next();
    }

    /**
     * Returns true if bond are to be matched.
     *
     * @return true if bond are to be matched
     */
    protected synchronized boolean isMatchBonds() {
        return matchBonds;
    }

    /**
     * Returns true if rings are to be matched.
     *
     * @return true if rings are to be matched
     */
    protected synchronized boolean isMatchRings() {
        return matchRings;
    }

    /**
     * Returns true if Query is a subgraph of the Target.
     *
     * @return true if Query is a subgraph of the Target
     */
    public synchronized boolean isSubgraph() {
        return this.subgraph;
    }

    public synchronized void clearMaps() {
        this.getMCSList().clear();
    }

    /**
     * @return the allBondMCS
     */
    public synchronized List<Map<IBond, IBond>> getAllBondMaps() {
        if (!getMCSList().isEmpty()) {
            return makeBondMapsOfAtomMaps(getQuery(), getTarget(), getMCSList());
        }
        return new ArrayList<>();
    }

    /**
     * @param subgraph the subgraph to set
     */
    public synchronized void setSubgraph(boolean subgraph) {
        this.subgraph = subgraph;
    }

    /**
     * Returns bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mappings mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond maps between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    public synchronized List<Map<IBond, IBond>> makeBondMapsOfAtomMaps(IAtomContainer ac1,
            IAtomContainer ac2, List<AtomAtomMapping> mappings) {
        List<Map<IBond, IBond>> bondMaps = Collections.synchronizedList(new ArrayList<Map<IBond, IBond>>());
        for (AtomAtomMapping mapping : mappings) {
            bondMaps.add(makeBondMapOfAtomMap(ac1, ac2, mapping));
        }
        return bondMaps;
    }

    /**
     *
     * Returns bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     *
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mapping mappings between sourceAtomCount and targetAtomCount
     * molecule atoms
     * @return bond map between sourceAtomCount and targetAtomCount molecules
     * based on the atoms
     */
    private synchronized Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2,
            AtomAtomMapping mapping) {

        Map<IBond, IBond> bondbondMappingMap = Collections.synchronizedMap(new HashMap<IBond, IBond>());

        for (Map.Entry<IAtom, IAtom> map1 : mapping.getMappingsByAtoms().entrySet()) {
            for (Map.Entry<IAtom, IAtom> map2 : mapping.getMappingsByAtoms().entrySet()) {
                if (map1.getKey() != map2.getKey()) {
                    IBond bond1 = ac1.getBond(map1.getKey(), map2.getKey());
                    IBond bond2 = ac2.getBond(map1.getValue(), map2.getValue());
                    if (bond1 != null && bond2 != null && !bondbondMappingMap.containsKey(bond1)) {
                        bondbondMappingMap.put(bond1, bond2);
                    }
                }
            }
        }
//        System.out.println("Mol Map size:" + bondbondMappingMap.size());
        return bondbondMappingMap;
    }

    /**
     * @return the matchAtomType
     */
    public boolean isMatchAtomType() {
        return matchAtomType;
    }
}
