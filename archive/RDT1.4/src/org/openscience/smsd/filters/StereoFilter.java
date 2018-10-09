/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.subGraph.uk>
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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.smsd.AtomAtomMapping;

/**
 * Filter on stereo and bond matches.
 *
 * @author Syed Asad Rahman <asad@ebi.subGraph.uk>
 * 
 */
public final class StereoFilter extends Sotter implements IChemicalFilter<Double> {

    private final List<Double> stereoScore;
    private final ChemicalFilters chemfilter;

    StereoFilter(ChemicalFilters chemfilter) {
        this.chemfilter = chemfilter;
        stereoScore = Collections.synchronizedList(new ArrayList<Double>());
    }

    @Override
    public synchronized Double sortResults(
            Map<Integer, AtomAtomMapping> allStereoAtomMCS,
            Map<Integer, Double> stereoScoreMap) throws CDKException {

        getStereoBondChargeMatch(stereoScoreMap, allStereoAtomMCS);

        Map<Integer, Double> sortedStereoScoreMap = sortMapByValueInDescendingOrder(stereoScoreMap);
        double highestStereoScore
                = sortedStereoScoreMap.isEmpty() ? 0
                : sortedStereoScoreMap.values().iterator().next();
        return highestStereoScore;
    }

    @Override
    public synchronized List<Double> getScores() {
        return Collections.unmodifiableList(stereoScore);
    }

    @Override
    public synchronized void clearScores() {
        stereoScore.clear();
    }

    @Override
    public synchronized void addScore(int counter, Double score) {
        stereoScore.add(counter, score);
    }

    @Override
    public synchronized void fillMap(Map<Integer, Double> stereoScoreMap) {
        int Index = 0;
        for (Double score : stereoScore) {
            stereoScoreMap.put(Index, score);
            Index++;
        }
    }

    private synchronized boolean getStereoBondChargeMatch(Map<Integer, Double> stereoScoreMap,
            Map<Integer, AtomAtomMapping> allStereoAtomMCS) throws CDKException {

        boolean stereoMatchFlag = false;
        for (Integer Key : allStereoAtomMCS.keySet()) {
            try {
                double score = 0.0;
                //            System.out.println("\nStart score " + score);
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
                Logger.getLogger(StereoFilter.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return stereoMatchFlag;
    }

    private synchronized Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2,
            AtomAtomMapping mappings) {

        Map<IBond, IBond> bondbondMappingMap = new HashMap<>();

        for (Map.Entry<IAtom, IAtom> map1 : mappings.getMappingsByAtoms().entrySet()) {
            for (Map.Entry<IAtom, IAtom> map2 : mappings.getMappingsByAtoms().entrySet()) {
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

    private synchronized double getAtomScore(double scoreGlobal, AtomAtomMapping atomMapMCS, IAtomContainer reactant,
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

    private synchronized double getBondScore(double scoreGlobal, Map<IBond, IBond> bondMaps) {
        double score = scoreGlobal;
        for (Map.Entry<IBond, IBond> matchedBonds : bondMaps.entrySet()) {

            IBond RBond = matchedBonds.getKey();
            IBond PBond = matchedBonds.getValue();

            score += getBondTypeMatches(RBond, PBond);
        }
        return score;
    }

    private synchronized double getBondTypeMatches(IBond queryBond, IBond targetBond) {
        double score = 0;

        if (targetBond instanceof IQueryBond && queryBond instanceof IBond) {
            IQueryBond bond = (IQueryBond) targetBond;
            IQueryAtom atom1 = (IQueryAtom) (targetBond.getAtom(0));
            IQueryAtom atom2 = (IQueryAtom) (targetBond.getAtom(1));
            if (bond.matches(queryBond)) {
                // ok, bonds match
                if (atom1.matches(queryBond.getAtom(0)) && atom2.matches(queryBond.getAtom(1))
                        || atom1.matches(queryBond.getAtom(1)) && atom2.matches(queryBond.getAtom(0))) {
                    // ok, queryAtom match in either order
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
                // ok, bonds match
                if (atom1.matches(targetBond.getAtom(0)) && atom2.matches(targetBond.getAtom(1))
                        || atom1.matches(targetBond.getAtom(1)) && atom2.matches(targetBond.getAtom(0))) {
                    // ok, queryAtom match in either order
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
            if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                    && (reactantBondType == productBondType)) {
                score += 8;
            } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
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
     *
     * @param bond
     * @return
     */
    public synchronized static int convertBondStereo(IBond bond) {
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
     *
     * @param bond
     * @return
     */
    public synchronized static int convertBondOrder(IBond bond) {
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

    private synchronized double getRingMatchScore(List<IAtomContainer> list) throws CloneNotSupportedException {
        double lScore = 0;
        IAtomContainer listMap = list.get(0).clone();
        IAtomContainer subGraph = list.get(1).clone();
        try {
            Cycles cycles = Cycles.all(subGraph);
            lScore = getRingMatch(cycles.toRingSet(), listMap);
        } catch (Intractable ex) {
            Logger.getLogger(StereoFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return lScore;
    }

    private synchronized double getRingMatch(IRingSet rings, IAtomContainer atoms) {
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

    private synchronized List<IAtomContainer> getMappedFragment(IAtomContainer molecule, Collection<IAtom> atomsMCS) throws CloneNotSupportedException {
        IAtomContainer subgraphContainer;

        if (molecule instanceof IAtomContainer) {
            subgraphContainer = molecule.getBuilder().newInstance(IAtomContainer.class, molecule);
        } else {
            return new ArrayList<>(2);
        }
        List<IAtom> list = new ArrayList<>(atomsMCS.size());
        for (IAtom atom : atomsMCS) {
            int post = molecule.indexOf(atom);
            list.add(subgraphContainer.getAtom(post));
        }

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
