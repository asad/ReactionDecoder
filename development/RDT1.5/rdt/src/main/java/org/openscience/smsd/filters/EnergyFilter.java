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
import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.tools.BondEnergies;

/**
 * Filter based on energies.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * 
 */
public final class EnergyFilter extends Sotter implements IChemicalFilter<Double> {

    public static final Double MAX_ENERGY = Double.MAX_VALUE;
    private final List<Double> bEnergies;
    private final ChemicalFilters chemfilter;

    EnergyFilter(ChemicalFilters chemfilter) {
        this.chemfilter = chemfilter;
        bEnergies = Collections.synchronizedList(new ArrayList<Double>());

    }

    @Override
    public synchronized Double sortResults(
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
    public synchronized List<Double> getScores() {
        return Collections.unmodifiableList(bEnergies);
    }

    @Override
    public synchronized void clearScores() {
        bEnergies.clear();
    }

    @Override
    public synchronized void addScore(int counter, Double value) {
        bEnergies.add(counter, value);
    }

    @Override
    public synchronized void fillMap(Map<Integer, Double> energySelectionMap) {
        int Index = 0;
        for (Double score : bEnergies) {
            energySelectionMap.put(Index, score);
            Index++;
        }
    }

    private synchronized Double getMappedMoleculeEnergies(AtomAtomMapping mcsAtomSolution) throws CDKException {

//        System.out.println("\nSort By Energies");
        double totalBondEnergy = -9999.0;

        IAtomContainer educt = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getQuery());
        IAtomContainer product = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, chemfilter.getTarget());

        for (int i = 0; i < educt.getAtomCount(); i++) {
            educt.getAtom(i).setFlag(999, false);
        }

        for (int i = 0; i < product.getAtomCount(); i++) {
            product.getAtom(i).setFlag(999, false);
        }

        if (mcsAtomSolution != null) {
            Map<IAtom, IAtom> mappingsByAtoms = mcsAtomSolution.getMappingsByAtoms();
            mappingsByAtoms.entrySet().stream().map((mapping) -> {
                mapping.getKey().setFlag(999, true);
                return mapping;
            }).forEach((mapping) -> {
                mapping.getValue().setFlag(999, true);
            });
            totalBondEnergy = getEnergy(educt, product);
        }

        /*
         * Reset the flag
         */
        for (int i = 0; i < educt.getAtomCount(); i++) {
            educt.getAtom(i).setFlag(999, false);
        }

        for (int i = 0; i < product.getAtomCount(); i++) {
            product.getAtom(i).setFlag(999, false);
        }

        return totalBondEnergy;
    }

    private synchronized static double getEnergy(IAtomContainer educt, IAtomContainer product) throws CDKException {
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

    private synchronized static double getBondEnergy(IBond bond, BondEnergies bondEnergy) {
        double energy = 0.0;
        if ((bond.getAtom(0).getFlag(999) == true && bond.getAtom(1).getFlag(999) == false)
                || (bond.getAtom(0).getFlag(999) == false && bond.getAtom(1).getFlag(999) == true)) {
            int val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
            energy = val;
        }
        return energy;
    }
}
