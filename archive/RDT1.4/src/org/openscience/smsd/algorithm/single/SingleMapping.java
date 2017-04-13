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
package org.openscience.smsd.algorithm.single;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.tools.BondEnergies;

/**
 * This class handles single atom mapping. Either query and/or target molecule
 * with single atom is mapped by this class.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class SingleMapping {

    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private final Map<Integer, Double> connectedBondOrder;

    /**
     * Default
     */
    public SingleMapping() {
        connectedBondOrder = new TreeMap<>();
    }

    /**
     * Returns single mapping solutions.
     *
     * @param source
     * @param target
     * @return Mappings
     * @throws CDKException
     */
    protected synchronized List<Map<IAtom, IAtom>> getOverLaps(
            IAtomContainer source,
            IAtomContainer target) throws CDKException {
        List<Map<IAtom, IAtom>> mappings = new ArrayList<>();
        this.source = source;
        this.target = target;

        if (source.getAtomCount() == 1
                || (source.getAtomCount() > 0 && source.getBondCount() == 0)) {
            setSourceSingleAtomMap(mappings);
        }
        if (target.getAtomCount() == 1
                || (target.getAtomCount() > 0 && target.getBondCount() == 0)) {
            setTargetSingleAtomMap(mappings);
        }
        return postFilter(mappings);
    }

    /**
     * Returns single mapping solutions.
     *
     * @param source
     * @param target
     * @return Mappings
     * @throws CDKException
     */
    protected synchronized List<Map<IAtom, IAtom>> getOverLaps(
            IQueryAtomContainer source, IAtomContainer target) throws CDKException {
        List<Map<IAtom, IAtom>> mappings = new ArrayList<>();
        this.source = source;
        this.target = target;

        if (source.getAtomCount() == 1
                || (source.getAtomCount() > 0 && source.getBondCount() == 0)) {
            setSourceSingleAtomMap(mappings);
        }
        if (target.getAtomCount() == 1
                || (target.getAtomCount() > 0 && target.getBondCount() == 0)) {
            setTargetSingleAtomMap(mappings);
        }
        return postFilter(mappings);
    }

    private synchronized void setSourceSingleAtomMap(List<Map<IAtom, IAtom>> mappings) throws CDKException {
        int counter = 0;
        BondEnergies be = BondEnergies.getInstance();
        for (IAtom sourceAtom : source.atoms()) {
            for (IAtom targetAtom : target.atoms()) {
                Map<IAtom, IAtom> mapAtoms = new HashMap<>();
                if (sourceAtom instanceof IQueryAtom) {
                    if (((IQueryAtom) sourceAtom).matches(targetAtom)) {
                        mapAtoms.put(sourceAtom, targetAtom);
                        List<IBond> Bonds = target.getConnectedBondsList(targetAtom);

                        double totalOrder = 0;
                        for (IBond bond : Bonds) {
                            Order bondOrder = bond.getOrder();
                            if (bondOrder == null) {
                                continue;
                            }
                            totalOrder += bondOrder.numeric() + be.getEnergies(bond);
                        }

                        if (!Objects.equals(targetAtom.getFormalCharge(), sourceAtom.getFormalCharge())) {
                            totalOrder += 0.5;
                        }

                        connectedBondOrder.put(counter, totalOrder);
                        mappings.add(counter++, mapAtoms);
                    }
                } else if (sourceAtom.getSymbol().equalsIgnoreCase(targetAtom.getSymbol())) {

                    mapAtoms.put(sourceAtom, targetAtom);
                    List<IBond> Bonds = target.getConnectedBondsList(targetAtom);

                    double totalOrder = 0;
                    for (IBond bond : Bonds) {
                        Order bondOrder = bond.getOrder();
                        if (bondOrder == null) {
                            continue;
                        }
                        totalOrder += bondOrder.numeric() + be.getEnergies(bond);
                    }

                    if (!Objects.equals(targetAtom.getFormalCharge(), sourceAtom.getFormalCharge())) {
                        totalOrder += 0.5;
                    }

                    connectedBondOrder.put(counter, totalOrder);
                    mappings.add(counter, mapAtoms);
                    counter++;
                }
            }
        }
    }

    private synchronized void setTargetSingleAtomMap(List<Map<IAtom, IAtom>> mappings) throws CDKException {
        int counter = 0;
        BondEnergies be = BondEnergies.getInstance();
        for (IAtom targetAtom : target.atoms()) {
            for (IAtom sourceAtoms : source.atoms()) {
                Map<IAtom, IAtom> mapAtoms = new HashMap<>();
                if (targetAtom instanceof IQueryAtom) {
                    if (((IQueryAtom) targetAtom).matches(sourceAtoms)) {
                        if (targetAtom.getSymbol().equalsIgnoreCase(sourceAtoms.getSymbol())) {
                            mapAtoms.put(sourceAtoms, targetAtom);
                            List<IBond> Bonds = source.getConnectedBondsList(sourceAtoms);

                            double totalOrder = 0;
                            for (IBond bond : Bonds) {
                                Order bondOrder = bond.getOrder();
                                if (bondOrder == null) {
                                    continue;
                                }
                                totalOrder += bondOrder.numeric() + be.getEnergies(bond);
                            }
                            if (!Objects.equals(sourceAtoms.getFormalCharge(), targetAtom.getFormalCharge())) {
                                totalOrder += 0.5;
                            }
                            connectedBondOrder.put(counter, totalOrder);
                            mappings.add(counter, mapAtoms);
                            counter++;
                        }
                    }
                } else if (targetAtom.getSymbol().equalsIgnoreCase(sourceAtoms.getSymbol())) {
                    mapAtoms.put(sourceAtoms, targetAtom);
                    List<IBond> Bonds = source.getConnectedBondsList(sourceAtoms);

                    double totalOrder = 0;
                    for (IBond bond : Bonds) {
                        Order bondOrder = bond.getOrder();
                        if (bondOrder == null) {
                            continue;
                        }
                        totalOrder += bondOrder.numeric() + be.getEnergies(bond);
                    }
                    if (!Objects.equals(sourceAtoms.getFormalCharge(), targetAtom.getFormalCharge())) {
                        totalOrder += 0.5;
                    }
                    connectedBondOrder.put(counter, totalOrder);
                    mappings.add(counter, mapAtoms);
                    counter++;
                }
            }
        }
    }

    private synchronized List<Map<IAtom, IAtom>> postFilter(List<Map<IAtom, IAtom>> mappings) {
        List<Map<IAtom, IAtom>> sortedMap = new ArrayList<>();

        if (mappings.isEmpty()) {
            return sortedMap;
        }

        Map<Integer, Double> sortedMapByValue = sortByValue(connectedBondOrder);
        for (Integer key : sortedMapByValue.keySet()) {
            Map<IAtom, IAtom> mapToBeMoved = mappings.get(key);
            sortedMap.add(mapToBeMoved);
        }
        return sortedMap;
    }

    private <K, V> Map<K, V> sortByValue(Map<K, V> map) {
        List list = new LinkedList(map.entrySet());
        Collections.sort(list, new Comparator() {
            @Override
            public int compare(Object object1, Object object2) {
                return ((Comparable) ((Map.Entry<K, V>) (object1)).getValue()).compareTo(
                        ((Map.Entry<K, V>) (object2)).getValue());
            }
        });
        Map<K, V> result = new LinkedHashMap<>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Map.Entry<K, V> entry = (Map.Entry<K, V>) it.next();
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }
}
