/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.tools;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import uk.ac.ebi.reactionblast.mechanism.helper.BondChange;

/**
 * @refer for valency http://en.wikipedia.org/wiki/Periodic_table_(valence)
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ValencyCalculator {

    private static Map<String, Integer> valencElectronMap = null;
    private static boolean isInitialized = false;
    private static final Logger LOG = Logger.getLogger(ValencyCalculator.class.getName());

    private static void initialize() {
        if (isInitialized) {
            return;
        }
        valencElectronMap = Collections.synchronizedSortedMap(new TreeMap<String, Integer>());
        for (int i = 1; i < PeriodicTable.getElementCount(); i++) {
            String symbol = PeriodicTable.getSymbol(i);
            if (PeriodicTable.getGroup(symbol) != null
                    && (PeriodicTable.getGroup(symbol) < 3 || PeriodicTable.getGroup(symbol) > 12)) {

                switch (PeriodicTable.getGroup(symbol)) {
                    case 1:
                        valencElectronMap.put(symbol, 1);
                        break;
                    case 2:
                        valencElectronMap.put(symbol, 2);
                        break;
                    case 13:
                        valencElectronMap.put(symbol, 3);
                        break;
                    case 14:
                        valencElectronMap.put(symbol, 4);
                        break;
                    case 15:
                        valencElectronMap.put(symbol, 5);
                        break;
                    case 16:
                        valencElectronMap.put(symbol, 6);
                        break;
                    case 17:
                        valencElectronMap.put(symbol, 7);
                        break;
                    case 18:
                        valencElectronMap.put(symbol, 8);
                        break;
                    default:
                        valencElectronMap.put(symbol, 0);
                        break;
                }
            } else {
                valencElectronMap.put(symbol, 99);
            }
//            System.out.println("Atom " + symbol + " grp: " + PeriodicTable.getGroup(symbol) + " P " + PeriodicTable.getPeriod(symbol));
        }
        /* 
         * Metal
         */
        valencElectronMap.put("Sc", 3);
        valencElectronMap.put("Ti", 4);
        valencElectronMap.put("V", 5);
        valencElectronMap.put("Cr", 6);
        valencElectronMap.put("Mn", 4);
        valencElectronMap.put("Ni", 2);
        valencElectronMap.put("Cu", 2);
        valencElectronMap.put("Zn", 2);
        valencElectronMap.put("Fe", 3);
        valencElectronMap.put("Co", 3);
        /* 
         * Generics
         */
        valencElectronMap.put("*", 1);
        valencElectronMap.put("R", 1);
        valencElectronMap.put("A", 1);
        valencElectronMap.put("X", 8);
        valencElectronMap.put("PsH", 1);
//        System.out.println("valencElectronMap Size " + valencElectronMap.size());
        isInitialized = true;
    }

    /**
     * This method calculates the valence of an atom.
     *
     * @param atom The IAtom for which the DescriptorValue is requested
     * @return atomValence The valency for the given atom
     * @throws CDKException
     */
    public static Integer getValenceElectron(IAtom atom) throws CDKException {
        initialize();
        Integer atomValence;
        String symbol = atom.getSymbol();
        if (valencElectronMap.containsKey(symbol)) {
            atomValence = valencElectronMap.get(symbol);
        } else {
            Logger.getLogger(ValencyCalculator.class.getName()).log(Level.WARNING, "Element {0} not found. Valence assigned 99.", symbol);
            atomValence = 99;
        }
        return atomValence;
    }

    /**
     *
     * @param m
     * @param atom
     * @param skipHydrogen
     * @return
     * @throws CDKException
     */
    public static Integer getFreeValenceElectrons(IAtomContainer m, IAtom atom, boolean skipHydrogen) throws CDKException {
        initialize();
        Integer totalConnectedBondOrder = 0;
        List<IAtom> connectedAtoms = m.getConnectedAtomsList(atom);
        int counterH = 0;
        for (IAtom connAtom : connectedAtoms) {
            if (skipHydrogen && connAtom.getSymbol().equalsIgnoreCase("H")) {
                counterH++;
            }
            IBond bond = m.getBond(atom, connAtom);
            totalConnectedBondOrder += BondChange.convertBondOrder(bond);
        }
        Integer charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
        return skipHydrogen ? (getValenceElectron(atom) - totalConnectedBondOrder + counterH - charge)
                : (getValenceElectron(atom) - totalConnectedBondOrder - charge);
    }

    /**
     * @return Elements
     */
    public static String[] getElements() {
        initialize();
        String[] st = new String[valencElectronMap.size()];
        int i = 0;
        for (String s : valencElectronMap.keySet()) {
            st[i++] = s;
        }
        return st;
    }

    /**
     * @return Element Map
     */
    public static Map<String, Integer> getElementMap() {
        initialize();
        return Collections.unmodifiableMap(valencElectronMap);
    }

    public static int getSize() {
        initialize();
        return valencElectronMap.size();
    }

    public static Iterable<String> getKeySet() {
        initialize();
        return valencElectronMap.keySet();
    }

    /**
     *
     * @param key
     * @return the valence
     */
    public static int getValue(String key) {
        initialize();
        return valencElectronMap.get(key);
    }
}
