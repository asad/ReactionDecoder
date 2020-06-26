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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.helper;

import java.util.ArrayList;
import java.util.Collections;
import static java.util.Collections.sort;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import org.openscience.cdk.CDKConstants;
import static org.openscience.cdk.CDKConstants.UNSET;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.aromatizeCDK;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.aromatizeDayLight;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;

/**
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MoleculeInitializer {

    /**
     * Prepare the molecule for analysis.
     * <p/>
     * We perform ring perception and aromaticity detection and set up the
     * appropriate properties. Right now, this function is called each time we
     * need to do a query and this is inefficient.
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MoleculeInitializer.class);

    /**
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    public synchronized static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        if (atomContainer == null) {
            return;
        }
        try {
            try {
                // figure out which atoms are in aromatic rings:
                percieveAtomTypesAndConfigureAtoms(atomContainer);
                aromatizeCDK(atomContainer);
            } catch (CDKException e) {
                aromatizeDayLight(atomContainer);
            }
        } catch (CDKException e) {
            LOGGER.error(Level.WARNING, "Error in aromaticity dectection. ", atomContainer.getID());
        }

        String SMALLEST_RING_SIZE = "SMALLEST_RING_SIZE";
        if (!(atomContainer instanceof IQueryAtomContainer)) {
            Map<String, Integer> valencesTable = new HashMap<>();
            valencesTable.put("H", 1);
            valencesTable.put("Li", 1);
            valencesTable.put("Be", 2);
            valencesTable.put("B", 3);
            valencesTable.put("C", 4);
            valencesTable.put("N", 5);
            valencesTable.put("O", 6);
            valencesTable.put("F", 7);
            valencesTable.put("Na", 1);
            valencesTable.put("Mg", 2);
            valencesTable.put("Al", 3);
            valencesTable.put("Si", 4);
            valencesTable.put("P", 5);
            valencesTable.put("S", 6);
            valencesTable.put("Cl", 7);
            valencesTable.put("K", 1);
            valencesTable.put("Ca", 2);
            valencesTable.put("Ga", 3);
            valencesTable.put("Ge", 4);
            valencesTable.put("As", 5);
            valencesTable.put("Se", 6);
            valencesTable.put("Br", 7);
            valencesTable.put("Rb", 1);
            valencesTable.put("Sr", 2);
            valencesTable.put("In", 3);
            valencesTable.put("Sn", 4);
            valencesTable.put("Sb", 5);
            valencesTable.put("Te", 6);
            valencesTable.put("I", 7);
            valencesTable.put("Cs", 1);
            valencesTable.put("Ba", 2);
            valencesTable.put("Tl", 3);
            valencesTable.put("Pb", 4);
            valencesTable.put("Bi", 5);
            valencesTable.put("Po", 6);
            valencesTable.put("At", 7);
            valencesTable.put("Fr", 1);
            valencesTable.put("Ra", 2);
            valencesTable.put("Cu", 2);
            valencesTable.put("Mn", 2);
            valencesTable.put("Co", 2);

            // do all ring perception
            IRingSet allRings = null;
            CycleFinder cycleFinder = Cycles.or(Cycles.all(),
                    Cycles.or(Cycles.relevant(),
                            Cycles.essential()));

            Cycles cycles = cycleFinder.find(atomContainer);
            allRings = cycles.toRingSet();
            /*
             * Mark aromatic rings
             */
            RingSetManipulator.markAromaticRings(allRings);

            // sets SSSR information
            //IRingSet sssr = new SSSRFinder(atomContainer).findEssentialRings();
            //New Method
            CycleFinder cf = Cycles.essential();
            cycles = cf.find(atomContainer); // ignore error - essential cycles do not check tractability
            IRingSet sssr = cycles.toRingSet();

            for (IAtom atom : atomContainer.atoms()) {
                if (atom == null) {
                    continue;
                }

                // add a property to each ring atom that will be an array of
                // Integers, indicating what size ring the given atom belongs to
                // Add SSSR ring counts
                if (allRings != null && allRings.contains(atom)) { // it's in a ring
                    atom.setIsInRing(true);
                    atom.setIsAromatic(true);
                    atom.setFlag(CDKConstants.ISINRING, true);
                    atom.setFlag(CDKConstants.ISALIPHATIC, false);
                    // lets find which ring sets it is a part of
                    List<Integer> ringsizes = new ArrayList<>();
                    IRingSet currentRings = allRings.getRings(atom);
                    int min = 0;
                    for (int i = 0; i < currentRings.getAtomContainerCount(); i++) {
                        int size = currentRings.getAtomContainer(i).getAtomCount();
                        if (min > size) {
                            min = size;
                        }
                        ringsizes.add(size);
                    }
                    Collections.sort(ringsizes);
                    atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
                    atom.setProperty(CDKConstants.SMALLEST_RINGS, sssr.getRings(atom));
                    atom.setProperty(SMALLEST_RING_SIZE, min);
                } else {
                    atom.setIsInRing(false);
                    atom.setIsAromatic(false);
                    atom.setFlag(CDKConstants.ISINRING, false);
                    atom.setFlag(CDKConstants.ISALIPHATIC, true);
                    atom.setProperty(SMALLEST_RING_SIZE, 0);
                }

                // determine how many rings bonds each atom is a part of
                int hCount;
                if (Objects.equals(atom.getImplicitHydrogenCount(), CDKConstants.UNSET)) {
                    hCount = 0;
                } else {
                    hCount = atom.getImplicitHydrogenCount();
                }

                List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
                int total = hCount + connectedAtoms.size();
                hCount = connectedAtoms.stream().filter((connectedAtom)
                        -> (connectedAtom.getSymbol().equals("H"))).map((_item) -> 1)
                        .reduce(hCount, Integer::sum);
                atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
                atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);

                if (valencesTable.get(atom.getSymbol()) != null) {
                    int formalCharge = Objects.equals(atom.getFormalCharge(), CDKConstants.UNSET)
                            ? 0 : atom.getFormalCharge();
                    atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
                }
            }

            for (IBond bond : atomContainer.bonds()) {
                if (allRings != null && allRings.getRings(bond).getAtomContainerCount() > 0) {
                    bond.setFlag(CDKConstants.ISINRING, true);
                    bond.setFlag(CDKConstants.ISALIPHATIC, false);

                    bond.setIsInRing(true);
                    bond.setIsAromatic(true);
                } else {
                    bond.setIsInRing(false);
                    bond.setIsAromatic(false);
                }
            }

            for (IAtom atom : atomContainer.atoms()) {
                List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);

                int counter = 0;
                IAtom any;
                for (IAtom connectedAtom : connectedAtoms) {
                    any = connectedAtom;
                    if (any.getFlag(CDKConstants.ISINRING)) {
                        counter++;
                    }
                }
                atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
            }

            ExtAtomContainerManipulator.aromatizeMolecule(atomContainer);
        }

    }

    /*
     * Checks some simple heuristics for whether the subgraph query can
     * realistically be atom subgraph of the supergraph. If, for example, the
     * number of nitrogen atoms in the query is larger than that of the
     * supergraph it cannot be part of it.
     *
     *
     *
     * @param q Query
     * @param t Target
     * @return true if subgraph else false
     */
    public static boolean testIsSubgraphHeuristics(IAtomContainer q, IAtomContainer t) {

        /*
         a={c,c,c,o,n}
         b={c,c,c,p}
       
         expectedMaxGraphmatch=3;
         */
        List<String> atomUniqueCounter1 = new ArrayList<>();
        List<String> atomUniqueCounter2 = new ArrayList<>();

        for (IAtom a : q.atoms()) {
            String hyb = a.getHybridization() == UNSET
                    ? a.getSymbol() : a.getAtomTypeName();
            atomUniqueCounter1.add(hyb);
        }

        for (IAtom b : t.atoms()) {
            String hyb = b.getHybridization() == UNSET
                    ? b.getSymbol() : b.getAtomTypeName();
            atomUniqueCounter2.add(hyb);
        }

        sort(atomUniqueCounter1);
        sort(atomUniqueCounter2);

        if (atomUniqueCounter1.isEmpty()) {
            return false;
        }
        List<String> common = new LinkedList<>(atomUniqueCounter1);
        common.retainAll(atomUniqueCounter2);

        atomUniqueCounter1.clear();
        atomUniqueCounter2.clear();
        return common.size() == atomUniqueCounter1.size() ? true : common.size() == atomUniqueCounter2.size();
    }

    /**
     * Checks some simple heuristics for whether the subgraph query can
     * realistically be atom subgraph of the supergraph. If, for example, the
     * number of nitrogen atoms in the query is larger than that of the
     * supergraph it cannot be part of it.
     *
     * @param ac1 the supergraph to be checked.
     * @param ac2 the subgraph to be tested for. Must not be an
     * IQueryAtomContainer.
     * @param shouldMatchBonds
     * @return true if the subgraph ac1 has atom chance to be atom subgraph of
     * ac2
     */
    public synchronized static boolean testIsSubgraphHeuristics(
            IAtomContainer ac1,
            IAtomContainer ac2,
            AtomMatcher am,
            BondMatcher bm) {

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;

        IBond bond;

        if (false) {
            for (int i = 0; i < ac1.getBondCount(); i++) {
                bond = ac1.getBond(i);
                if (bond == null) {
                    continue;
                }
                if (bond instanceof IQueryBond) {
                    continue;
                }
                if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                    ac1AromaticBondCount++;
                } else if (bond.getOrder() == IBond.Order.SINGLE) {
                    ac1SingleBondCount++;
                } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                    ac1DoubleBondCount++;
                } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                    ac1TripleBondCount++;
                }
            }
            for (int indexI = 0; indexI < ac2.getBondCount(); indexI++) {
                bond = ac2.getBond(indexI);
                if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                    ac2AromaticBondCount++;
                } else if (bond.getOrder() == IBond.Order.SINGLE) {
                    ac2SingleBondCount++;
                } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                    ac2DoubleBondCount++;
                } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                    ac2TripleBondCount++;
                }
            }

            if (ac2SingleBondCount < ac1SingleBondCount) {
                return false;
            }
            if (ac2AromaticBondCount < ac1AromaticBondCount) {
                return false;
            }
            if (ac2DoubleBondCount < ac1DoubleBondCount) {
                return false;
            }
            if (ac2TripleBondCount < ac1TripleBondCount) {
                return false;
            }
        }

        IAtom atom;
        Map<String, Integer> map = new HashMap<>();
        for (int i = 0; i < ac1.getAtomCount(); i++) {
            atom = ac1.getAtom(i);
            if (atom == null) {
                continue;
            }
            if (atom instanceof IQueryAtom) {
                continue;
            }
            if (map.containsKey(atom.getSymbol())) {
                int val = map.get(atom.getSymbol()) + 1;
                map.put(atom.getSymbol(), val);
            } else {
                map.put(atom.getSymbol(), 1);
            }
        }
        for (int i = 0; i < ac2.getAtomCount(); i++) {
            atom = ac2.getAtom(i);
            if (map.containsKey(atom.getSymbol())) {
                int val = map.get(atom.getSymbol()) - 1;
                if (val > 0) {
                    map.put(atom.getSymbol(), val);
                } else {
                    map.remove(atom.getSymbol());
                }
            }
        }
        return map.isEmpty();
    }
}
