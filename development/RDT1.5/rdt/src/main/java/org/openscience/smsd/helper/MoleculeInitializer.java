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
package org.openscience.smsd.helper;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class MoleculeInitializer {

    /**
     * Defines which set of rings to define rings in the target.
     */
    private enum RingSet {

        /**
         * Smallest Set of Smallest Rings (or Minimum Cycle Basis - but not
         * strictly the same). Defines what is typically thought of as a 'ring'
         * however the non-uniqueness leads to ambiguous matching.
         */
        SmallestSetOfSmallestRings {
                    @Override
                    IRingSet ringSet(IAtomContainer m) {
                        return new SSSRFinder(m).findSSSR();
                    }
                },
        /**
         * Intersect of all Minimum Cycle Bases (or SSSR) and thus is a subset.
         * The set is unique but may excludes rings (e.g. from bridged systems).
         */
        EssentialRings {
                    @Override
                    IRingSet ringSet(IAtomContainer m) {
                        return new SSSRFinder(m).findEssentialRings();
                    }
                },
        /**
         * Union of all Minimum Cycle Bases (or SSSR) and thus is a superset.
         * The set is unique but may include more rings then is necessary.
         */
        RelevantRings {
                    @Override
                    IRingSet ringSet(IAtomContainer m) {
                        return new SSSRFinder(m).findRelevantRings();
                    }
                };

        /**
         * Compute a ring set for a molecule.
         *
         * @param m molecule
         * @return the ring set for the molecule
         */
        abstract IRingSet ringSet(IAtomContainer m);
    }

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
    private static final ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(MoleculeInitializer.class);

    /**
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    public synchronized static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
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
            AllRingsFinder arf = new AllRingsFinder();
            IRingSet allRings = null;
            try {
                allRings = arf.findAllRings(atomContainer);
            } catch (CDKException e) {
                Logger.warn(e.toString());
            }

            // sets SSSR information
            IRingSet sssr = new SSSRFinder(atomContainer).findEssentialRings();

            for (IAtom atom : atomContainer.atoms()) {

                // add a property to each ring atom that will be an array of
                // Integers, indicating what size ring the given atom belongs to
                // Add SSSR ring counts
                if (allRings != null && allRings.contains(atom)) { // it's in a ring
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
                for (IAtom connectedAtom : connectedAtoms) {
                    if (connectedAtom.getSymbol().equals("H")) {
                        hCount++;
                    }
                }
                atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
                atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);

                if (valencesTable.get(atom.getSymbol()) != null) {
                    int formalCharge = Objects.equals(atom.getFormalCharge(), CDKConstants.UNSET) ? 0 : atom.getFormalCharge();
                    atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
                }
            }

            for (IBond bond : atomContainer.bonds()) {
                if (allRings != null && allRings.getRings(bond).getAtomContainerCount() > 0) {
                    bond.setFlag(CDKConstants.ISINRING, true);
                    bond.setFlag(CDKConstants.ISALIPHATIC, false);
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
    public synchronized static boolean testIsSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2, boolean shouldMatchBonds) {

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;

        IBond bond;

        if (shouldMatchBonds) {
            for (int i = 0; i < ac1.getBondCount(); i++) {
                bond = ac1.getBond(i);
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
