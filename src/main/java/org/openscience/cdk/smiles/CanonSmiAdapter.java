/*
 * Copyright (c) 2014. EMBL, European Bioinformatics Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.openscience.cdk.smiles;

import java.io.IOException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import static org.openscience.cdk.graph.GraphUtil.EdgeToBondMap.withSpaceFor;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.ringsearch.RingSearch;
import uk.ac.ebi.beam.Graph;
import static java.lang.Math.abs;
import static org.openscience.cdk.graph.GraphUtil.toAdjList;
import static org.openscience.cdk.graph.invariant.Canon.label;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * @author John May
 */
public class CanonSmiAdapter {

    // convert to Beam excluding stereo (not canonicalised) and aromaticity
    static final CDKToBeam CDK2BEAM = new CDKToBeam();
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CanonSmiAdapter.class);

    /**
     *
     * @param ac
     * @return
     * @throws CDKException
     * @throws IOException
     */
    public static String create(IAtomContainer ac) throws CDKException, IOException {

        EdgeToBondMap bonds = withSpaceFor(ac);
        int[][] graph = toAdjList(ac, bonds);

        long[] labels = label(ac, graph, betterInvariants(ac, graph, bonds));

        Graph g = CDK2BEAM.toBeamGraph(ac);
        g = g.permute(toPermutation(labels));
        g = g.resonate(); // ensure consistent Kekule form 

        return g.toSmiles();
    }

    private static int[] toPermutation(long[] labels) throws CDKException {
        int[] cpy = new int[labels.length];
        for (int i = 0; i < labels.length; i++) {
            cpy[i] = (int) labels[i] - 1;
        }
        return cpy;
    }

    /**
     *
     * @param container
     * @param graph
     * @param bonds
     * @return
     */
    public static long[] betterInvariants(IAtomContainer container, int[][] graph, EdgeToBondMap bonds) {
        long[] labels = new long[graph.length];

        RingSearch ringSearch = new RingSearch(container, graph);

        for (int v = 0; v < graph.length; v++) {
            IAtom atom = container.getAtom(v);

            int deg = graph[v].length;
            int impH = implH(atom);
            int expH = 0;
            int elem = atomicNumber(atom);
            int chg = charge(atom);

            int valence = impH;

            // count non-suppressed (explicit) hydrogens and valence
            for (int w : graph[v]) {
                if (atomicNumber(container.getAtom(w)) == 1) {
                    expH++;
                }
                if (bonds.get(v, w).getOrder() == null) {
                    continue;
                }
                valence += bonds.get(v, w).getOrder().numeric();
            }

            long label = 0; // connectivity (first in)
            label |= deg + impH & 0xf;
            label <<= 4;    // connectivity (heavy) <= 15 (4 bits)
            label |= deg - expH & 0xf;
            label <<= 7;   // atomic number <= 127 (7 bits)
            label |= elem & 0x7f;
            label <<= 4;    // valence <= 15 (4 bits)
            label |= valence & 0xf;
            label <<= 1;   // charge sign == 1 (1 bit)
            label |= chg >> 31 & 0x1;
            label <<= 2;   // charge <= 3 (2 bits)
            label |= abs(chg) & 0x3;
            label <<= 4;   // hydrogen count <= 15 (4 bits)
            label |= impH + expH & 0xf;
            label <<= 1;   // ring membership (1 bit)
            label |= ringSearch.cyclic(v) ? 0 : 1;

            labels[v] = label;
        }
        return labels;
    }

    private static int atomicNumber(IAtom atom) {
        Integer elem = atom.getAtomicNumber();
        if (elem != null) {
            return elem;
        }
        if (atom instanceof IPseudoAtom) {
            return 0;
        }
        throw new NullPointerException("a non-psuedo atom had unset atomic number");
    }

    private static int implH(IAtom atom) {
        Integer h = atom.getImplicitHydrogenCount();
        if (h != null) {
            return h;
        }
        if (atom instanceof IPseudoAtom) {
            return 0;
        }
        throw new NullPointerException("a non-psuedo atom had unset hydrogen count");
    }

    private static int charge(IAtom atom) {
        Integer charge = atom.getFormalCharge();
        if (charge != null) {
            return charge;
        }
        return 0;
    }

    private CanonSmiAdapter() {
    }
}
