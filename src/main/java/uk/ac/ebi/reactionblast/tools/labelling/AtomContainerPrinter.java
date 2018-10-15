/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.tools.labelling;

import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class AtomContainerPrinter {

    /**
     *
     * @param atomContainer
     * @return
     */
    public String toString(IAtomContainer atomContainer) {
        StringBuilder sb = new StringBuilder();
        for (IAtom atom : atomContainer.atoms()) {
            sb.append(atom.getSymbol());
        }
        sb.append(" ");
        List<Edge> edges = new ArrayList<>();
        for (IBond bond : atomContainer.bonds()) {
            if (bond.getAtomCount() < 2) {
                edges.add(new Edge(-1, -1, -1, "!", "!"));
                continue;
            }
            IAtom a0 = bond.getAtom(0);
            IAtom a1 = bond.getAtom(1);
            int a0N = atomContainer.indexOf(a0);
            int a1N = atomContainer.indexOf(a1);
            String a0S = a0.getSymbol();
            String a1S = a1.getSymbol();
            int o = bond.getOrder().numeric();
            if (a0N < a1N) {
                edges.add(new Edge(a0N, a1N, o, a0S, a1S));
            } else {
                edges.add(new Edge(a1N, a0N, o, a1S, a0S));
            }
        }
        sort(edges);
        sb.append(edges.toString());
        return sb.toString();
    }

    private class Edge implements Comparable<Edge> {

        public String firstString;
        public String lastString;
        public int first;
        public int last;
        public int order;

        Edge(int first, int last, int order, String firstString, String lastString) {
            this.first = first;
            this.last = last;
            this.order = order;
            this.firstString = firstString;
            this.lastString = lastString;
        }

        @Override
        public int compareTo(Edge o) {
            if (first < o.first || (first == o.first && last < o.last)) {
                return -1;
            } else {
                return 1;
            }
        }

        @Override
        public String toString() {
            return firstString + first + ":" + lastString + last + "(" + order + ")";
        }
    }

}
