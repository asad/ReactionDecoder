/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.blocks;

import static java.lang.System.out;
import java.util.List;
import java.util.Map;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import signature.AbstractVertexSignature;

/**
 *
 * @author asad
 */
public class SubgraphAtomSignature extends AbstractVertexSignature {

    private final IAtomContainer atomContainer;
    private Map<Integer, int[]> subgraphAdjacencyLists;
    // XXX - tmp fix to get around the architectural bugs of using adj lists.
    private boolean useAdjLists = false;

    /**
     *
     * @param atomContainer
     * @param atomIndex
     */
    public SubgraphAtomSignature(IAtomContainer atomContainer, int atomIndex) {
        super();
        this.atomContainer = atomContainer;
        super.createMaximumHeight(atomIndex, atomContainer.getAtomCount());
    }

    /**
     *
     * @param atomContainer
     * @param atomIndex
     * @param subgraphAdjacencyLists
     */
    public SubgraphAtomSignature(IAtomContainer atomContainer, int atomIndex,
            Map<Integer, int[]> subgraphAdjacencyLists) {
        super();

        useAdjLists = true; // NOT RECOMMENDED YET...

        this.atomContainer = atomContainer;
        this.subgraphAdjacencyLists = subgraphAdjacencyLists;
        int graphVertexCount = atomContainer.getAtomCount();
        super.createMaximumHeight(atomIndex, graphVertexCount);
    }

    /**
     *
     * @param atomIndex
     * @return
     */
    @Override
    protected int[] getConnected(int atomIndex) {
        if (useAdjLists) {
            if (subgraphAdjacencyLists.containsKey(atomIndex)) {
                return subgraphAdjacencyLists.get(atomIndex);
            } else {
                return new int[]{};
            }
        } else {

            // XXX TMP
            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(
                    atomContainer.getAtom(atomIndex));
            int[] connected = new int[connectedAtoms.size()];
            int i = 0;
            for (IAtom connectedAtom : connectedAtoms) {
                connected[i] = atomContainer.indexOf(connectedAtom);
                i++;
            }
            return connected;
        }
    }

    /**
     *
     * @param atomIndexA
     * @param atomIndexB
     * @return
     */
    @Override
    protected String getEdgeLabel(int atomIndexA, int atomIndexB) {
        IAtom atomA = atomContainer.getAtom(atomIndexA);
        IAtom atomB = atomContainer.getAtom(atomIndexB);
        IBond bond = atomContainer.getBond(atomA, atomB);
        if (bond.getFlag(ISAROMATIC)) {
            return "@";
        } else if (bond.getFlag(ISINRING)) {
            return "%";
        }
        switch (bond.getOrder()) {
            case SINGLE:
                return "";
            case DOUBLE:
                return "=";
            case TRIPLE:
                return "#";
            case QUADRUPLE:
                return "$";
            default:
                return "";
        }
    }

    /**
     *
     * @param arg0
     * @return
     */
    @Override
    protected int getIntLabel(int arg0) {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     *
     * @param atomIndex
     * @return
     */
    @Override
    protected String getVertexSymbol(int atomIndex) {
        IAtom atom = atomContainer.getAtom(atomIndex);
        try {
            Integer charge = atom.getFormalCharge();
            if (charge == null || charge == 0) {
                return atom.getSymbol();
            } else {
                return atom.getSymbol() + charge;
            }
        } catch (NullPointerException npe) {
            //XXX remove!
            out.println("npe getting atom " + atomIndex
                    + " from " + atomContainer.getID());
            npe.printStackTrace();
            return "C";
        }
    }

    /**
     *
     * @param edgeLabel
     * @return
     */
    @Override
    protected int convertEdgeLabelToColor(String edgeLabel) {
        switch (edgeLabel) {
            case "":
                return 1;
            case "=":
                return 2;
            case "#":
                return 3;
            case "$":
                return 4;
            case "@":
                return 5;
            case "%":
                return 6;
        }
        return 0;
    }
}
