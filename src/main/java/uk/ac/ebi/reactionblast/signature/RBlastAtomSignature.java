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
package uk.ac.ebi.reactionblast.signature;

import java.util.List;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import signature.AbstractVertexSignature;

/**
 *
 * @author maclean
 *
 */
public class RBlastAtomSignature extends AbstractVertexSignature {

    /**
     *
     */
    public static final String CHARGE_SEPARATOR = ":";
    private final IAtomContainer atomContainer;
    private boolean useAromatics = true;
    private boolean useCharge = true;
    private boolean isBondSensitive = true;

    /**
     *
     * @param atomIndex
     * @param atomContainer
     */
    public RBlastAtomSignature(int atomIndex, IAtomContainer atomContainer) {
        super();
        this.atomContainer = atomContainer;
        super.createMaximumHeight(atomIndex, atomContainer.getAtomCount());
    }

    /**
     *
     * @param atomIndex
     * @param atomContainer
     * @param height
     */
    public RBlastAtomSignature(
            int atomIndex, IAtomContainer atomContainer, int height) {
        super();
        this.atomContainer = atomContainer;
        super.create(atomIndex, atomContainer.getAtomCount(), height);
    }

    /**
     *
     * @return
     */
    public boolean isUseAromatics() {
        return useAromatics;
    }

    /**
     *
     * @param useAromatics
     */
    public void setUseAromatics(boolean useAromatics) {
        this.useAromatics = useAromatics;
    }

    /**
     *
     * @return
     */
    public boolean isUseCharge() {
        return useCharge;
    }

    /**
     *
     * @param useCharge
     */
    public void setUseCharge(boolean useCharge) {
        this.useCharge = useCharge;
    }

    /**
     *
     * @return
     */
    public boolean isBondSensitive() {
        return isBondSensitive;
    }

    /**
     *
     * @param isBondSensitive
     */
    public void setBondSensitive(boolean isBondSensitive) {
        this.isBondSensitive = isBondSensitive;
    }

    /**
     *
     * @param atomIndex
     * @return
     */
    @Override
    public String getVertexSymbol(int atomIndex) {
        IAtom atom = atomContainer.getAtom(atomIndex);
        Integer charge = atom.getFormalCharge();
        if (!useCharge || charge == null || charge == 0) {
            return atom.getSymbol();
        } else {
            return atom.getSymbol() + CHARGE_SEPARATOR + charge;
        }
    }

    /**
     *
     * @param atomIndex
     * @return
     */
    @Override
    protected int[] getConnected(int atomIndex) {
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
        if (useAromatics && bond.getFlag(ISAROMATIC)) {
            return "@";
        } else if (useAromatics && bond.getFlag(ISINRING)) {
            return "%";
        }
        if (!isBondSensitive) {
            return "";
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
     * @param atomIndex
     * @return
     */
    @Override
    protected int getIntLabel(int atomIndex) {
        // TODO Auto-generated method stub
        return 0;
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
