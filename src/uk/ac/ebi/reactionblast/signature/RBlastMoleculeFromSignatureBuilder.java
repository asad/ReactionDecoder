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
package uk.ac.ebi.reactionblast.signature;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

import signature.AbstractGraphBuilder;

/**
 *
 * @author maclean
 *
 */
public class RBlastMoleculeFromSignatureBuilder extends AbstractGraphBuilder {

    /**
     * The chem object builder
     */
    private final IChemObjectBuilder builder;
    /**
     * The container that is being constructed
     */
    private IAtomContainer container;
    private boolean useAromatics = true;
    private boolean useCharge = true;
    private boolean isBondSensitive = true;

    /**
     * Uses the chem object builder for making molecules.
     *
     * @param builder a builder for CDK molecules.
     */
    public RBlastMoleculeFromSignatureBuilder(IChemObjectBuilder builder) {
        this.builder = builder;
    }

    @Override
    /**
     * {@inheritDoc}
     */
    public void makeEdge(int vertexIndex1, int vertexIndex2,
            String vertexSymbol1, String vertexSymbol2, String edgeLabel) {
        if (edgeLabel.equals("") || !isBondSensitive) {
            container.addBond(vertexIndex1, vertexIndex2, IBond.Order.SINGLE);
        } else if (edgeLabel.equals("=")) {
            container.addBond(vertexIndex1, vertexIndex2, IBond.Order.DOUBLE);
        } else if (edgeLabel.equals("#")) {
            container.addBond(vertexIndex1, vertexIndex2, IBond.Order.TRIPLE);
        } else if (edgeLabel.equals("@")) {
            container.addBond(vertexIndex1, vertexIndex2, IBond.Order.SINGLE);
            if (useAromatics) {
                IBond bond = container.getBond(container.getBondCount() - 1);
                bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
                bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);
                bond.setFlag(CDKConstants.ISAROMATIC, true);
            }
        } else if (edgeLabel.equals("%")) {
            container.addBond(vertexIndex1, vertexIndex2, IBond.Order.SINGLE);
            if (useAromatics) {
                IBond bond = container.getBond(container.getBondCount() - 1);
                bond.getAtom(0).setFlag(CDKConstants.ISINRING, true);
                bond.getAtom(1).setFlag(CDKConstants.ISINRING, true);
                bond.setFlag(CDKConstants.ISINRING, true);
            }
        }
    }

    @Override
    /**
     * {@inheritDoc}
     */
    public void makeGraph() {
        this.container = this.builder.newInstance(IAtomContainer.class);
    }

    @Override
    /**
     * {@inheritDoc}
     */
    public void makeVertex(String label) {
        IAtom atom;
        if (label.contains(RBlastAtomSignature.CHARGE_SEPARATOR)) {
            String[] parts = label.split("\\" + RBlastAtomSignature.CHARGE_SEPARATOR);
            atom = this.builder.newInstance(IAtom.class, parts[0]);
            if (useCharge) {
                int charge = Integer.parseInt(parts[1]);
                atom.setFormalCharge(charge);
            }
        } else {
            atom = this.builder.newInstance(IAtom.class, label);
        }

        this.container.addAtom(atom);
    }

    /**
     * Gets the atom container.
     *
     * @return the constructed atom container
     */
    public IAtomContainer getAtomContainer() {
        return this.container;
    }

    /**
     * @return true if aromatic symbols will be used on aromatic bonds
     */
    public boolean isUseAromatics() {
        return useAromatics;
    }

    /**
     * @param useAromatics if true, will use aromatic symbols
     */
    public void setUseAromatics(boolean useAromatics) {
        this.useAromatics = useAromatics;
    }

    /**
     * @return true if charge symbols will be output
     */
    public boolean isUseCharge() {
        return useCharge;
    }

    /**
     * @param useCharge if true, charge symbols will be output
     */
    public void setUseCharge(boolean useCharge) {
        this.useCharge = useCharge;
    }

    /**
     * @return true if bond order symbols are used
     */
    public boolean isBondSensitive() {
        return isBondSensitive;
    }

    /**
     * @param isBondSensitive if true, bond order symbols will be used
     */
    public void setBondSensitive(boolean isBondSensitive) {
        this.isBondSensitive = isBondSensitive;
    }
}
