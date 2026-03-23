/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism.helper;

import java.io.Serializable;

import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class BondChange implements Serializable {

    private static final String NEW_LINE = System.lineSeparator();
    private static final long serialVersionUID = 9890766688070991L;

    /**
     *
     * @param bond
     * @return
     */
    public static int convertBondOrder(IBond bond) {
        if (bond.getOrder() == null) {
            return bond.isAromatic() ? 2 : 1;
        }
        switch (bond.getOrder()) {
            case QUADRUPLE:
                return 4;
            case TRIPLE:
                return 3;
            case DOUBLE:
                return 2;
            case SINGLE:
                return 1;
            default:
                // Handle UNSET or other cases — check aromaticity flag
                return bond.isAromatic() ? 2 : 1;
        }
    }

    /**
     *
     * @param bond
     * @return
     */
    public static int convertBondStereo(IBond bond) {
        int value;
        switch (bond.getStereo()) {
            case UP:
                value = 1;
                break;
            case UP_INVERTED:
                value = 1;
                break;
            case DOWN:
                value = 6;
                break;
            case DOWN_INVERTED:
                value = 6;
                break;
            case UP_OR_DOWN:
                value = 4;
                break;
            case UP_OR_DOWN_INVERTED:
                value = 4;
                break;
            case E_OR_Z:
                value = 3;
                break;
            default:
                value = 0;
        }
        return value;
    }

    private final IBond reactantBond;
    private final IBond productBond;
    private final float bondChangeDelta;

    /**
     *
     * @param reactantBond
     * @param productBond
     */
    public BondChange(IBond reactantBond, IBond productBond) {
        this.reactantBond = reactantBond;
        this.productBond = productBond;
        if (this.reactantBond != null && this.productBond != null) {
            this.bondChangeDelta = convertBondOrder(this.productBond) - convertBondOrder(this.reactantBond);
        } else if (this.reactantBond == null && this.productBond != null) {
            this.bondChangeDelta = convertBondOrder(this.productBond);
        } else if (this.reactantBond != null && this.productBond == null) {
            this.bondChangeDelta = convertBondOrder(this.reactantBond);
        } else {
            this.bondChangeDelta = 0;
        }
    }

    /**
     * @return the reactantBond
     */
    public IBond getReactantBond() {
        return reactantBond;
    }

    /**
     * @return the productBond
     */
    public IBond getProductBond() {
        return productBond;
    }

    /**
     * @return the bondChangeDelta
     */
    public float getBondChangeDelta() {
        return bondChangeDelta;
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append("\t");
        result.append(NEW_LINE);
        if (reactantBond != null) {
            result.append("R: ").append(reactantBond.getAtom(0).getSymbol());
            result.append("(").append(reactantBond.getAtom(0).getID()).append(")");
            result.append("[").append(convertBondOrder(reactantBond)).append("]");
            result.append(reactantBond.getAtom(1).getSymbol());
            result.append("(").append(reactantBond.getAtom(1).getID()).append(")");

        } else {
            result.append("NA");
        }

        if (productBond != null) {
            result.append(", P: ").append(productBond.getAtom(0).getSymbol());
            result.append("(").append(productBond.getAtom(0).getID()).append(")");
            result.append("[").append(convertBondOrder(productBond)).append("]");
            result.append(productBond.getAtom(1).getSymbol());
            result.append("(").append(productBond.getAtom(1).getID()).append(")");
            result.append(NEW_LINE);
        } else {
            result.append(", NA");
            result.append(NEW_LINE);
        }
        return result.toString();
    }

}
