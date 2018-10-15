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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.Serializable;
import static java.lang.System.getProperty;

import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondChange implements Serializable {

    private static final long serialVersionUID = 9890766688070991L;

    /**
     *
     * @param bond
     * @return
     */
    public static synchronized int convertBondOrder(IBond bond) {
        int value;
        switch (bond.getOrder()) {
            case QUADRUPLE:
                value = 4;
                break;
            case TRIPLE:
                value = 3;
                break;
            case DOUBLE:
                value = 2;
                break;
            case SINGLE:
                value = 1;
                break;
            default:
                value = 1;
        }
        return value;
    }

    /**
     *
     * @param bond
     * @return
     */
    public static synchronized int convertBondStereo(IBond bond) {
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
    private float bondChangeDelta;

    /**
     *
     * @param reactantBond
     * @param productBond
     */
    public BondChange(IBond reactantBond, IBond productBond) {
        this.bondChangeDelta = 0;
        this.reactantBond = reactantBond;
        this.productBond = productBond;
        if (this.reactantBond != null && this.productBond != null) {
            float change = convertBondOrder(this.productBond) - convertBondOrder(this.reactantBond);
            if (change != 0) {
                setBondChangeInformation(change);
            } else if (change == 0) {
                setBondChangeInformation(0);
            }
        } else if (this.reactantBond == null && this.productBond != null) {
            setBondChangeInformation(convertBondOrder(this.productBond));
        } else if (this.reactantBond != null && this.productBond == null) {
            setBondChangeInformation(convertBondOrder(this.reactantBond));
        }
    }

    /**
     * @return the reactantBond
     */
    public synchronized IBond getReactantBond() {
        return reactantBond;
    }

    /**
     * @return the productBond
     */
    public synchronized IBond getProductBond() {
        return productBond;
    }

    /**
     * @return the bondChangeDelta
     */
    public synchronized float getBondChangeDelta() {
        return bondChangeDelta;
    }

    /**
     * @param bondChangeDelta the bondChangeDelta to set
     */
    private synchronized void setBondChangeInformation(float bondChangeInformation) {
        this.bondChangeDelta = bondChangeInformation;
    }

    @Override
    public synchronized String toString() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");
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
