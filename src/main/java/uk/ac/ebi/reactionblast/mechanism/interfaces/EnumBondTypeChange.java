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
package uk.ac.ebi.reactionblast.mechanism.interfaces;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
/**
 * BOND TYPE CHANGE (FORMED/BROKEN, ORDER CHANGE, STEREOCHANGE, NO CHANGE,
 * UNKNOWN)
 */
public enum EnumBondTypeChange {

    /**
     *
     */
    SINGLE_BOND_BROKEN_FORMED(0, "SINGLE_BOND_BROKEN_FORMED"),
    /**
     *
     */
    DOUBLE_BOND_BROKEN_FORMED(1, "DOUBLE_BOND_BROKEN_FORMED"),
    /**
     *
     */
    TRIPLE_BOND_BROKEN_FORMED(2, "TRIPLE_BOND_BROKEN_FORMED"),
    /**
     *
     */
    SINGLE_BOND_ORDER_CHANGE(3, "SINGLE_BOND_ORDER_CHANGE"),
    /**
     *
     */
    DOUBLE_BOND_ORDER_CHANGE(4, "DOUBLE_BOND_ORDER_CHANGE"),
    /**
     *
     */
    TRIPLE_BOND_ORDER_CHANGE(5, "TRIPLE_BOND_ORDER_CHANGE"),
    /**
     *
     */
    BOND_STEREO_CHANGE(6, "BOND_STEREO_CHANGE"),
    /**
     *
     */
    NO_CHANGE(7, "NO_CHANGE"),
    /**
     *
     */
    UNKNOWN_CHANGE(8, "UNKNOWN_CHANGE");
    private final int type;
    private final String description;

    EnumBondTypeChange(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     *
     * @return
     */
    public synchronized int type() {
        return this.type;
    }

    /**
     *
     * @return
     */
    public synchronized String description() {
        return this.description;
    }

    /**
     *
     * @param <status>
     * @param obj
     * @return
     */
    public synchronized <status> int compareTo(EnumBondTypeChange obj) {
        return 0;
    }
}
