/*
 * Copyright (C) 2007-2017 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
 * Copyright (C) 2006-2017 Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public enum ECBLAST_BOND_CHANGE_FLAGS {

    /**
     * BondChangeAnnotation BOND_FORMED_OR_CLEAVED
     */
    BOND_FORMED_OR_CLEAVED(4, "BOND_CHANGE_FORMED_OR_CLEAVED"),
    /**
     * BOND_ORDER
     */
    BOND_ORDER(8, "BOND_CHANGE_ORDER"),
    /**
     * BOND_STEREO
     */
    BOND_STEREO(9, "BOND_CHANGE_STEREO"),
    /**
     * BOND_FORMED
     */
    BOND_FORMED(10, "BOND_CHANGE_FORMED"),
    /**
     * BOND_CLEAVED
     */
    BOND_CLEAVED(11, "BOND_CHANGE_CLEAVED"),
    /**
     * BONDCHANGEORDERGAIN
     */
    BOND_ORDER_GAIN(12, "BOND_CHANGE_ORDER_GAIN"),
    /**
     * BONDCHANGEORDERREDUCED
     */
    BOND_ORDER_REDUCED(13, "BOND_CHANGE_ORDER_REDUCED"),
    /**
     * BONDCHANGEPSEUDOBOND
     */
    PSEUDO_BOND(99, "PSEUDO_BOND_CHANGE");
    private final int type;
    private final String description;

    ECBLAST_BOND_CHANGE_FLAGS(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of algorithm.
     *
     * @return type of algorithm
     */
    public synchronized int type() {
        return this.type;
    }

    /**
     * Returns short description of the algorithm.
     *
     * @return description of the algorithm
     */
    public synchronized String description() {
        return this.description;
    }

    /**
     * Compares algorithm types.
     *
     * @param <status>
     * @param obj
     * @return status
     */
    public synchronized <status> int compareTo(ECBLAST_BOND_CHANGE_FLAGS obj) {
        return 0;
    }
}
