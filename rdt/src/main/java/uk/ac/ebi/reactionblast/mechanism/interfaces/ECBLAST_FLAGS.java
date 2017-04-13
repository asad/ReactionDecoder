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
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public enum ECBLAST_FLAGS {

    /**
     * BOND_CHANGE_INFORMATION
     */
    BOND_CHANGE_INFORMATION(0, "Mark Bond Changes"), /**
     * MAX Model
     */
    ATOM_STEREO_CHANGE_INFORMATION(1, "Mark Stereo Changes");
    private final int type;
    private final String description;

    ECBLAST_FLAGS(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of algorithm.
     *
     * @return type of algorithm
     */
    public int type() {
        return this.type;
    }

    /**
     * Returns short description of the algorithm.
     *
     * @return description of the algorithm
     */
    public String description() {
        return this.description;
    }

    /**
     * Compares algorithm types.
     *
     * @param <status>
     * @param obj
     * @return status
     */
    public <status> int compareTo(ECBLAST_FLAGS obj) {
        return 0;
    }
}
