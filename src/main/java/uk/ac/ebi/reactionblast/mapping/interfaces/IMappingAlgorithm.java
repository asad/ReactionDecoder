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
package uk.ac.ebi.reactionblast.mapping.interfaces;

/**
 *
 * This class represents various algorithm type supported. Presently Atom-Atom
 * Mapping supports 4 different kinds of algorithms:
 *
 * <OL> <lI>0: MIN Minimization Model, <lI>1: MAX Minimization Model, <lI>2: MAX
 * MIXTURE Model, <lI>3: ASSIMILATION Model <lI>4: MIN MIXTURE Model </OL>
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public enum IMappingAlgorithm {

    /**
     * MIN Model
     */
    MIN(0, "Local Minimization Model"),
    /**
     * MAX Model
     */
    MAX(1, "Global Maximization Model"),
    /**
     * MAX-MIN Model
     */
    MIXTURE(2, "Max-Mixture Model"),
    /**
     * RING Conservation Model
     */
    RINGS(3, "Ring Conservation Model"),
    /**
     * RING Conservation Model
     */
    USER_DEFINED(4, "Customised Mapping");
    private final int type;
    private final String description;

    IMappingAlgorithm(int aStatus, String desc) {
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
    public <status> int compareTo(IMappingAlgorithm obj) {
        return 0;
    }
}
