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
package uk.ac.ebi.reactionblast.stereo;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public enum IStereoAndConformation implements Comparable<IStereoAndConformation> {

    /**
     * Default
     */
    NONE(0, "CHIRALITY NONE"),
    /**
     * R.
     */
    R(1, "CHIRALITY Rectus"),
    /**
     * S.
     */
    S(2, "CHIRALITY Sinister"),
    /**
     * Either.
     */
    EITHER(3, "CHIRALITY R or S"),
    /**
     * M
     */
    M(4, "CHIRALITY M Configuration"),
    /**
     * P
     */
    P(5, "CHIRALITY P Configuration"),
    /*
     * as in Z-but-2-ene TOGETHER
     */
    /**
     *
     */
    Z(6, "TOGETHER atom Configuration"),
    /*
     * as in E-but-2-ene OPPOSITE
     */
    /**
     *
     */
    E(7, "OPPOSITE atom Configuration");  // 
    private final int type;
    private final String description;

    IStereoAndConformation(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of stereo.
     *
     * @return type of stereo
     */
    public int type() {
        return this.type;
    }

    /**
     * Returns short description of the stereo.
     *
     * @return description of the stereo
     */
    public String description() {
        return this.description;
    }
}
