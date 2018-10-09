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
package uk.ac.ebi.reactionblast.stereo.inchi;

/**
 * Tool for comparing chiralities.
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public enum InchiStereoAndConformation {

    /**
     * Default
     */
    UNKNOWN(0, "CHIRALITY NONE"),
    /**
     * R.
     */
    MINUS(1, "Negative tetrahedral volume"),
    /**
     * S.
     */
    PLUS(2, "Positive tetrahedral volume"),
    /**
     * NONE
     */
    NONE(3, "NOT A Chiral atom"),
    /*
     * as in Z-but-2-ene
     */
    /**
     *
     */
    TOGETHER(4, "Z-but-2-ene atom"),
    /*
     * as in Z-but-2-ene atom
     */
    /**
     *
     */
    OPPOSITE(5, "E-but-2-ene atom");  //  
    private final int type;
    private final String description;

    InchiStereoAndConformation(int aStatus, String desc) {
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
