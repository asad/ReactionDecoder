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
package uk.ac.ebi.reactionblast.stereo.wedge;

import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.NONE;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import static org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo.ANTI_CLOCKWISE;

/**
 *
 * @author Gilleain Torrance
 */
public class ACWTetrahedralWedgeRule extends AbstractTetrahedralWedgeRule {

    private final IBond.Stereo[] pattern = {DOWN, UP, NONE, NONE};

    /**
     *
     * @return
     */
    @Override
    public IBond.Stereo[] getPattern() {
        return pattern;
    }

    /**
     *
     * @return
     */
    @Override
    public ITetrahedralChirality.Stereo getStereo() {
        return ANTI_CLOCKWISE;
    }
}
