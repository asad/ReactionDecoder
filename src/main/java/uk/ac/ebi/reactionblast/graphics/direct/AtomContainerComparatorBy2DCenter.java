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
package uk.ac.ebi.reactionblast.graphics.direct;

import java.util.Comparator;

import javax.vecmath.Point2d;
import static org.openscience.cdk.geometry.GeometryTools.get2DCenter;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class AtomContainerComparatorBy2DCenter implements Comparator<IAtomContainer> {

    /**
     * Compare two AtomContainers based on their 2D position.
     *
     * @param atCont1
     * @param atCont2
     * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
     */
    @Override
    public int compare(IAtomContainer atCont1, IAtomContainer atCont2) {
        if (atCont1 != null && atCont2 != null) {
            Point2d p1 = get2DCenter(atCont1);
            Point2d p2 = get2DCenter(atCont2);
            if (p1 != null && p2 != null) {
                if (p1.x != p2.x) {
                    return new Double(p1.x).compareTo(p2.x);
                } else {
                    return new Double(p1.y).compareTo(p2.y);
                }
            }
        }
        return 0;
    }

}
