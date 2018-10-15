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
package uk.ac.ebi.reactionblast.graphics.direct.layout;

import java.awt.geom.Rectangle2D;
import static java.lang.String.format;

/**
 *
 * @author asad
 */
public class BoundsPrinter {

    /**
     * DEBUG method for printing readable rectangle 2Ds
     *
     * @param b
     * @return
     */
    public static String toString(Rectangle2D b) {
        return format("[(%2.0f, %2.0f), (%2.0f, %2.0f)] "
                + "= (%2.0f x %2.0f) @ [%2.0f, %2.0f]",
                b.getMinX(), b.getMinY(), b.getMaxX(), b.getMaxY(),
                b.getWidth(), b.getHeight(), b.getCenterX(), b.getCenterY());
    }

    private BoundsPrinter() {
    }

}
