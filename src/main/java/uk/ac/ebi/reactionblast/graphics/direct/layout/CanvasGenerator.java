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

import java.awt.Dimension;
import java.awt.geom.Rectangle2D;
import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author asad
 */
public interface CanvasGenerator {

    /**
     *
     * @param atomContainers
     * @param cellCanvas
     */
    public void layout(List<IAtomContainer> atomContainers, Dimension cellCanvas);

    /**
     *
     * @param atomContainer
     * @return
     */
    public Rectangle2D getCanvasForAtomContainer(IAtomContainer atomContainer);

    /**
     *
     * @return
     */
    public Dimension getSize();

}
