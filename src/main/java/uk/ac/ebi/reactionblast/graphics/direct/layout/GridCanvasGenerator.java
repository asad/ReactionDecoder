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
import java.util.List;
import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 *
 * @author asad
 */
public class GridCanvasGenerator extends AbstractCanvasGenerator implements CanvasGenerator {

    private final static ILoggingTool LOGGER
            = createLoggingTool(GridCanvasGenerator.class);

    private int rows;

    private int cols;

    private Dimension size;

    /**
     *
     */
    public GridCanvasGenerator() {
        this(1, 1); // hmmm...
    }

    /**
     *
     * @param rows
     * @param cols
     */
    public GridCanvasGenerator(int rows, int cols) {
        super();
        this.rows = rows;
        this.cols = cols;
    }

    /**
     *
     * @param atomContainers
     * @param cellCanvas
     */
    @Override
    public void layout(List<IAtomContainer> atomContainers, Dimension cellCanvas) {
        double w = cellCanvas.width;
        double h = cellCanvas.height;
        double centerX = w / 2;
        double centerY = h / 2;
        int colCounter = 0;
        int rowCounter = 0;
        for (IAtomContainer atomContainer : atomContainers) {
            createCanvas(atomContainer, new Point2d(centerX, centerY), cellCanvas);
            colCounter++;
            if (colCounter < cols) {
                centerX += w;
            } else {
                centerY += h;
                centerX = w / 2;
                colCounter = 0;
                rowCounter++;
            }

            if (rowCounter > rows) {
                LOGGER.debug("WARNING : Row limit exceeded");
            }
        }
        size = new Dimension(cols * cellCanvas.width, rows * cellCanvas.height);
    }

    /**
     *
     * @return
     */
    @Override
    public Dimension getSize() {
        return size;
    }
}
