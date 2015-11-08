/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.logging.Logger;
import javax.vecmath.Point2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;

public class ZoomToFitGridLayout {
    private static final Logger LOG = Logger.getLogger(ZoomToFitGridLayout.class.getName());

    private DirectMoleculeDrawer drawer;

    private int rows;

    private int cols;

    public ZoomToFitGridLayout() {
        drawer = new DirectMoleculeDrawer();
    }

    public ZoomToFitGridLayout(int rows, int cols) {
        this(new DirectMoleculeDrawer(), rows, cols);
    }

    public ZoomToFitGridLayout(DirectMoleculeDrawer drawer, int rows, int cols) {
        this.drawer = drawer;
        this.rows = rows;
        this.cols = cols;
    }

    public void layout(
            List<IAtomContainer> mols, Dimension cellCanvas, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();
        double w = cellCanvas.width;
        double h = cellCanvas.height;

        double centerX = w / 2;
        double centerY = h / 2;
        int colCounter = 1;
        for (IAtomContainer mol : mols) {
            double zoom = calculateZoom(mol, cellCanvas);
            g.translate(centerX, centerY);
            g.scale(zoom, zoom);
            drawer.drawMolecule(mol, g);
            g.setTransform(originalTransform);
            if (colCounter < cols) {
                centerX += w;
                colCounter++;
            } else {
                centerY += h;
                centerX = w / 2;
                colCounter = 1;
            }
        }
    }

    private double calculateZoom(IAtomContainer ac, Dimension canvas) {
        double scaleFactor
                = GeometryTools.getScaleFactor(ac, drawer.getParams().bondLength);
        GeometryTools.translate2DCenterTo(ac, new Point2d(0, 0));
        GeometryTools.scaleMolecule(ac, scaleFactor);
        Rectangle2D r2D = GeometryTools.getRectangle2D(ac);
        double canvasWidth = canvas.width;
        double canvasHeight = canvas.height;
        double borderX = drawer.getParams().borderX;
        double borderY = drawer.getParams().borderY;
        double objectWidth = r2D.getWidth() + (borderX * 2);
        double objectHeight = r2D.getHeight() + (borderY * 2);
        return Math.min(canvasWidth / objectWidth, canvasHeight / objectHeight);
    }

}
