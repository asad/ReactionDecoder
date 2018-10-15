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
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import static java.lang.Math.min;

import javax.vecmath.Vector2d;
import static org.openscience.cdk.geometry.GeometryTools.getRectangle2D;
import static org.openscience.cdk.geometry.GeometryTools.getScaleFactor;
import static org.openscience.cdk.geometry.GeometryTools.scaleMolecule;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

/**
 * Slightly broken class to layout molecules at different scales, but constant
 * canvas 'cell' size.
 *
 * It doesn't properly implement layout(T, Vector2d) as it needs a Dimension
 * (the canvas) and a Graphics object to do the scaling...
 *
 * @author maclean
 *
 */
public class ZoomToFitLayout extends AbstractDirectLayout<IAtomContainer> {

    private final DirectMoleculeDrawer drawer;

    /**
     *
     * @param drawer
     */
    public ZoomToFitLayout(DirectMoleculeDrawer drawer) {
        this.drawer = drawer;
    }

    /**
     *
     * @param mol
     * @param cellCanvas
     * @param g
     */
    public void layout(IAtomContainer mol, Dimension cellCanvas, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();
        double w = cellCanvas.width;
        double h = cellCanvas.height;

        double zoom = calculateZoom(mol, cellCanvas.width, cellCanvas.height);

        double centerX = w / 2;
        double centerY;
        Params params = drawer.getParams();
        if (params.drawMoleculeID) {
            centerY = (h / 2) - params.labelYGap;
        } else {
            centerY = h / 2;
        }

        g.translate(centerX, centerY);
        g.scale(zoom, zoom);
        drawer.drawMolecule(mol, g);
        g.setTransform(originalTransform);
    }

    private double calculateZoom(IAtomContainer ac, double w, double h) {
        double borderX = drawer.getParams().borderX;
        double borderY = drawer.getParams().borderY;
//        System.out.println("border " + borderX + " " + borderY);
        double canvasWidth = w;
        double canvasHeight = h;
        double scaleFactor
                = getScaleFactor(ac, drawer.getParams().bondLength);
        Rectangle2D r2D = getRectangle2D(ac);
//        Rectangle2D tmp = new Rectangle2D.Double(r2D.getMinX(), r2D.getMinY(), r2D.getWidth(), r2D.getHeight());
        translateTo(ac, 0, 0, r2D);
//        translateTo(ac, 0, 0, tmp);
        scaleMolecule(ac, scaleFactor);
        double objectWidth = r2D.getWidth() + (borderX * 2);
        double objectHeight = r2D.getHeight() + (borderY * 2);

        return min(canvasWidth / objectWidth, canvasHeight / objectHeight);
    }

    /**
     *
     * @param obj
     * @param axis
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainer obj, Vector2d axis) {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     *
     * @return
     */
    @Override
    public Vector2d getAxis() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     *
     * @return
     */
    @Override
    public double getAxisPosition() {
        // TODO Auto-generated method stub
        return 0;
    }

}
