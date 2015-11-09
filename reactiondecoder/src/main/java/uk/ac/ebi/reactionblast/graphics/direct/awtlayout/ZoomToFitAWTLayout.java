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
package uk.ac.ebi.reactionblast.graphics.direct.awtlayout;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import javax.vecmath.Vector2d;

import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.layout.AbstractDirectLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

public class ZoomToFitAWTLayout extends AbstractDirectLayout<IAtomContainer> {

    private DirectMoleculeDrawer drawer;

    public ZoomToFitAWTLayout(DirectMoleculeDrawer drawer) {
        this.drawer = drawer;
        super.setParams(drawer.getParams());
    }

    public void layout(IAtomContainer mol, Dimension cellCanvas, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();
        double cW = cellCanvas.width;
        double cH = cellCanvas.height;

        if (shouldInvert) {
            super.invert(mol);
        }

        BoundsTree tree = getBoundsTree(mol, g);
        double tW = tree.getWidth();
        double tH = tree.getHeight();

        // work out the label height, to get the zoom correctly
        Rectangle2D stringBounds = null;
        String label = mol.getID();
        Font labelFont = new Font(
                params.labelPanelFont, Font.PLAIN, params.labelPanelFontSize);
        if (params.drawLabelPanel) {
            g.setFont(labelFont);
            FontMetrics metrics = g.getFontMetrics();
            stringBounds = metrics.getStringBounds(label, g);
            double labelHeight = stringBounds.getHeight();
            cH += labelHeight;
        }

        double zoom = calculateZoom(tW, tH, cW, cH);

        // adjust the center to correct for labels
        double centerX = cW / 2;
        double centerY;
        Params params = drawer.getParams();
        if (params.drawMoleculeID) {
            centerY = (cH / 2) - params.labelYGap;
        } else if (params.drawLabelPanel) {
            double labelHeight = stringBounds.getHeight();
            double scaledLabelHeight = (labelHeight / 2) * (1 / zoom);
            centerY = ((double) cellCanvas.height / 2) - scaledLabelHeight;
//            centerY = cH / 2;
        } else {
            centerY = cH / 2;
        }

        g.translate(centerX, centerY);
        g.scale(zoom, zoom);

        drawer.drawMolecule(mol, g);

        // DEBUG
//        g.setColor(Color.red);
//        g.draw(tree.getRoot());
        g.setTransform(originalTransform);

        if (params.drawLabelPanel) {
            double cX = cW / 2;
            double cY = cH / 2;
            g.setFont(labelFont);
            FontMetrics metrics = g.getFontMetrics();
            double halfWidth = stringBounds.getWidth() / 2;
            double halfHeight = stringBounds.getHeight() / 2;
            double halfScaledTreeWidth = (tH * zoom) / 2;
            double lY = cY + halfScaledTreeWidth - params.borderY;
            double ascent = metrics.getAscent();

            // DEBUG
//            g.setColor(Color.BLUE);
//            g.draw(new Rectangle2D.Double(
//                    cX - halfWidth, lY - halfHeight, halfWidth * 2, halfHeight * 2));
            float x = (float) (cX - halfWidth);
            float y = (float) (lY - halfHeight + ascent);

            g.setColor(Color.BLACK);
//            System.out.println("drawing label " + label + " at " + x + " " + y);
            g.drawString(label, x, y);
        }
    }

    private BoundsTree getBoundsTree(IAtomContainer mol, Graphics2D g) {
        Rectangle2D bb = GeometryTools.getRectangle2D(mol);
        GeometryTools.translate2D(mol, -bb.getCenterX(), -bb.getCenterY());
        GeometryTools.scaleMolecule(mol,
                GeometryTools.getScaleFactor(mol, params.bondLength));
        MoleculeLayout exactLayout = new MoleculeLayout(params);
        return exactLayout.layout(mol, g);
    }

    private double calculateZoom(double tw, double th, double cw, double ch) {
        Params params = drawer.getParams();
        double borderX = params.borderX;
        double borderY = params.borderY;
//        System.out.println("border " + borderX + " " + borderY);
        double rW = tw + (borderX * 2);
        double rH = th + (borderY * 2);
        return Math.min(cw / rW, ch / rH);
    }

    @Override
    public BoundsTree layout(IAtomContainer obj, Vector2d axis) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Vector2d getAxis() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public double getAxisPosition() {
        // TODO Auto-generated method stub
        return 0;
    }

}
