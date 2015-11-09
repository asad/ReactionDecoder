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

package uk.ac.ebi.reactionblast.graphics.direct;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.geom.Ellipse2D;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import org.openscience.cdk.interfaces.IAtomContainer;

public class DirectHullDrawer extends AbstractDirectDrawer {
    
    private DirectMoleculeDrawer moleculeDrawer;
    
    public DirectHullDrawer() {
        moleculeDrawer = new DirectMoleculeDrawer();
        setParams(moleculeDrawer.getParams());
    }
    
    public Image drawHull(IAtomContainer atomContainer, int w, int h) {
        Image image = super.makeBlankImage(w, h);
        Graphics2D g = (Graphics2D) image.getGraphics();
        if (params.useAntialias) {
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                               RenderingHints.VALUE_ANTIALIAS_ON);
        }
        for (int i = 0; i < atomContainer.getAtomCount(); i++) {
            atomContainer.getAtom(i).setID(String.valueOf(i));
        }
        params.drawAtomID = true;
        drawHull(atomContainer, g);
        return image;
    }
    
    public void drawHull(IAtomContainer atomContainer, Graphics2D g) {
        DirectArrowDrawer arrowDrawer = new DirectArrowDrawer(getParams());
        ConvexHull hull = new ConvexHull(atomContainer);
        moleculeDrawer.drawMolecule(atomContainer, g);
        Point2d prev = null;
        Point2d first = null;
        
        for (Point2d hullPoint : hull) {
            if (prev == null) {
                prev = hullPoint;
                first = prev;
            } else {
                g.setColor(Color.RED);
                drawLine(prev, hullPoint, g);
                g.setColor(Color.BLACK);
                Point2d midPoint = new Point2d(prev);
                midPoint.interpolate(hullPoint, 0.5);
                Vector2d direction = new Vector2d(hullPoint);
                direction.sub(prev);
                direction.normalize();
                arrowDrawer.drawArrow(g, midPoint, direction);
                
                prev = hullPoint;
            }
        }
        g.setColor(Color.RED);
        drawLine(first, prev, g);
        g.setColor(Color.BLACK);
        Point2d midPoint = new Point2d(prev);
        midPoint.interpolate(first, 0.5);
        Vector2d direction = new Vector2d(first);
        direction.sub(prev);
        direction.normalize();
        arrowDrawer.drawArrow(g, midPoint, direction);
//        ConvexHull.Rectangle r = hull.getMinimumAreaBoundingRectangle();
        ConvexHull.Rectangle r = hull.getMinimumAreaBoundingRectangleBruteForce();
        Vector2d majorAxis = r.getMajorAxis();
        majorAxis.normalize();
        Point2d center = hull.getCenter();
        arrowDrawer.drawArrow(g, center, majorAxis);
        g.setColor(Color.BLACK);
        drawLine(r.cornerA, r.cornerB, g);
        drawLine(r.cornerB, r.cornerC, g);
        drawLine(r.cornerC, r.cornerD, g);
        drawLine(r.cornerD, r.cornerA, g);
        g.setColor(Color.BLUE);
        g.fill(new Ellipse2D.Double(r.cornerA.x - 3, r.cornerA.y - 3, 6, 6));
        g.setColor(Color.MAGENTA);
        g.fill(new Ellipse2D.Double(r.cornerB.x - 3, r.cornerB.y - 3, 6, 6));
        g.setColor(Color.YELLOW);
        g.fill(new Ellipse2D.Double(r.cornerC.x - 3, r.cornerC.y - 3, 6, 6));
        g.setColor(Color.CYAN);
        g.fill(new Ellipse2D.Double(r.cornerD.x - 3, r.cornerD.y - 3, 6, 6));
        
        g.setColor(Color.GREEN);
        g.fill(new Ellipse2D.Double(r.pointX.x - 2, r.pointX.y - 2, 4, 4));
        g.setColor(Color.PINK);
        g.fill(new Ellipse2D.Double(r.pointY.x - 2, r.pointY.y - 2, 4, 4));
        g.setColor(Color.ORANGE);
        g.fill(new Ellipse2D.Double(r.pointZ.x - 2, r.pointZ.y - 2, 4, 4));

    }

}
