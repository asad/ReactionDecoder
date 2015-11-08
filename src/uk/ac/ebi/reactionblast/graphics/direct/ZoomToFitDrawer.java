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

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.logging.Logger;
import javax.vecmath.Point2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.CanvasGenerator;
import uk.ac.ebi.reactionblast.graphics.direct.layout.GridCanvasGenerator;

public class ZoomToFitDrawer {
    private static final Logger LOG = Logger.getLogger(ZoomToFitDrawer.class.getName());
    
    private DirectMoleculeDrawer moleculeDrawer;
    
    private CanvasGenerator canvasGenerator;
    
    private Params params;
    
    /**
     *
     */
    public ZoomToFitDrawer() {
        this(new DirectMoleculeDrawer(), new GridCanvasGenerator());
    }
    
    /**
     *
     * @param moleculeDrawer
     * @param canvasGenerator
     */
    public ZoomToFitDrawer(DirectMoleculeDrawer moleculeDrawer,
            CanvasGenerator canvasGenerator) {
        this.moleculeDrawer = moleculeDrawer;
        this.params = moleculeDrawer.getParams();
        this.canvasGenerator = canvasGenerator;
    }
    
    /**
     *
     * @param mols
     * @param cellCanvas
     * @param g
     */
    public void draw(List<IAtomContainer> mols,  Dimension cellCanvas, Graphics2D g) {
        canvasGenerator.layout(mols, cellCanvas);
        AffineTransform originalTransform = g.getTransform();
        for (IAtomContainer mol : mols) {
            Rectangle2D canvas = canvasGenerator.getCanvasForAtomContainer(mol);
            g.translate(canvas.getCenterX(), canvas.getCenterY());
            double zoom = calculateZoom(mol, canvas);
            g.scale(zoom, zoom);
            
            moleculeDrawer.drawMolecule(mol, g);
            g.setTransform(originalTransform);
        }
    }
    
    private double calculateZoom(IAtomContainer ac, Rectangle2D canvas) {
        double scaleFactor = GeometryTools.getScaleFactor(ac, params.bondLength);
        GeometryTools.translate2DCenterTo(ac, new Point2d(0,0));
        GeometryTools.scaleMolecule(ac, scaleFactor);
        Rectangle2D r2D = GeometryTools.getRectangle2D(ac);
        double canvasWidth = canvas.getWidth();
        double canvasHeight = canvas.getHeight();
        double objectWidth = r2D.getWidth() + (params.borderX * 2);
        double objectHeight = r2D.getHeight() + (params.borderY * 2);
        return Math.min(canvasWidth / objectWidth, canvasHeight / objectHeight);
    }

}
