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


package uk.ac.ebi.reactionblast.graphics.direct.awtlayout;

import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import static org.openscience.cdk.geometry.GeometryTools.translate2D;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.graphics.direct.Axis;
import static uk.ac.ebi.reactionblast.graphics.direct.Axis.X;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

/**
 *
 * @author asad
 */
public abstract class AbstractAWTReactionLayout extends AbstractAWTLayout<IReaction> {
    
    /**
     *
     */
    protected BoundsTree reactantBoundsTree;
    
    /**
     *
     */
    protected BoundsTree productBoundsTree;
    
    /**
     *
     */
    protected Vector2d axis;
    
    /**
     *
     */
    protected Axis arrowAxis;
    
    /**
     *
     */
    protected double arrowPos;
    
    /**
     *
     */
    protected LinearAtomContainerSetLayout reactantLayout;
    
    /**
     *
     */
    protected LinearAtomContainerSetLayout productLayout;
   
    /**
     *
     * @param molSet
     * @param molSetBoundsTree
     * @param dx
     * @param dy
     */
    public void shiftMoleculeSet(IAtomContainerSet molSet, 
            BoundsTree molSetBoundsTree, double dx, double dy) {
//        System.out.println(molSetBoundsTree);
        int counter = 0;
        for (IAtomContainer molecule : molSet.atomContainers()) {
            String molLabel = molSet.getID() + "_" + molecule.getID() + ":" + counter; 
//            System.out.println("shifting " + molLabel + " from " + BoundsPrinter.toString(GeometryTools.getRectangle2D(molecule)));
            Rectangle2D bounds = molSetBoundsTree.get(molLabel);
            bounds.setFrame(bounds.getMinX() + dx, bounds.getMinY() + dy,
                    bounds.getWidth(), bounds.getHeight());
            translate2D(molecule, dx, dy);
//            System.out.println("shifting " + molecule.getID() + " to " + BoundsPrinter.toString(GeometryTools.getRectangle2D(molecule)));
            counter++;
        }
    }
    
    /**
     *
     * @return
     */
    public Axis getArrowAxis() {
        return arrowAxis;
    }
    
    /**
     *
     * @param pos
     */
    public void setArrowPos(double pos) {
        arrowPos = pos;
    }
    
    /**
     *
     * @return
     */
    public Point2d getArrowCenter() {
        Rectangle2D bounds = getBoundsTree().getRoot();
        if (arrowAxis == X) {
            return new Point2d(arrowPos, bounds.getCenterY());
        } else {
            return new Point2d(bounds.getCenterX(), arrowPos);
        }
    }

    /**
     *
     * @return
     */
    public abstract Vector2d getAxis();
  
}
