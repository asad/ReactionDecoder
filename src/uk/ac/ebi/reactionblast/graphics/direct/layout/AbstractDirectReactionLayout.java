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

import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.graphics.direct.Axis;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

public abstract class AbstractDirectReactionLayout extends AbstractDirectLayout<IReaction> {
    
    protected BoundsTree reactantBoundsTree;
    
    protected BoundsTree productBoundsTree;
    
    protected Axis arrowAxis;
    
    protected double arrowPos;
    
    protected LinearMoleculeSetLayout reactantLayout;
    
    protected LinearMoleculeSetLayout productLayout;
    
    public AbstractDirectReactionLayout() {
        this(true);
    }
    
    public AbstractDirectReactionLayout(boolean shouldInvert) {
        this(shouldInvert, new Vector2d(1, 0));
    }
    
    public AbstractDirectReactionLayout(boolean shouldInvert, Vector2d moleculeAxis) {
        super(shouldInvert);
        reactantLayout = new LinearMoleculeSetLayout(null, shouldInvert, moleculeAxis);
        productLayout = new LinearMoleculeSetLayout(null, shouldInvert, moleculeAxis);
    }

    
    public void setParams(Params params) {
        this.params = params;
        reactantLayout.setParams(params);
        productLayout.setParams(params);
    }
    
    public BoundsTree getReactantBounds() {
        return reactantBoundsTree;
    }

    public BoundsTree getProductBounds() {
        return productBoundsTree;
    }
    
    public Vector2d getReactantAxis() {
        return reactantLayout.getAxis();
    }
    
    public Vector2d getProductAxis() {
        return productLayout.getAxis();
    }
    
    public double getReactantAxisPos() {
        return reactantLayout.getAxisPosition();
    }
    
    public double getProductAxisPos() {
        return productLayout.getAxisPosition();
    }

    public Axis getArrowAxis() {
        return arrowAxis;
    }

    public void shiftReaction(IReaction reaction, Vector2d axis, double x, double y) {
        shiftMoleculeSet(reaction.getReactants(), reactantBoundsTree, x, y);
        shiftMoleculeSet(reaction.getProducts(), productBoundsTree, x + (x * axis.x), y + (y * axis.y));
    }
    
    public void shiftMoleculeSet(IAtomContainerSet molSet, 
                        BoundsTree molSetBoundsTree, double dx, double dy) {
        int counter = 0;
        String rootLabel = molSet.getID();
        for (IAtomContainer molecule : molSet.atomContainers()) {
            String label = rootLabel + "_" + molecule.getID() + ":" + counter;
            Rectangle2D bounds = molSetBoundsTree.get(label);
            bounds.setFrame(bounds.getCenterX() + dx, bounds.getCenterY() + dy,
                            bounds.getWidth(), bounds.getHeight());
            GeometryTools.translate2D(molecule, dx, dy);
            counter++;
        }
    }
    
    public Point2d getArrowCenter() {
        if (arrowAxis == Axis.X) {
            return new Point2d(arrowPos, getAxisPosition());
        } else {
            return new Point2d(getAxisPosition(), arrowPos);
        }
    }
    
    public double getArrowPos() {
        return arrowPos;
    }

    public void setArrowPos(double pos) {
        arrowPos = pos;
    }
}
