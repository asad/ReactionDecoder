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
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import static org.openscience.cdk.geometry.GeometryTools.translate2D;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.graphics.direct.Axis;
import static uk.ac.ebi.reactionblast.graphics.direct.Axis.X;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

/**
 *
 * @author asad
 */
public abstract class AbstractDirectReactionLayout extends AbstractDirectLayout<IReaction> {
    
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
    protected Axis arrowAxis;
    
    /**
     *
     */
    protected double arrowPos;
    
    /**
     *
     */
    protected LinearMoleculeSetLayout reactantLayout;
    
    /**
     *
     */
    protected LinearMoleculeSetLayout productLayout;
    
    /**
     *
     */
    public AbstractDirectReactionLayout() {
        this(true);
    }
    
    /**
     *
     * @param shouldInvert
     */
    public AbstractDirectReactionLayout(boolean shouldInvert) {
        this(shouldInvert, new Vector2d(1, 0));
    }
    
    /**
     *
     * @param shouldInvert
     * @param moleculeAxis
     */
    public AbstractDirectReactionLayout(boolean shouldInvert, Vector2d moleculeAxis) {
        super(shouldInvert);
        reactantLayout = new LinearMoleculeSetLayout(null, shouldInvert, moleculeAxis);
        productLayout = new LinearMoleculeSetLayout(null, shouldInvert, moleculeAxis);
    }

    /**
     *
     * @param params
     */
    public void setParams(Params params) {
        this.params = params;
        reactantLayout.setParams(params);
        productLayout.setParams(params);
    }
    
    /**
     *
     * @return
     */
    public BoundsTree getReactantBounds() {
        return reactantBoundsTree;
    }

    /**
     *
     * @return
     */
    public BoundsTree getProductBounds() {
        return productBoundsTree;
    }
    
    /**
     *
     * @return
     */
    public Vector2d getReactantAxis() {
        return reactantLayout.getAxis();
    }
    
    /**
     *
     * @return
     */
    public Vector2d getProductAxis() {
        return productLayout.getAxis();
    }
    
    /**
     *
     * @return
     */
    public double getReactantAxisPos() {
        return reactantLayout.getAxisPosition();
    }
    
    /**
     *
     * @return
     */
    public double getProductAxisPos() {
        return productLayout.getAxisPosition();
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
     * @param reaction
     * @param axis
     * @param x
     * @param y
     */
    public void shiftReaction(IReaction reaction, Vector2d axis, double x, double y) {
        shiftMoleculeSet(reaction.getReactants(), reactantBoundsTree, x, y);
        shiftMoleculeSet(reaction.getProducts(), productBoundsTree, x + (x * axis.x), y + (y * axis.y));
    }
    
    /**
     *
     * @param molSet
     * @param molSetBoundsTree
     * @param dx
     * @param dy
     */
    public void shiftMoleculeSet(IAtomContainerSet molSet, 
                        BoundsTree molSetBoundsTree, double dx, double dy) {
        int counter = 0;
        String rootLabel = molSet.getID();
        for (IAtomContainer molecule : molSet.atomContainers()) {
            String label = rootLabel + "_" + molecule.getID() + ":" + counter;
            Rectangle2D bounds = molSetBoundsTree.get(label);
            bounds.setFrame(bounds.getCenterX() + dx, bounds.getCenterY() + dy,
                            bounds.getWidth(), bounds.getHeight());
            translate2D(molecule, dx, dy);
            counter++;
        }
    }
    
    /**
     *
     * @return
     */
    public Point2d getArrowCenter() {
        if (arrowAxis == X) {
            return new Point2d(arrowPos, getAxisPosition());
        } else {
            return new Point2d(getAxisPosition(), arrowPos);
        }
    }
    
    /**
     *
     * @return
     */
    public double getArrowPos() {
        return arrowPos;
    }

    /**
     *
     * @param pos
     */
    public void setArrowPos(double pos) {
        arrowPos = pos;
    }
}
