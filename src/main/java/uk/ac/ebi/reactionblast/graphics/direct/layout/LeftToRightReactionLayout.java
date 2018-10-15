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

import static java.lang.Math.max;

import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static uk.ac.ebi.reactionblast.graphics.direct.Axis.X;

/**
 * Layout atoms containers by modifying the points in the molecules.
 *
 * @author maclean
 *
 */
public class LeftToRightReactionLayout extends AbstractDirectReactionLayout {

    /**
     *
     */
    public LeftToRightReactionLayout() {
        this(true);
    }

    /**
     *
     * @param shouldLayout
     */
    public LeftToRightReactionLayout(boolean shouldLayout) {
        this(shouldLayout, new Vector2d(0, 1));
    }

    /**
     *
     * @param shouldLayout
     * @param moleculeAxis
     */
    public LeftToRightReactionLayout(boolean shouldLayout, Vector2d moleculeAxis) {
        super(shouldLayout, moleculeAxis);
        arrowAxis = X;
    }

    /**
     *
     * @param reaction
     * @param axis
     * @return
     */
    @Override
    public BoundsTree layout(IReaction reaction, Vector2d axis) {
        IAtomContainerSet reactants = reaction.getReactants();
        reactants.setID("r");
        reactantBoundsTree = reactantLayout.layout(reactants, axis);

        IAtomContainerSet products = reaction.getProducts();
        products.setID("p");
        productBoundsTree = productLayout.layout(products, axis);

        int borderX = params.borderX;
        int borderY = params.borderY;
        int arrowGap = params.arrowGap;
        int arrowLength = params.arrowLength;

        double rbH = reactantBoundsTree.getHeight();
        double pbH = productBoundsTree.getHeight();
        double rbW = reactantBoundsTree.getWidth();
        double maxH = max(rbH, pbH);

        double dx = borderX;
        double dy = borderY + (maxH / 2);
        shiftMoleculeSet(reaction.getReactants(), reactantBoundsTree, dx, dy);
        dx = borderX + rbW + arrowLength + (2 * arrowGap);
        shiftMoleculeSet(reaction.getProducts(), productBoundsTree, dx, dy);

        boundsTree = new BoundsTree(
                reaction.getID(), productBoundsTree, reactantBoundsTree);
        double arrowCenterX = borderX + rbW + arrowGap + (arrowLength / 2);
        arrowPos = arrowCenterX;

        return boundsTree;
    }

    /**
     *
     * @return
     */
    @Override
    public Vector2d getAxis() {
        return new Vector2d(1, 0);
    }

    /**
     *
     * @return
     */
    @Override
    public double getAxisPosition() {
        return (boundsTree.getHeight() / 2) + params.borderY;
    }

}
