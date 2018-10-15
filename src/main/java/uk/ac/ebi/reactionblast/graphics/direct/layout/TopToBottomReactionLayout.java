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

import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IReaction;
import static uk.ac.ebi.reactionblast.graphics.direct.Axis.Y;

/**
 *
 * @author asad
 */
public class TopToBottomReactionLayout extends AbstractDirectReactionLayout {

    /**
     *
     */
    public TopToBottomReactionLayout() {
        this(true);
    }

    /**
     *
     * @param shouldInvert
     */
    public TopToBottomReactionLayout(boolean shouldInvert) {
        this(shouldInvert, new Vector2d(1, 0));
    }

    /**
     *
     * @param shouldInvert
     * @param moleculeAxis
     */
    public TopToBottomReactionLayout(boolean shouldInvert, Vector2d moleculeAxis) {
        super(shouldInvert, moleculeAxis);
        arrowAxis = Y;
    }

    /**
     *
     * @param reaction
     * @param axis
     * @return
     */
    @Override
    public BoundsTree layout(IReaction reaction, Vector2d axis) {
        Vector2d molSetAxis = new Vector2d(1, 0);
        productBoundsTree = productLayout.layout(reaction.getProducts(), molSetAxis);
        reactantBoundsTree = reactantLayout.layout(reaction.getReactants(), molSetAxis);

        int borderX = params.borderX;
        int borderY = params.borderY;
        int arrowGap = params.arrowGap;
        int arrowLength = params.arrowLength;

        double rbH = reactantBoundsTree.getHeight();
        double pbH = productBoundsTree.getHeight();

        double dx = borderX;
        double dy = borderY + (rbH / 2);
        shiftMoleculeSet(reaction.getReactants(), reactantBoundsTree, dx, dy);
        dy = borderY + rbH + arrowLength + (2 * arrowGap) + (pbH / 2);
        shiftMoleculeSet(reaction.getProducts(), productBoundsTree, dx, dy);

        boundsTree = new BoundsTree(
                "reaction", productBoundsTree, reactantBoundsTree);
        double arrowCenterY = borderY + rbH + arrowGap + (arrowLength / 2);
        arrowPos = arrowCenterY;
        return boundsTree;
    }

    /**
     *
     * @return
     */
    @Override
    public Vector2d getAxis() {
        return new Vector2d(0, 1);
    }

    /**
     *
     * @return
     */
    @Override
    public double getAxisPosition() {
        return (boundsTree.getWidth() / 2) + params.borderX;
    }
}
