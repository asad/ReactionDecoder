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

import java.awt.Font;
import static java.awt.Font.PLAIN;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import static java.lang.String.valueOf;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

/**
 *
 * @author asad
 */
public class LinearAtomContainerSetLayout extends AbstractAWTLayout<IAtomContainerSet> {

    private Vector2d moleculeSetAxis;
    private MoleculeLayout moleculeLayout;

    /**
     *
     * @param moleculeSetAxis
     */
    public LinearAtomContainerSetLayout(Vector2d moleculeSetAxis) {
        this(moleculeSetAxis, new Params());
    }

    /**
     *
     * @param moleculeSetAxis
     * @param params
     */
    public LinearAtomContainerSetLayout(Vector2d moleculeSetAxis, Params params) {
        this.params = params;
        moleculeLayout = new MoleculeLayout(params);
        this.moleculeSetAxis = moleculeSetAxis;
    }

    /**
     *
     * @param atomContainerSet
     * @param graphics
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainerSet atomContainerSet, Graphics2D graphics) {
        // TODO : plusBounds.getHeight() for T2B layout
        double molGap = 2 * params.plusGap;

        Font plusFont = new Font("ROMAN", PLAIN, params.plusFontSize);
        graphics.setFont(plusFont);

        if (atomContainerSet.getAtomContainerCount() > 1) {
            Rectangle2D plusBounds;
            plusBounds = super.getTextBounds(graphics, "+");
            molGap += (plusBounds.getWidth());
        }

        String atomContainerSetID = atomContainerSet.getID();
        boundsTree = new BoundsTree(atomContainerSetID);

        Point2d curr = new Point2d(0, 0);
        int moleculeCounter = 0;
        for (IAtomContainer molecule : atomContainerSet.atomContainers()) {
//            System.out.println("curr pos = " + curr.x + " " + curr.y);
            String label = molecule.getID();
            if (label == null || label.isEmpty()) {
                label = "mol" + valueOf(moleculeCounter);
            } else {
                label += ":" + valueOf(moleculeCounter);
            }

            BoundsTree molBounds = moleculeLayout.layout(molecule, label, graphics);

            double boundsWidth = molBounds.getWidth();
            double halfBoundsWidth = boundsWidth / 2;

            curr.scaleAdd(halfBoundsWidth, moleculeSetAxis, curr);
            translateTo(molecule, curr.x, curr.y, molBounds);
            curr.scaleAdd(halfBoundsWidth, moleculeSetAxis, curr);

            curr.scaleAdd(molGap, moleculeSetAxis, curr);

            boundsTree.add(atomContainerSetID, molBounds);
            moleculeCounter++;
        }

        return boundsTree;
    }

    /**
     *
     * @param obj
     * @param rootLabel
     * @param graphics
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainerSet obj, String rootLabel,
            Graphics2D graphics) {
        // TODO Auto-generated method stub
        return null;
    }
}
