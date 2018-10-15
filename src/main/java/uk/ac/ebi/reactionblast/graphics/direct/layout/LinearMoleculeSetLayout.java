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
import static java.lang.String.valueOf;
import static java.util.logging.Level.SEVERE;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.geometry.GeometryTools.getRectangle2D;
import static org.openscience.cdk.geometry.GeometryTools.getScaleFactor;
import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import static org.openscience.cdk.geometry.GeometryTools.scaleMolecule;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

/**
 *
 * @author asad
 */
public class LinearMoleculeSetLayout extends AbstractDirectLayout<IAtomContainerSet> {

    private final static ILoggingTool LOGGER
            = createLoggingTool(LinearMoleculeSetLayout.class);

    /**
     * This is an axis for the individual molecules to be aligned to
     */
    private Vector2d moleculeAxis;

    /**
     *
     * @param params
     */
    public LinearMoleculeSetLayout(Params params) {
        this(params, true);
    }

    /**
     *
     * @param params
     * @param shouldInvert
     */
    public LinearMoleculeSetLayout(Params params, boolean shouldInvert) {
        this(params, shouldInvert, new Vector2d(1, 0));
    }

    /**
     *
     * @param params
     * @param shouldInvert
     * @param moleculeAxis
     */
    public LinearMoleculeSetLayout(
            Params params, boolean shouldInvert, Vector2d moleculeAxis) {
        super(shouldInvert);
        this.moleculeAxis = moleculeAxis;
        setParams(params);
    }

    /**
     *
     * @param atomContainerSet
     * @param moleculeSetAxis
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainerSet atomContainerSet, Vector2d moleculeSetAxis) {
        int bondLength = params.bondLength;
        int molGap = 2 * params.plusGap;

        // if the molecules don't have labels, need to label them
        int molLabel = 0;

        String rootLabel = atomContainerSet.getID();
        boundsTree = new BoundsTree(rootLabel);
        Point2d curr = new Point2d(0, 0);
        int i = 0;
        for (IAtomContainer molecule : atomContainerSet.atomContainers()) {
            if (!has2DCoordinates(molecule)) {
                //Added by Asad for 3D to 2D

                StructureDiagramGenerator sdg
                        = new StructureDiagramGenerator(new AtomContainer(molecule));
                try {
                    sdg.generateCoordinates();
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                molecule = sdg.getMolecule();

            }
            invert(molecule);
            if (params.alignMolecules && moleculeAxis != null) {
                align(molecule, moleculeAxis);
            }
            scaleMolecule(molecule,
                    getScaleFactor(molecule, bondLength));
            Rectangle2D bounds = getRectangle2D(molecule);

            double boundsWidth = bounds.getWidth();
            double halfBoundsWidth = boundsWidth / 2;

            curr.scaleAdd(halfBoundsWidth, moleculeSetAxis, curr);
            translateTo(molecule, curr.x, curr.y, bounds);
            curr.scaleAdd(halfBoundsWidth, moleculeSetAxis, curr);
            curr.scaleAdd(molGap, moleculeSetAxis, curr);

            String moleculeLabel = molecule.getID();
            if (moleculeLabel == null || moleculeLabel.isEmpty()) {
                moleculeLabel = "mol" + valueOf(molLabel);
                molLabel++;
            } else {
                moleculeLabel += ":" + i;
            }

            boundsTree.add(rootLabel + "_" + moleculeLabel, bounds);
            i++;
            shouldInvert = true;
        }
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
        return (boundsTree.getWidth() / 2) + params.borderX;
    }
}
