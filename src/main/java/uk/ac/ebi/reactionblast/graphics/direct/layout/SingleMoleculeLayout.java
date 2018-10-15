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
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.geometry.GeometryTools.getRectangle2D;
import static org.openscience.cdk.geometry.GeometryTools.getScaleFactor;
import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import static org.openscience.cdk.geometry.GeometryTools.scaleMolecule;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

/**
 *
 * @author asad
 */
public class SingleMoleculeLayout extends AbstractDirectLayout<IAtomContainer> {
     private final static ILoggingTool LOGGER
            = createLoggingTool(SingleMoleculeLayout.class);

    private StructureDiagramGenerator sdg;
    private boolean forceRelayout;

    /**
     *
     * @param params
     */
    public SingleMoleculeLayout(Params params) {
        this(params, false);
    }

    /**
     *
     * @param params
     * @param forceRelayout
     */
    public SingleMoleculeLayout(Params params, boolean forceRelayout) {
        setParams(params);
        sdg = new StructureDiagramGenerator();
        this.forceRelayout = forceRelayout;
    }

    /**
     *
     * @param atomContainer
     * @param axis
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainer atomContainer, Vector2d axis) {
        // XXX axis is used here to mean center point! :( bad design....
        Point2d center = new Point2d(axis);

        if (forceRelayout || !has2DCoordinates(atomContainer)) {
            sdg.setMolecule(new AtomContainer(atomContainer), false);
            try {
                if (isConnected(atomContainer)) {
                    sdg.generateCoordinates();
                } else {
                    LOGGER.debug("Disconnected components needs to be layout separately");
                }
            } catch (CDKException e) {
                e.printStackTrace();
            }
        }
        double scale = getScaleFactor(atomContainer, params.bondLength);
        Rectangle2D bounds = getRectangle2D(atomContainer);
        scaleMolecule(atomContainer, scale);
        translateTo(atomContainer, center.x, center.y, bounds);
        String label = atomContainer.getID();
        return new BoundsTree(label, label, bounds);
    }

    /**
     *
     * @return
     */
    @Override
    public Vector2d getAxis() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     *
     * @return
     */
    @Override
    public double getAxisPosition() {
        // TODO Auto-generated method stub
        return 0;
    }
}
