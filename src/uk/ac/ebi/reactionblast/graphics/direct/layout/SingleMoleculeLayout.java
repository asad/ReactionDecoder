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
import java.util.logging.Logger;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import uk.ac.ebi.reactionblast.graphics.direct.Params;

public class SingleMoleculeLayout extends AbstractDirectLayout<IAtomContainer> {
    private static final Logger LOG = Logger.getLogger(SingleMoleculeLayout.class.getName());

    private StructureDiagramGenerator sdg;
    private boolean forceRelayout;

    public SingleMoleculeLayout(Params params) {
        this(params, false);
    }

    public SingleMoleculeLayout(Params params, boolean forceRelayout) {
        setParams(params);
        sdg = new StructureDiagramGenerator();
        this.forceRelayout = forceRelayout;
    }

    @Override
    public BoundsTree layout(IAtomContainer atomContainer, Vector2d axis) {
        // XXX axis is used here to mean center point! :( bad design....
        Point2d center = new Point2d(axis);

        if (forceRelayout || !GeometryTools.has2DCoordinates(atomContainer)) {
            sdg.setMolecule(new AtomContainer(atomContainer), false);
            try {
                if (ConnectivityChecker.isConnected(atomContainer)) {
                    sdg.generateCoordinates();
                } else {
                    System.err.println("Disconnected components needs to be layout separately");
                }
            } catch (CDKException e) {
                e.printStackTrace();
            }
        }
        double scale = GeometryTools.getScaleFactor(atomContainer, params.bondLength);
        Rectangle2D bounds = GeometryTools.getRectangle2D(atomContainer);
        GeometryTools.scaleMolecule(atomContainer, scale);
        translateTo(atomContainer, center.x, center.y, bounds);
        String label = atomContainer.getID();
        return new BoundsTree(label, label, bounds);
    }

    @Override
    public Vector2d getAxis() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public double getAxisPosition() {
        // TODO Auto-generated method stub
        return 0;
    }
}
