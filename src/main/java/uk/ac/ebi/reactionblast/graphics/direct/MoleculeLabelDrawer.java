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
package uk.ac.ebi.reactionblast.graphics.direct;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.Map;
import javax.vecmath.Point2f;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

/**
 * Draws a set of text labels for one or more molecules.
 *
 * @author maclean
 *
 */
public class MoleculeLabelDrawer extends AbstractDirectDrawer {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MoleculeLabelDrawer.class);

    // TODO: remove axis
    /**
     *
     * @param axis
     * @param params
     */
    public MoleculeLabelDrawer(Axis axis, Params params) {
        setParams(params);
    }

    /**
     *
     * @param labelMap
     * @param labelBounds
     * @param g
     */
    public void draw(Map<String, String> labelMap, BoundsTree labelBounds, Graphics2D g) {
        // layout labels according to the positions in the moleculeBounds
        labelMap.keySet().forEach((boundsLabel) -> {
            String label = labelMap.get(boundsLabel);
            Rectangle2D bounds = labelBounds.get(boundsLabel);
            double x = bounds.getCenterX();
            double y = bounds.getCenterY();
            Point2f p = super.getTextPoint(g, label, x, y);
//            System.out.println("drawing string at " + x + " " + y);
            g.drawString(label, p.x, p.y);
        });
    }
}
