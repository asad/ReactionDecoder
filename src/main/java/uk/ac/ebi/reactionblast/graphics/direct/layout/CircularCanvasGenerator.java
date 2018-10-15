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

import java.awt.Dimension;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.Math.toRadians;
import java.util.List;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author asad
 */
public class CircularCanvasGenerator extends
        AbstractCanvasGenerator implements CanvasGenerator {

    private Vector2d vectorToStart;

    // TODO : better way to do this?
    private boolean putFirstInCenter;

    private Dimension size;

    /**
     *
     */
    public CircularCanvasGenerator() {
        this(new Vector2d(-1, 0));
    }

    /**
     *
     * @param vectorToStart
     */
    public CircularCanvasGenerator(Vector2d vectorToStart) {
        this.vectorToStart = vectorToStart;
    }

    /**
     *
     * @param putFirstInCenter
     */
    public CircularCanvasGenerator(boolean putFirstInCenter) {
        this(new Vector2d(-1, 0), putFirstInCenter);
    }

    /**
     *
     * @param vectorToStart
     * @param putFirstInCenter
     */
    public CircularCanvasGenerator(Vector2d vectorToStart, boolean putFirstInCenter) {
        this.vectorToStart = vectorToStart;
        this.putFirstInCenter = putFirstInCenter;
    }

    /**
     *
     * @param atomContainers
     * @param cellCanvas
     */
    @Override
    public void layout(List<IAtomContainer> atomContainers, Dimension cellCanvas) {
        int n;
        if (putFirstInCenter) {
            n = atomContainers.size() - 1;
        } else {
            n = atomContainers.size();
        }
        if (n < 1) {
            return;
        }

        double maxDim = max(cellCanvas.width, cellCanvas.height);
        double alpha = toRadians(360 / n);
        double cosA = cos(alpha);
        double sinA = sin(alpha);
        double circleRadius = (maxDim / 2) / sin(alpha / 2);

        double totalDim = (2 * circleRadius) + maxDim;
        size = new Dimension((int) totalDim, (int) totalDim);
        Point2d center = new Point2d(totalDim / 2, totalDim / 2);

        Vector2d v = new Vector2d(vectorToStart);
        v.normalize();

        int index;
        if (putFirstInCenter) {
            createCanvas(atomContainers.get(0), center, cellCanvas);
            index = 1;
        } else {
            index = 0;
        }
        while (index < atomContainers.size()) {
            IAtomContainer atomContainer = atomContainers.get(index);
            Point2d canvasCenter = new Point2d(center);
            canvasCenter.scaleAdd(circleRadius, v, canvasCenter);
            createCanvas(atomContainer, canvasCenter, cellCanvas);
            Vector2d w = new Vector2d();
            w.x = (cosA * v.x) + (sinA * v.y);
            w.y = (-sinA * v.x) + (cosA * v.y);
            v = w;
            index++;
        }
    }

    /**
     *
     * @return
     */
    @Override
    public Dimension getSize() {
        return size;
    }

}
