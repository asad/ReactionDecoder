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

import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import javax.vecmath.Point2f;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

/**
 * A layout object that uses a java.awt.Graphics2D instance to work out text sizes, which makes the layout more precise
 * than just using atom positions.
 *
 * @author maclean
 * @param <T>
 *
 */
public abstract class AbstractAWTLayout<T> {

    /**
     *
     */
    protected Graphics2D graphics;

    /**
     *
     */
    protected AbstractAWTLayout parent;

    /**
     *
     */
    protected T currentObject;

    /**
     *
     */
    protected Params params;

    /**
     *
     */
    protected BoundsTree boundsTree;

    /**
     *
     * @return
     */
    public Params getParams() {
        return params;
    }

    /**
     *
     * @param params
     */
    public void setParams(Params params) {
        this.params = params;
    }

    /**
     *
     * @param obj
     * @param graphics
     * @return
     */
    public abstract BoundsTree layout(T obj, Graphics2D graphics);

    /**
     *
     * @param obj
     * @param rootLabel
     * @param graphics
     * @return
     */
    public abstract BoundsTree layout(T obj, String rootLabel, Graphics2D graphics);

    /**
     *
     * @return
     */
    public BoundsTree getBoundsTree() {
        return boundsTree;
    }

    /**
     *
     * @return
     */
    public T getCurrentObject() {
        return currentObject;
    }

    /**
     *
     * @return
     */
    public Graphics2D getGraphics() {
        return graphics;
    }

    /**
     *
     * @param graphics
     */
    public void setGraphics(Graphics2D graphics) {
        this.graphics = graphics;
    }

    /**
     *
     * @param g
     * @param text
     * @param cX
     * @param cY
     * @return
     */
    public Point2f getTextPoint(Graphics g, String text, double cX, double cY) {
        FontMetrics metrics = g.getFontMetrics();
        Rectangle2D stringBounds = metrics.getStringBounds(text, g);
        double halfWidth = stringBounds.getWidth() / 2;
        double halfHeight = stringBounds.getHeight() / 2;
        double ascent = metrics.getAscent();
        float x = (float) (cX - halfWidth);
        float y = (float) (cY - halfHeight + ascent);
        return new Point2f(x, y);
    }

    /**
     *
     * @param g
     * @param text
     * @return
     */
    public Rectangle2D getTextBounds(Graphics g, String text) {
        FontMetrics fontMetrics;
        fontMetrics = g.getFontMetrics();
        return fontMetrics.getStringBounds(text, g);
    }

    /**
     *
     * @param ac
     * @param x
     * @param y
     * @param boundsTree
     */
    public void translateTo(IAtomContainer ac, double x, double y, BoundsTree boundsTree) {
        Rectangle2D bounds = boundsTree.getRoot();
        double dx = x - bounds.getCenterX();
        double dy = y - bounds.getCenterY();
        for (IAtom atom : ac.atoms()) {
            atom.getPoint2d().x += dx;
            atom.getPoint2d().y += dy;
        }
        boundsTree.shift(dx, dy);
    }
}
