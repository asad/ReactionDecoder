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

import java.awt.BasicStroke;
import static java.awt.Color.BLACK;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import static uk.ac.ebi.reactionblast.graphics.direct.Params.ArrowType.BIDIRECTIONAL;

/**
 *
 * @author asad
 */
public class DirectArrowDrawer extends AbstractDirectDrawer {

    private final static Vector2d X_AXIS = new Vector2d(1, 0);
    private final static Vector2d Y_AXIS = new Vector2d(0, 1);

    private final Params params;

    /**
     *
     * @param params
     */
    public DirectArrowDrawer(Params params) {
        this.params = params;
    }

    /**
     *
     * @param g
     * @param c
     * @param v
     */
    public void drawArrow(Graphics2D g, Point2d c, Vector2d v) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(new BasicStroke());
        if (params.drawFatArrow) {
            if (params.arrowType == BIDIRECTIONAL) {
                drawDoubleHeadedFatArrow(g, c, v, null);
            } else {
                drawFatArrow(g, c, v, null);
            }
        } else {
            drawThinArrow(g, c, v, null);
        }
        g.setStroke(savedStroke);
    }

    /**
     *
     * @param g
     * @param c
     * @param v
     * @param text
     */
    public void drawFatArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = params.arrowLength;
        int arrowHeadLength = params.arrowHeadLength;
        int arrowHeadIndent = params.arrowHeadIndent;
        int arrowBodyWidth = params.arrowBodyWidth;

        double arrowHeadAngleRad = toRadians(params.arrowHeadAngle);
        double arrowHeadAngleRadPrime = toRadians(360 - params.arrowHeadAngle);
        double cosA = cos(arrowHeadAngleRad);
        double sinA = sin(arrowHeadAngleRad);
        double cosAPrime = cos(arrowHeadAngleRadPrime);
        double sinAPrime = sin(arrowHeadAngleRadPrime);

        int halfLength = arrowLength / 2;
        int halfBodyWidth = arrowBodyWidth / 2;

        g.setColor(BLACK);

        Vector2d nV = new Vector2d(v);
        nV.negate();
        Vector2d p = new Vector2d(v.y, -v.x);
        Vector2d nP = new Vector2d(-v.y, v.x);

        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd(halfLength, nV, c);

        Point2d upperTail = new Point2d(tail);
        upperTail.scaleAdd(halfBodyWidth, nP, upperTail);

        Point2d lowerTail = new Point2d(tail);
        lowerTail.scaleAdd(halfBodyWidth, p, lowerTail);

        Point2d head = new Point2d(c);
        head.scaleAdd(halfLength, v, c);

        // CW = clockwise (around the head point) CCW = counterclockwise
        Vector2d ccwVec = new Vector2d((cosA * nV.x) + (sinA * nV.y), (cosA * nV.y) - (sinA * nV.x));
        Vector2d cwVec = new Vector2d((cosAPrime * nV.x) + (sinAPrime * nV.y), (cosAPrime * nV.y) - (sinAPrime * nV.x));

        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd(arrowHeadLength, ccwVec, head);

        Point2d indentCCW = new Point2d(headCCW);
        indentCCW.scaleAdd(arrowHeadIndent, p, indentCCW);

        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd(arrowHeadLength, cwVec, head);

        Point2d indentCW = new Point2d(headCW);
        indentCW.scaleAdd(arrowHeadIndent, nP, indentCW);

        Path2D polygon = new Path2D.Double();
        polygon.moveTo(head.x, head.y);
        polygon.lineTo(headCCW.x, headCCW.y);
        polygon.lineTo(indentCCW.x, indentCCW.y);
        polygon.lineTo(upperTail.x, upperTail.y);
        polygon.lineTo(lowerTail.x, lowerTail.y);
        polygon.lineTo(indentCW.x, indentCW.y);
        polygon.lineTo(headCW.x, headCW.y);
        polygon.closePath();

        if (params.drawArrowFilled) {
            g.fill(polygon);
        } else {
            g.draw(polygon);
        }

        if (text != null) {
            drawText(g, text, c, v, nV);
        }
    }

    /**
     *
     * @param g
     * @param c
     * @param v
     * @param text
     */
    public void drawDoubleHeadedFatArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = params.arrowLength;
        int arrowHeadLength = params.arrowHeadLength;
        int arrowHeadIndent = params.arrowHeadIndent;

        double arrowHeadAngleRad = toRadians(params.arrowHeadAngle);
        double arrowHeadAngleRadPrime = toRadians(360 - params.arrowHeadAngle);
        double cosA = cos(arrowHeadAngleRad);
        double sinA = sin(arrowHeadAngleRad);
        double cosAPrime = cos(arrowHeadAngleRadPrime);
        double sinAPrime = sin(arrowHeadAngleRadPrime);

        int halfLength = arrowLength / 2;

        g.setColor(BLACK);

        Vector2d nV = new Vector2d(v);
        nV.negate();
        Vector2d p = new Vector2d(v.y, -v.x);
        Vector2d nP = new Vector2d(-v.y, v.x);

        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd(halfLength, nV, c);

        Point2d head = new Point2d(c);
        head.scaleAdd(halfLength, v, c);

        // CW = clockwise (around the head point) CCW = counterclockwise
        Vector2d ccwVec = new Vector2d((cosA * nV.x) + (sinA * nV.y), (cosA * nV.y) - (sinA * nV.x));
        Vector2d nCCWVec = new Vector2d(ccwVec);
        nCCWVec.negate();

        Vector2d cwVec = new Vector2d((cosAPrime * nV.x) + (sinAPrime * nV.y), (cosAPrime * nV.y) - (sinAPrime * nV.x));
        Vector2d nCWVec = new Vector2d(cwVec);
        nCWVec.negate();

        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd(arrowHeadLength, ccwVec, head);

        Point2d headIndentCCW = new Point2d(headCCW);
        headIndentCCW.scaleAdd(arrowHeadIndent, p, headIndentCCW);

        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd(arrowHeadLength, cwVec, head);

        Point2d headIndentCW = new Point2d(headCW);
        headIndentCW.scaleAdd(arrowHeadIndent, nP, headIndentCW);

        Point2d tailCCW = new Point2d(tail);
        tailCCW.scaleAdd(arrowHeadLength, nCWVec, tailCCW);

        Point2d tailCW = new Point2d(tail);
        tailCW.scaleAdd(arrowHeadLength, nCCWVec, tailCW);

        Point2d upperTail = new Point2d(tailCCW);
        upperTail.scaleAdd(arrowHeadIndent, p, upperTail);

        Point2d lowerTail = new Point2d(tailCW);
        lowerTail.scaleAdd(arrowHeadIndent, nP, lowerTail);

        Path2D polygon = new Path2D.Double();
        polygon.moveTo(head.x, head.y);
        polygon.lineTo(headCCW.x, headCCW.y);
        polygon.lineTo(headIndentCCW.x, headIndentCCW.y);
        polygon.lineTo(upperTail.x, upperTail.y);
        polygon.lineTo(tailCCW.x, tailCCW.y);
        polygon.lineTo(tail.x, tail.y);
        polygon.lineTo(tailCW.x, tailCW.y);
        polygon.lineTo(lowerTail.x, lowerTail.y);
        polygon.lineTo(headIndentCW.x, headIndentCW.y);
        polygon.lineTo(headCW.x, headCW.y);
        polygon.closePath();

        if (params.drawArrowFilled) {
            g.fill(polygon);
        } else {
            g.draw(polygon);
        }

        if (text != null) {
            drawText(g, text, c, v, nV);
        }

    }

    /**
     *
     * @param g
     * @param c
     * @param v
     * @param text
     */
    public void drawThinArrow(Graphics2D g, Point2d c, Vector2d v, String text) {
        int arrowLength = params.arrowLength;
        int arrowHeadLength = params.arrowHeadLength;
        double arrowHeadAngleRad = toRadians(params.arrowHeadAngle);
        double arrowHeadAngleRadPrime = toRadians(360 - params.arrowHeadAngle);
        double cosA = cos(arrowHeadAngleRad);
        double sinA = sin(arrowHeadAngleRad);
        double cosAPrime = cos(arrowHeadAngleRadPrime);
        double sinAPrime = sin(arrowHeadAngleRadPrime);

        int halfLength = arrowLength / 2;

        g.setColor(BLACK);

        Vector2d nV = new Vector2d(v);
        nV.negate();

        Point2d tail = new Point2d(c.x, c.y);
        tail.scaleAdd(halfLength, nV, c);
        Point2d head = new Point2d(c);
        head.scaleAdd(halfLength, v, c);

        // CW = clockwise (around the head point) CCW = counterclockwise
        Vector2d ccwVec = new Vector2d((cosA * nV.x) + (sinA * nV.y), (cosA * nV.y) - (sinA * nV.x));
        Vector2d cwVec = new Vector2d((cosAPrime * nV.x) + (sinAPrime * nV.y), (cosAPrime * nV.y) - (sinAPrime * nV.x));

        Point2d headCCW = new Point2d(head.x, head.y);
        headCCW.scaleAdd(arrowHeadLength, ccwVec, head);
        Point2d headCW = new Point2d(head.x, head.y);
        headCW.scaleAdd(arrowHeadLength, cwVec, head);

        g.draw(new Line2D.Double(tail.x, tail.y, head.x, head.y));
        g.draw(new Line2D.Double(head.x, head.y, headCCW.x, headCCW.y));
        g.draw(new Line2D.Double(head.x, head.y, headCW.x, headCW.y));

        if (text != null) {
            drawText(g, text, c, v, nV);
        }
    }

    private void drawText(Graphics2D g, String text, Point2d c, Vector2d v, Vector2d nV) {
        AffineTransform originalTransform = g.getTransform();
        double angle = getAngle(v);

        Rectangle2D textBounds = getTextBounds(g, text);
        double distance = textBounds.getWidth() / 2;
        Point2d start = new Point2d(c.x, c.y);
        if (angle < toRadians(90) || angle > toRadians(270)) {
            start.scaleAdd(distance, nV, c);
        } else {
            start.scaleAdd(distance, v, c);
            double angDeg = (180 + toDegrees(angle)) % 360;
            angle = toRadians(angDeg);
        }

        g.translate(start.x, start.y);
        g.rotate(angle);

        Font font = g.getFont();
        FontRenderContext frc = g.getFontRenderContext();
        GlyphVector gv = font.createGlyphVector(frc, text);
        int length = gv.getNumGlyphs();
        for (int i = 0; i < length; i++) {
            g.fill(gv.getGlyphOutline(i));
        }

        g.rotate((2 * PI) - angle);
        g.setTransform(originalTransform);
    }

    private double getAngle(Vector2d v) {
        double xAngle = X_AXIS.angle(v);
        double yAngle = Y_AXIS.angle(v);
        if (xAngle < toRadians(90)) {          // Q1 or Q2
            if (yAngle < toRadians(90)) {      // Q2
                return xAngle;
            } else {                                // Q1
                return toRadians(360) - xAngle;
            }
        } else {                                    // Q3 or Q4
            return toRadians(90) + yAngle;
//            if (yAngle < Math.toRadians(90)) {      // Q3
//                
//            } else {                                // Q4
//                return Math.toRadians(90) + yAngle;
//            }
        }
    }

}
