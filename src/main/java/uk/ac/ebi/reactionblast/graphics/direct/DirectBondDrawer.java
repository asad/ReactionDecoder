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
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import static java.lang.Math.min;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import org.openscience.cdk.exception.Intractable;
import static org.openscience.cdk.geometry.GeometryTools.get2DCenter;
import static org.openscience.cdk.geometry.GeometryTools.getRectangle2D;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import org.openscience.cdk.interfaces.IBond.Stereo;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN_INVERTED;
import static org.openscience.cdk.interfaces.IBond.Stereo.NONE;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_INVERTED;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_OR_DOWN;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import uk.ac.ebi.reactionblast.graphics.direct.Params.BondStrokeCap;
import uk.ac.ebi.reactionblast.graphics.direct.Params.BondStrokeJoin;

/**
 *
 * @author asad
 */
public class DirectBondDrawer extends AbstractDirectDrawer {

    private final LabelManager labelManager;
    private final Stroke dashedWedgeStroke;
    private Stroke bondStroke;

    /**
     *
     * @param params
     * @param labelManager
     */
    public DirectBondDrawer(Params params, LabelManager labelManager) {
        setParams(params);
        this.labelManager = labelManager;
        dashedWedgeStroke = new BasicStroke(params.dashedWedgeStroke);
    }

    private void setBondStroke() {
        int cap = BasicStroke.CAP_BUTT;
        int join = BasicStroke.JOIN_BEVEL;
        if (params.bondStrokeCap == BondStrokeCap.ROUND) {
            cap = BasicStroke.CAP_ROUND;
        } else if (params.bondStrokeCap == BondStrokeCap.SQUARE) {
            cap = BasicStroke.CAP_SQUARE;
        }

        if (params.bondStrokeJoin == BondStrokeJoin.BEVEL) {
            join = BasicStroke.JOIN_BEVEL;
        } else if (params.bondStrokeJoin == BondStrokeJoin.ROUND) {
            join = BasicStroke.JOIN_ROUND;
        }
        bondStroke = new BasicStroke(params.bondStrokeWidth, cap, join);
    }

    /**
     *
     * @param molecule
     * @param g
     * @throws org.openscience.cdk.exception.Intractable
     */
    public void drawBonds(IAtomContainer molecule, Graphics2D g) throws Intractable {
        setBondStroke();
        g.setStroke(bondStroke);

        //IRingSet ringSet = new SSSRFinder(molecule).findSSSR();
        //New Method
        CycleFinder cf = Cycles.mcb();
        Cycles cycles = cf.find(molecule); // ignore error - essential cycles do not check tractability
        IRingSet ringSet = cycles.toRingSet();

        ringSet.sortAtomContainers(new AtomContainerComparatorBy2DCenter());
        addRingCentersToAtomAnnotationPositions(molecule, ringSet);
        Map<IBond, IAtomContainer> bondRingMap = fillBondRingMap(ringSet);

        g.setColor(BLACK);
        for (IBond bond : molecule.bonds()) {
            if (shouldDraw(bond)) {
                drawBond(bond, bondRingMap, g);
            }
            labelManager.addBondToAtomAnnotationPositions(bond);
        }

        if (params.drawAromaticCircles) {
            for (IAtomContainer ring : ringSet.atomContainers()) {
                if (ringIsAromatic(ring)) {
                    drawRingCircle(ring, g);
                }
            }
        }
        List<IBond> drawnRingBonds = new ArrayList<>();
        for (IAtomContainer ring : ringSet.atomContainers()) {
            Point2d c = get2DCenter(ring);
            for (IBond bond : ring.bonds()) {
                if (drawnRingBonds.contains(bond)) {
                } else if (bond.getFlag(ISAROMATIC) && params.drawAromaticCircles) {
                    Point2d p1 = bond.getAtom(0).getPoint2d();
                    Point2d p2 = bond.getAtom(1).getPoint2d();
                    drawOffsetBond(p1, p2, c, g);
                    drawnRingBonds.add(bond);
                } else if (bond.getOrder() == SINGLE) {
                } else {
                    Point2d p1 = bond.getAtom(0).getPoint2d();
                    Point2d p2 = bond.getAtom(1).getPoint2d();
                    drawOffsetBond(p1, p2, c, g);
                    drawnRingBonds.add(bond);
                }
            }
        }
    }

    private void addRingCentersToAtomAnnotationPositions(IAtomContainer mol, IRingSet ringSet) {
        for (IAtomContainer ring : ringSet.atomContainers()) {
            for (IAtom atom : ring.atoms()) {
                List<IAtom> connectedAtoms = mol.getConnectedAtomsList(atom);
                List<IAtom> connectedAtomsInRing = new ArrayList<>();
                for (IAtom connectedAtom : connectedAtoms) {
                    if (ring.contains(connectedAtom)) {
                        connectedAtomsInRing.add(connectedAtom);
                    }
                }
                labelManager.addRingCenterToAtomAnnotationPosition(
                        atom, connectedAtomsInRing);
            }
        }
    }

    private Map<IBond, IAtomContainer> fillBondRingMap(IRingSet ringSet) {
        Map<IBond, IAtomContainer> bondRingMap
                = new HashMap<>();
        for (IAtomContainer ringAsAtomContainer : ringSet.atomContainers()) {
            IRing ring = (IRing) ringAsAtomContainer;
            for (IBond bond : ring.bonds()) {
                bondRingMap.put(bond, ring);
            }
        }
        return bondRingMap;
    }

    /**
     *
     * @param bond
     * @param bondRingMap
     * @param g
     */
    public void drawBond(
            IBond bond, Map<IBond, IAtomContainer> bondRingMap, Graphics2D g) {
        Point2d p1 = bond.getAtom(0).getPoint2d();
        Point2d p2 = bond.getAtom(1).getPoint2d();
        IBond.Order order = bond.getOrder();
        IBond.Stereo stereo = bond.getStereo();
        if (stereo == NONE
                && (order == SINGLE || bond.getFlag(ISAROMATIC))) {
            drawLine(p1, p2, g);
        } else if (order == DOUBLE) {
            if (bondRingMap.containsKey(bond)) {
                drawLine(p1, p2, g);
            } else {
                drawDoubleBond(p1, p2, g);
            }
        } else if (order == TRIPLE) {
            drawTripleBond(p1, p2, g);
        } else if (stereo != NONE) {
            drawStereo(p1, p2, stereo, g);
        }
    }

    private void drawTripleBond(Point2d p1, Point2d p2, Graphics2D g) {
        Vector2d perpendicular = makePerpendicular(p1, p2);
        perpendicular.scale(params.tripleBondGap);
        Vector2d negativePerp = new Vector2d(perpendicular);
        negativePerp.negate();

        drawLine(displace(p1, perpendicular), displace(p2, perpendicular), g);
        drawLine(p1, p2, g);
        drawLine(displace(p1, negativePerp), displace(p2, negativePerp), g);
    }

    private void drawStereo(Point2d p1, Point2d p2, Stereo stereo, Graphics2D g) {
        if (null != stereo) {
            switch (stereo) {
                case UP_OR_DOWN:
                    drawWigglyLine(p1, p2, g);
                    break;
                case DOWN:
                    drawWedge(p1, p2, false, g);
                    break;
                case DOWN_INVERTED:
                    drawWedge(p2, p1, false, g);
                    break;
                case UP:
                    drawWedge(p1, p2, true, g);
                    break;
                case UP_INVERTED:
                    drawWedge(p2, p1, true, g);
                    break;
                // ?
                default:
                    break;
            }
        }
    }

    private void drawWedge(Point2d p1, Point2d p2, boolean isFilled, Graphics2D g) {
        Vector2d halfWidthVector = new Vector2d(p2.y - p1.y, p1.x - p2.x);
        halfWidthVector.normalize();
        halfWidthVector.scale(params.filledWedgeWidth / 2);
        Vector2d negHalfWidthVector = new Vector2d(halfWidthVector);
        negHalfWidthVector.negate();
        Point2d p2a = displace(p2, halfWidthVector);
        Point2d p2b = displace(p2, negHalfWidthVector);

        if (isFilled) {
            drawFilledWedge(p1, p2a, p2b, g);
        } else {
            drawDashedWedge2(p1, p2a, p2b, g);
        }

    }

    /**
     *
     * @param a
     * @param b
     * @param c
     * @param g
     */
    public void drawDashedWedge(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(dashedWedgeStroke);
        double distance = b.distance(a);
        double gapFactor = params.dashedGapFactor;
        double gap = distance * gapFactor;
        double numberOfDashes = distance / gap;
        double d = 0;

        // draw by interpolating along the edges of the triangle
        for (int i = 0; i < numberOfDashes; i++) {
            Point2d p1 = new Point2d();
            p1.interpolate(a, b, d);
            Point2d p2 = new Point2d();
            p2.interpolate(a, c, d);

            drawLine(p1, p2, g);
            if (distance * (d + gapFactor) >= distance) {
                break;
            } else {
                d += gapFactor;
            }
        }
        g.setStroke(savedStroke);
    }

    /**
     *
     * @param a
     * @param b
     * @param c
     * @param g
     */
    public void drawDashedWedge2(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Stroke savedStroke = g.getStroke();
        g.setStroke(dashedWedgeStroke);
        double distance = b.distance(a);
        double gapFactor = params.dashedGapFactor;
        double gap = distance * gapFactor;
        double numberOfDashes = distance / gap;
        double currentDistance = 0;
        Point2d d = new Point2d(b);
        d.interpolate(c, 0.5);
        Vector2d perp = makePerpendicular(a, d);
        Vector2d nPerp = new Vector2d(perp);
        nPerp.negate();
        double maxWidth = params.dashedWedgeWidth / 4;
        double currentWidth = maxWidth * params.dashedWidthFactor;
        // draw like a ladder with increasing rung length
        for (int i = 0; i < numberOfDashes; i++) {
            Point2d rungCenter = new Point2d(a);
            rungCenter.interpolate(d, currentDistance);

            Point2d p1 = new Point2d(rungCenter);
            p1.scaleAdd(currentWidth, perp, p1);

            Point2d p2 = new Point2d(rungCenter);
            p2.scaleAdd(currentWidth, nPerp, p2);

            drawLine(p1, p2, g);
            if (distance * (currentDistance + gapFactor) >= distance) {
                break;
            } else {
                currentDistance += gapFactor;
                currentWidth += maxWidth * (params.dashedWidthFactor);
            }
        }
        g.setStroke(savedStroke);
    }

    private void drawFilledWedge(Point2d a, Point2d b, Point2d c, Graphics2D g) {
        Path2D path = new Path2D.Double();
        path.moveTo(a.x, a.y);
        path.lineTo(b.x, b.y);
        path.lineTo(c.x, c.y);
        path.closePath();
        g.fill(path);
    }

    /**
     *
     * @param p1
     * @param p2
     * @param g
     */
    public void drawWigglyLine(Point2d p1, Point2d p2, Graphics2D g) {
        double gapProportion = 0.1;
        double wiggleWidth = params.wiggleLineWidth;

        Vector2d line = new Vector2d(p2);
        line.sub(p1);
        double length = line.length();
        double gap = length * gapProportion;
        int numberOfSegments = 10;

        line.normalize();
        Vector2d perpendicular = makePerpendicular(line);
        Vector2d negPerp = new Vector2d(perpendicular);
        negPerp.negate();
        Point2d centerLinePoint = new Point2d(p1);

        Path2D path = new Path2D.Double();
        path.moveTo(p1.x, p1.y);

        // start at the first peak
        centerLinePoint.scaleAdd(gap / 2, line, centerLinePoint);
        Point2d tipPoint = new Point2d(centerLinePoint);
        tipPoint.scaleAdd(wiggleWidth / 2, perpendicular, tipPoint);
        for (int i = 0; i < numberOfSegments - 1; i++) {
            centerLinePoint.scaleAdd(gap / 2, line, centerLinePoint);

            path.quadTo(
                    tipPoint.x, tipPoint.y, centerLinePoint.x, centerLinePoint.y);
            centerLinePoint.scaleAdd(gap / 2, line, centerLinePoint);

            // alternate between up and down
            tipPoint = new Point2d(centerLinePoint);
            if (i % 2 == 0) {
                tipPoint.scaleAdd(wiggleWidth / 2, negPerp, tipPoint);
            } else {
                tipPoint.scaleAdd(wiggleWidth / 2, perpendicular, tipPoint);
            }
        }
        // finish the last curve

        g.draw(path);
    }

    private Vector2d makePerpendicular(Point2d p1, Point2d p2) {
        Vector2d line = new Vector2d(p1);
        line.sub(p2);
        line.normalize();
        return makePerpendicular(line);
    }

    private void drawDoubleBond(Point2d p1, Point2d p2, Graphics2D g) {
        Vector2d perpendicular = makePerpendicular(p1, p2);
        perpendicular.scale(params.doubleBondGap);
        Vector2d negativePerp = new Vector2d(perpendicular);
        negativePerp.negate();

        drawLine(displace(p1, perpendicular), displace(p2, perpendicular), g);
        drawLine(displace(p1, negativePerp), displace(p2, negativePerp), g);
    }

    private void drawOffsetBond(Point2d p1, Point2d p2, Point2d c, Graphics2D g) {
        double distanceProportion = params.offsetBondDistanceProportion;
        Point2d w = new Point2d();
        w.interpolate(c, p1, distanceProportion);

        Point2d u = new Point2d();
        u.interpolate(c, p2, distanceProportion);

        drawLine(w, u, g);
    }

    private Vector2d makePerpendicular(Vector2d line) {
        Vector2d perp = new Vector2d(-line.y, line.x);
        perp.normalize();
        return perp;
    }

    private Point2d displace(Point2d point, Vector2d vector) {
        Point2d displacedPoint = new Point2d(point);
        displacedPoint.add(vector);
        return displacedPoint;
    }

    private void drawRingCircle(IAtomContainer ring, Graphics2D g) {
        Point2d center = get2DCenter(ring);
        Rectangle2D bounds = getRectangle2D(ring);
        double diameter = min(bounds.getWidth(), bounds.getHeight());
        diameter *= params.ringProportion;
        double radius = diameter / 2;
        g.draw(new Ellipse2D.Double(center.x - radius, center.y - radius, diameter, diameter));
    }

    private boolean ringIsAromatic(IAtomContainer ring) {
        for (IAtom atom : ring.atoms()) {
            if (!atom.getFlag(ISAROMATIC)) {
                return false;
            }
        }
        for (IBond b : ring.bonds()) {
            if (!b.getFlag(ISAROMATIC)) {
                return false;
            }
        }
        return true;
    }

    private boolean shouldDraw(IBond bond) {
        boolean symbol0IsH = bond.getAtom(0).getSymbol().equals("H");
        boolean symbol1IsH = bond.getAtom(1).getSymbol().equals("H");
        boolean bothAreH = symbol0IsH && symbol1IsH;
        boolean atLeastOneIsH = symbol0IsH || symbol1IsH;
        boolean neitherAreH = !symbol0IsH && !symbol1IsH;
        if (bothAreH || neitherAreH) {
            return true;
        } else if (atLeastOneIsH) {
            return params.drawExplicitHydrogens;
        } else {
            return true;
        }
    }
}
