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

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import static org.openscience.cdk.geometry.GeometryTools.get2DCenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.ConvexHull;

/**
 * Aligns molecules to a line specified by a vector. There are different methods
 * for different scenarios - for example, alignToMaxWidth will use the line
 * through the two atoms with greatest distance, which tends to align a molecule
 * close to its actual axis.
 *
 * The alignToMinAreaBox, on the other hand, tries to find the rectangle with
 * the smallest area that still encloses all the atoms, and aligns along the
 * longest axis of this box. This often orients long molecules across the
 * diagonal of the area.
 *
 * @author maclean
 *
 */
public class MoleculeAligner {

    /**
     *
     */
    public static final Vector2d X_AXIS = new Vector2d(1, 0);

    /**
     *
     */
    public static final Vector2d Y_AXIS = new Vector2d(0, 1);

    /**
     * Finds the minimum-area bounding box of the atom container and uses the
     * longest side as the central axis. Note that this will tend to align long
     * molecules (such as hydrocarbon chains) at 45&deg; as this will be the
     * best solution for the minimum area box.
     *
     * @param atomContainer the atoms to align
     * @param axis the axis to align to
     */
    public static void alignToMinAreaBox(IAtomContainer atomContainer, Vector2d axis) {
        ConvexHull hull = new ConvexHull(atomContainer);
        alignToAxis(atomContainer, hull.getMajorAxis(), axis, hull.getCenter());
    }

    /**
     *
     * @param atomContainer
     * @return
     */
    public static Vector2d getMaxWidthVector(IAtomContainer atomContainer) {
        int nAtoms = atomContainer.getAtomCount();
        Vector2d widthVector = null;
        IAtom maxI = null;
        IAtom maxJ = null;
        double maxDistance = 0;
        for (int indexI = nAtoms - 1; indexI >= 0; indexI--) {
            IAtom atomI = atomContainer.getAtom(indexI);
            Point2d pointI = atomI.getPoint2d();
            if (pointI == null) {
                continue;
            }
            for (int indexJ = indexI - 1; indexJ >= 0; indexJ--) {
                IAtom atomJ = atomContainer.getAtom(indexJ);
                Point2d pointJ = atomJ.getPoint2d();
                if (pointJ == null) {
                    continue;
                }
                double distance = pointI.distance(pointJ);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    maxI = atomI;
                    maxJ = atomJ;
                }
            }
        }
//        System.out.println("maxI = " + atomContainer.indexOf(maxI)
//                         + "maxJ = " + atomContainer.indexOf(maxJ));
        if (maxI != null && maxJ != null) {
            widthVector = new Vector2d(maxI.getPoint2d());
            widthVector.sub(maxJ.getPoint2d());
        } else {
            return new Vector2d(0, 0);
        }
        return widthVector;
    }

    /**
     *
     * @param atomContainer
     * @param axis
     */
    public static void alignToMaxWidth(IAtomContainer atomContainer, Vector2d axis) {
        Vector2d widthVector = getMaxWidthVector(atomContainer);
        Point2d center = get2DCenter(atomContainer);
        alignToAxis(atomContainer, widthVector, axis, center);
    }

    private static double getPolarAngle(Vector2d vector) {
        double x = vector.x;
        double y = vector.y;
        if (x > 0) {
            return atan(y / x);
        } else if (x < 0) {
            if (y >= 0) {
                return atan(y / x) + PI;
            } else {
                return atan(y / x) - PI;
            }
        } else if (y > 0) {
            return PI / 2;
        } else if (y < 0) {
            return -(PI / 2);
        } else {
            return 0;
        }
    }

    /**
     * Given two axes as vectors, calculate the minimum angle needed to rotate
     * one onto the other. The <code>axisFrom</code> vector is treated as a line
     * - so both directions along it are considered - while the <code>
     * axisTo</code> vector is treated as a vector.
     *
     * @param axisFrom
     * @param axisTo
     * @return
     */
    public static double getMinAngle(Vector2d axisFrom, Vector2d axisTo) {
        // all angles converted to [0, 2*PI] from [-PI, PI]
        double polarAngleForwardFrom = atan2(axisFrom.y, axisFrom.x);
        double polarAngleBackwardFrom = atan2(-axisFrom.y, -axisFrom.x);
        double polarAngleTo = atan2(axisTo.y, axisTo.x);
        double forwardDiff = polarAngleForwardFrom - polarAngleTo;
        double backwardDiff = polarAngleBackwardFrom - polarAngleTo;
        double minAngleDiff;
        if (abs(forwardDiff) < abs(backwardDiff)) {
            minAngleDiff = forwardDiff;
        } else {
            minAngleDiff = backwardDiff;
        }

//        System.out.println(
//                "Axis from " + String.format("(%2.2f %2.2f)", axisFrom.x, axisFrom.y)
//                + "Axis to " + axisTo
//                + " angleForwardFrom "  + toStr(polarAngleForwardFrom)
//                + " angleBackwardFrom " + toStr(polarAngleBackwardFrom)
//                + " angleTo "           + toStr(polarAngleTo)
//                + " forwardDiff "       + toStr(forwardDiff)
//                + " backwardDiff "      + toStr(backwardDiff)
//                + " minAngleDiff "      + toStr(minAngleDiff)
//                );
        return -minAngleDiff;
    }

//    private static String toStr(double angle) {
//        return String.format("%2.2f", Math.toDegrees(angle));
//    }
    /**
     *
     * @param atomContainer
     * @param axisFrom
     * @param axisTo
     * @param center
     */
    public static void alignToAxis(IAtomContainer atomContainer,
            Vector2d axisFrom,
            Vector2d axisTo,
            Point2d center) {
        double angle = getMinAngle(axisFrom, axisTo);
        double cosA = cos(angle);
        double sinA = sin(angle);
        double minCosA = 1 - cosA;
        for (IAtom atom : atomContainer.atoms()) {
            Point2d p = atom.getPoint2d();
            double x = (cosA * p.x) - (sinA * p.y)
                    + (center.x * minCosA) + (center.y * sinA);
            double y = (sinA * p.x) + (cosA * p.y)
                    + (center.y * minCosA) - (center.x * sinA);
            p.x = x;
            p.y = y;
            atom.setPoint2d(p);
        }
    }

    private MoleculeAligner() {
    }
}
