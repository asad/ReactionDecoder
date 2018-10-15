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

import java.awt.geom.Rectangle2D;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.min;
import static java.lang.Math.sin;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.asList;
import static java.util.Arrays.sort;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author asad
 */
public class ConvexHull implements Iterable<Point2d> {

    private Point2d[] hull;

    private String[] hullIDs;
    private final Vector2d X_AXIS = new Vector2d(1, 0);

    /**
     *
     * @param atomContainer
     */
    public ConvexHull(IAtomContainer atomContainer) {
        Point2d[] points = new Point2d[atomContainer.getAtomCount()];
        int i = 0;
        for (IAtom atom : atomContainer.atoms()) {
            points[i] = atom.getPoint2d();
            i++;
        }
        if (i < atomContainer.getAtomCount()) {
            Point2d[] nonNullPoints = new Point2d[i];
            int k = 0;
            for (Point2d point : points) {
                if (point != null) {
                    nonNullPoints[k] = point;
                    k++;
                }
            }
            points = nonNullPoints;
        }
        makeFromPoints(points);
        hullIDs = new String[hull.length];
        for (IAtom atom : atomContainer.atoms()) {
            if (atom.getPoint2d() != null && atom.getID() != null) {
                Point2d point = atom.getPoint2d();
                String id = atom.getID();
                int hullIndex = 0;
                for (Point2d hullPoint : hull) {
                    if (hullPoint == point) {
                        hullIDs[hullIndex] = id;
                        break;
                    }
                    hullIndex++;
                }
            }
        }
    }

    /**
     *
     * @param points
     */
    public ConvexHull(Point2d[] points) {
        makeFromPoints(points);
    }

    /**
     *
     * @return
     */
    public Vector2d getMajorAxis() {
        Rectangle minimumAreaBoundingRectangle
                = //            getMinimumAreaBoundingRectangle();
                getMinimumAreaBoundingRectangleBruteForce();
        return minimumAreaBoundingRectangle.getMajorAxis();
    }

    /**
     *
     * @return
     */
    public Point2d getCenter() {
        Point2d center = new Point2d();
        for (Point2d hullPoint : hull) {
            center.x += hullPoint.x;
            center.y += hullPoint.y;
        }
        center.x /= hull.length;
        center.y /= hull.length;
        return center;
    }

    /**
     *
     * @return
     */
    public Rectangle getMinimumAreaBoundingRectangleBruteForce() {
        Rectangle minRect = null;
        double minArea = -1;
        int winnerIndex = -1;
        for (int index = 0; index < hull.length - 1; index++) {
            Vector2d edge = edgeVector(hull[index], hull[index + 1]);
            Rectangle rect = getRectangleBrute(edge, index, index + 1);
            double area = rect.area();
            if (minRect == null || area < minArea) {
                minRect = rect;
                minArea = area;
                winnerIndex = index;
            }
//            System.out.println("rect " + rect);
        }
        Vector2d edge = edgeVector(hull[hull.length - 1], hull[0]);
        Rectangle rect = getRectangleBrute(edge, hull.length - 1, 0);
        double area = rect.area();
        if (minRect == null || area < minArea) {
            minRect = rect;
            minArea = area;
            winnerIndex = hull.length;
        }
//        System.out.println("winner = " + winnerIndex);

        return minRect;
    }

    /**
     *
     * @return
     */
    public Rectangle getMinimumAreaBoundingRectangle() {
        assert hull != null;
        Point2d minY = null;
        Point2d maxY = null;
        int indexA = -1;
        int indexB = -1;
        for (int index = 0; index < hull.length; index++) {
            Point2d point = hull[index];
            if (minY == null || point.y < minY.y) {
                minY = point;
                indexA = index;
            }
            if (maxY == null || point.y > maxY.y) {
                maxY = point;
                indexB = index;
            }
        }
        Vector2d caliperA = new Vector2d(1, 0);
        Vector2d caliperB = new Vector2d(-1, 0);
        double rotatedAngle = 0;
        double minArea = MAX_VALUE;
        Rectangle minRect = null;
        while (rotatedAngle < PI) {
            if (indexA == hull.length - 1) {
                indexA = 0;
            }
            if (indexB == hull.length - 1) {
                indexB = 0;
            }
            Vector2d edgeA = edgeVector(hull[indexA], hull[indexA + 1]);
            Vector2d edgeB = edgeVector(hull[indexB], hull[indexB + 1]);
            double angleA = edgeA.angle(caliperA);
            double angleB = edgeB.angle(caliperB);
            double minAngle = min(angleA, angleB);
            caliperA = rotate(caliperA, minAngle);
            caliperB = rotate(caliperB, minAngle);
            Rectangle rectangle;
            if (angleA < angleB) {
                indexA++;
                rectangle = getRectangle(edgeA, indexA, indexA + 1);
            } else {
                indexB++;
                rectangle = getRectangle(edgeB, indexB, indexB + 1);
            }
            rotatedAngle += minAngle;
            double area = rectangle.area();
            if (area < minArea) {
                minArea = area;
                minRect = rectangle;
            }
//            System.out.println(
//                    "rotated angle = " + rotatedAngle
//                    + " min " + minArea
//                    + " r = " + rectangle
//                    + " caliper_a " + caliperA
//                    + " caliper_b " + caliperB);
        }
        return minRect;
    }

    private Rectangle getRectangleBrute(
            Vector2d vector, int tailPointIndex, int headPointIndex) {
        Point2d headPoint = hull[headPointIndex];
        Point2d tailPoint = hull[tailPointIndex];
        int index = headPointIndex;
        int visited = 0;
        int min = 0;
        Point2d vMax = null;
        Point2d thirdPoint = null;
        double thirdPointDist = 0.0;
        double minAngle = PI * 2;
        double maxAngle = 0;
        Vector2d vN = new Vector2d(vector);
        vN.normalize();
//        System.out.println("tailIndex " + tailPointIndex 
//                + " visiting points from " 
//                + hullIDs[headPointIndex]
//                + " tailPt" + toString(tailPoint)
//                + " headPt" + toString(headPoint));
        int max = 0;
        while (visited < hull.length) {
            if (index == hull.length) {
                index = 0;
            }
            if (vMax == null) {
                vMax = hull[index];
            } else {
                double angle = prj(tailPoint, headPoint, hull[index]);
//                System.out.println(index + " proj " + hullIDs[index] + " " + angle);
                if (angle < minAngle) {
                    min = index;
                    minAngle = angle;
                }
                if (angle > maxAngle) {
                    vMax = hull[index];
                    max = index;
                    maxAngle = angle;
                }
            }
            if (thirdPoint == null) {
                thirdPoint = hull[index];
            } else {
                double d = pointLineDistance(tailPoint, headPoint, hull[index]);
                if (d > thirdPointDist) {
                    thirdPointDist = d;
                    thirdPoint = hull[index];
                }
            }
            index++;
            visited++;
        }
        Point2d vMin = hull[min];
        Point2d tailProj = project(tailPoint, headPoint, vMax, true);
        Point2d headProj = project(tailPoint, headPoint, vMin, true);
//        System.out.println("vMax = " + hullIDs[max] + " vMin = " + hullIDs[min]);
        Rectangle r
                = new Rectangle(thirdPoint, tailProj, headProj, thirdPointDist);
        r.pointY = vMin;
        r.pointZ = vMax;
//        System.out.println(toString(tailPoint, headPoint, tailProj, headProj, vMin, vMax));
        return r;
    }

    private Rectangle getRectangle(
            Vector2d vector, int tailPointIndex, int headPointIndex) {

        Point2d headPoint = hull[headPointIndex];
        Point2d tailPoint = hull[tailPointIndex];

        // search backwards through the hull for an extremal point 
        int tailExPtIndex = tailPointIndex;
        Point2d tailExPt = hull[tailExPtIndex];
        boolean increasing = true;
        double proj = pointLineDistance(tailPoint, headPoint, tailExPt);
        while (increasing) {
            // get the next point, wrapping around if necessary
            int nextIndex;
            if (tailExPtIndex > 0) {
                nextIndex = tailExPtIndex - 1;
            } else {
                nextIndex = hull.length - 1;
            }
            Point2d nextPoint = hull[nextIndex];
            double nextProj = pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextProj > proj) {
                proj = nextProj;
                tailExPtIndex = nextIndex;
                tailExPt = nextPoint;
            } else {
                increasing = false;
            }
        }

        // convert the extremal point to a corner point by projecting it on V
        Vector2d negV = new Vector2d(vector);
        negV.negate();

        Point2d projTail = project(tailPoint, headPoint, tailExPt);

        // search forwards through the hull for an extremal point 
        int headExPtIndex = headPointIndex;
        Point2d headExPt = hull[headExPtIndex];
        increasing = true;
        proj = pointLineDistance(tailPoint, headPoint, headExPt);
        while (increasing) {
            // get the next point, wrapping around if necessary
            int nextIndex;
            if (headExPtIndex < hull.length - 1) {
                nextIndex = headExPtIndex + 1;
            } else {
                nextIndex = 0;
            }
            Point2d nextPoint = hull[nextIndex];
            double nextProj = pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextProj > proj) {
                proj = nextProj;
                headExPtIndex = nextIndex;
                headExPt = nextPoint;
            } else {
                increasing = false;
            }
        }

        // convert the extremal point to a corner point by projecting it on V
        Point2d projHead = project(tailPoint, headPoint, headExPt);

        // search forwards through the hull for the last extremal point
        int remainExPtIndex = headExPtIndex;
        Point2d remainExPoint = hull[remainExPtIndex];
        increasing = true;
        double dist = pointLineDistance(tailPoint, headPoint, remainExPoint);
        while (increasing) {
            int nextIndex;
            if (remainExPtIndex < hull.length - 1) {
                nextIndex = remainExPtIndex + 1;
            } else {
                nextIndex = 0;
            }
            Point2d nextPoint = hull[nextIndex];
            double nextDistance = pointLineDistance(tailPoint, headPoint, nextPoint);
            if (nextDistance > dist) {
                dist = nextDistance;
                remainExPtIndex = nextIndex;
                remainExPoint = nextPoint;
            } else {
                increasing = false;
            }
        }
//        System.out.println(toString(tailPoint, headPoint, remainExPoint, tailExPt, projTail, headExPt, projHead)); 

        return new Rectangle(remainExPoint, projTail, projHead,
                pointLineDistance(tailPoint, headPoint, remainExPoint));
    }

    /**
     *
     * @param points
     * @return
     */
    public String toString(Point2d... points) {
        String str = "[";
        for (Point2d point : points) {
            str += format("(%2.0f, %2.0f)", point.x, point.y);
        }
        return str + "]";
    }

    private Point2d project(Point2d p1, Point2d p2, Point2d p3) {
        return project(p1, p2, p3, false);
    }

    private Point2d project(Point2d p1, Point2d p2, Point2d p3, boolean outSeg) {
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        if (dx == 0 && dy == 0) {
            return new Point2d(p1);
        } else {
            double t = ((p3.x - p1.x) * dx + (p3.y - p1.y) * dy)
                    / (dx * dx + dy * dy);
            Point2d p;
            if (outSeg && (t > 0 && t < 1)) {
                if (t > 0.5) {
                    p = p2;
                } else {
                    p = p1;
                }
            } else {
                p = new Point2d(p1.x + (t * dx), p1.y + (t * dy));
            }
//            if (outSeg) System.out.println("projecting t = " + t + " " + toString(p1, p2, p3, p));
            return p;
        }
    }

    private double prj(Point2d p1, Point2d p2, Point2d p3) {
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        return ((p3.x - p1.x) * dx + (p3.y - p1.y) * dy) / (dx * dx + dy + dy);
    }

    private double pointLineDistance(Point2d p1, Point2d p2, Point2d p3) {
        Point2d p = project(p1, p2, p3);
        return p3.distance(p);
    }

    private Vector2d rotate(Vector2d vector, double angle) {
        Vector2d rotatedVector = new Vector2d();
        double cosTh = cos(angle);
        double sinTh = sin(angle);
        rotatedVector.x = cosTh * vector.x - sinTh * vector.y;
        rotatedVector.y = sinTh * vector.x + cosTh * vector.y;
        return rotatedVector;
    }

    private Vector2d edgeVector(Point2d fromPoint, Point2d toPoint) {
        Vector2d edge = new Vector2d(fromPoint);
        edge.sub(toPoint);
        return edge;
    }

    /**
     *
     * @return
     */
    public Rectangle2D getAxisAlignedMinimumBoundingRectangle() {
        double minX = MAX_VALUE;
        double minY = MAX_VALUE;
        double maxX = MIN_VALUE;
        double maxY = MIN_VALUE;
        for (Point2d point : hull) {
            if (point.x < minX) {
                minX = point.x;
            }
            if (point.y < minY) {
                minY = point.y;
            }
            if (point.x > maxX) {
                maxX = point.x;
            }
            if (point.y > maxY) {
                maxY = point.y;
            }
        }

        return new Rectangle2D.Double(minX, minY, maxX, maxY);
    }

    // use Graham's scan algorithm - implementation from wikipedia
    private void makeFromPoints(Point2d[] points) {
        if (points.length < 4) {
            hull = points;
            return;
        }
        int indexOfLowPoint = -1;
        Point2d lowPoint = null;
        for (int index = 0; index < points.length; index++) {
            Point2d current = points[index];
            if (indexOfLowPoint == -1 || current.y > lowPoint.y) {
                lowPoint = current;
                indexOfLowPoint = index;
            }
        }
        swap(points, 0, indexOfLowPoint);
        points = sortByPolarAngle(points);
//        sortByXCoord(points);
        int m = 2;
        for (int i = 3; i < points.length; i++) {
            while (ccw(points[m - 1], points[m], points[i]) <= 0) {
                if (m == 2) {
                    swap(points, m, i);
                    i++;
                } else {
                    m--;
                }
            }
            m++;
            swap(points, m, i);
        }
        hull = new Point2d[m];
        arraycopy(points, 0, hull, 0, m);
    }

    // allegedly, book 'Computational Geometry' has info on this
    // (Berkman & Schrieber, 2008)
//    private void sortByXCoord(Point2d[] points) {
//        Point2d ref = points[0];
//        Arrays.sort(points, new Comparator<Point2d>() {
//
//            @Override
//            public int compare(Point2d p0, Point2d p1) {
//                if (p0.x < p1.x) {
//                    return -1;
//                } else if (p0.x > p1.x) {
//                    return 1;
//                } else {
//                    return 0;
//                }
//            }
//            
//        });
//        points[0] = ref;
//    }
    private Point2d[] sortByPolarAngle(Point2d[] points) {
        Point2d ref = points[0];
        final Map<Point2d, Double> angles = new HashMap<>();
        angles.put(ref, 0.0);
        for (int pointIndex = 1; pointIndex < points.length; pointIndex++) {
            Point2d point = points[pointIndex];
            double angle = getAngle(ref, point);
            angles.put(point, angle);
        }
        sort(points, new Comparator<Point2d>() {

            @Override
            public int compare(Point2d p0, Point2d p1) {
                return angles.get(p0).compareTo(angles.get(p1));
            }

        });
        Point2d[] sortedPoints = new Point2d[points.length + 1];
        sortedPoints[0] = points[points.length - 1];
        arraycopy(points, 0, sortedPoints, 1, points.length);
        return sortedPoints;
    }

    private double getAngle(Point2d ref, Point2d point) {
//        double angle = Math.atan((point.y - ref.y) / (point.x - ref.x));
//        if (angle < 0) angle += Math.PI;
//        return angle;
        Vector2d rp = new Vector2d(ref);
        rp.sub(point);
        rp.normalize();
        return X_AXIS.angle(rp);
    }

    private void swap(Point2d[] points, int i, int j) {
        Point2d tmp = points[i];
        points[i] = points[j];
        points[j] = tmp;
    }

    private double ccw(Point2d p1, Point2d p2, Point2d p3) {
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    }

    @Override
    public Iterator<Point2d> iterator() {
        return asList(hull).iterator();
    }

    /**
     * A rectangle that may not be axis-aligned
     *
     */
    public class Rectangle {

        /**
         *
         */
        public Point2d pointX;

        /**
         *
         */
        public Point2d pointY;

        /**
         *
         */
        public Point2d pointZ;

        /**
         *
         */
        public Point2d cornerA;

        /**
         *
         */
        public Point2d cornerB;

        /**
         *
         */
        public Point2d cornerC;

        /**
         *
         */
        public Point2d cornerD;

        /**
         *
         * @param pointOnAB
         * @param cornerC
         * @param cornerD
         * @param distToCD
         */
        public Rectangle(Point2d pointOnAB, Point2d cornerC, Point2d cornerD, double distToCD) {
            pointX = new Point2d(pointOnAB);
            this.cornerC = new Point2d(cornerC);
            this.cornerD = new Point2d(cornerD);
            Vector2d cdVec = new Vector2d(cornerD);
            cdVec.sub(cornerC);
            Vector2d cdVecNormalized = new Vector2d(cdVec);
            if (cdVec.x != 0 && cdVec.y != 0) {
                cdVecNormalized.normalize();
            }
            Vector2d perp = new Vector2d(cdVecNormalized.y, -cdVecNormalized.x);
//            System.out.println(
//                    pointOnAB + " " +  cornerC + " " +  cornerD + " " +  distToCD
//                    + " " +  cdVec + " " +  perp);
            cornerA = new Point2d(cornerD);
            cornerA.scaleAdd(distToCD, perp, cornerA);
            cornerB = new Point2d(cornerC);
            cornerB.scaleAdd(distToCD, perp, cornerB);
        }

        /**
         *
         * @return
         */
        public double area() {
//           return getWidth() * getHeight();
            return new Point2d(cornerA).distance(new Point2d(cornerC))
                    * new Point2d(cornerC).distance(new Point2d(cornerD));
        }

        public String toString() {
            return format("[(%2.0f, %2.0f), (%2.0f, %2.0f), (%2.0f, %2.0f), (%2.0f, %2.0f)]",
                    cornerA.x, cornerA.y, cornerB.x, cornerB.y, cornerC.x, cornerC.y, cornerD.x, cornerD.y);
        }

        /**
         *
         * @return
         */
        public double getWidth() {
            Vector2d cd = new Vector2d(cornerC);
            cd.sub(cornerD);
            return cd.length();
        }

        /**
         *
         * @return
         */
        public Vector2d getMajorAxis() {
            Vector2d cd = new Vector2d(cornerC);
            cd.sub(cornerD);
            double cdLen = cd.length();
            Vector2d ad = new Vector2d(cornerA);
            ad.sub(cornerD);
            double adLen = ad.length();
            if (adLen > cdLen) {
                return ad;
            } else {
                return cd;
            }
        }

        /**
         *
         * @return
         */
        public double getHeight() {
            Vector2d ac = new Vector2d(cornerA);
            ac.sub(cornerC);
            return ac.length();
        }
    }

}
