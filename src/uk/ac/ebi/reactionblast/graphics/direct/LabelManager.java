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

package uk.ac.ebi.reactionblast.graphics.direct;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

public class LabelManager {


    private static final Vector2d POS_X = new Vector2d(1, 0);

    private static final Vector2d POS_Y = new Vector2d(0, 1);

    private static final Vector2d vN = new Vector2d(0, -1);

    private static final Vector2d vNE = new Vector2d(1, -1);

    private static final Vector2d vE = new Vector2d(1, 0);

    private static final Vector2d vSE = new Vector2d(1, 1);

    private static final Vector2d vS = new Vector2d(0, 1);

    private static final Vector2d vSW = new Vector2d(-1, 1);

    private static final Vector2d vW = new Vector2d(-1, 0);

    private static final Vector2d vNW = new Vector2d(-1, -1);
    private static final Logger LOG = Logger.getLogger(LabelManager.class.getName());

    private final Map<IAtom, BitSet> atomAnnotationPositions;

    public LabelManager() {
        atomAnnotationPositions = new HashMap<IAtom, BitSet>();
    }

    public String getAnnotationPositionsAsString(IAtom atom) {
        StringBuilder sb = new StringBuilder("|");
        BitSet positions = getAtomAnnotationPositions(atom);
        AnnotationPosition[] values = AnnotationPosition.values();
        for (int i = 0; i < values.length; i++) {
            if (positions.get(i)) {
                sb.append(values[i]);
                sb.append("|");
            }
        }
        return sb.toString();
    }

    public AnnotationPosition getNextSparePosition(IAtom atom) {
        return getNextSparePosition(getAtomAnnotationPositions(atom));
    }

    public AnnotationPosition getNextSparePosition(BitSet positions) {
        for (int i = 0; i < AnnotationPosition.values().length; i++) {
            if (positions.get(i)) {
            } else {
                return AnnotationPosition.values()[i];
            }
        }
        return null;
    }

    public Vector2d getVectorFromPosition(AnnotationPosition position) {
        switch (position) {
            case N:     return vN;
            case NE:    return vNE;
            case E:     return vE;
            case SE:    return vSE;
            case S:     return vS;
            case SW:    return vSW;
            case W:     return vW;
            case NW:    return vNW;
            default:    return vN;
        }
    }

    public Vector2d getLeftPerpendicularFromPosition(AnnotationPosition position) {
        switch (position) {
            case N:     return vW;
            case NE:    return vNW;
            case E:     return vN;
            case SE:    return vNE;
            case S:     return vE;
            case SW:    return vSE;
            case W:     return vS;
            case NW:    return vSW;
            default:    return vN;
        }
    }

    public Vector2d getRightPerpendicularFromPosition(AnnotationPosition position) {
        switch (position) {
            case N:     return vE;
            case NE:    return vSE;
            case E:     return vS;
            case SE:    return vSW;
            case S:     return vW;
            case SW:    return vNW;
            case W:     return vN;
            case NW:    return vNE;
            default:    return vS;
        }
    }

    public BitSet getAtomAnnotationPositions(IAtom atom) {
        if (atomAnnotationPositions.containsKey(atom)) {
            return atomAnnotationPositions.get(atom);
        } else {
            BitSet positions = new BitSet();
            atomAnnotationPositions.put(atom, positions);
            return positions;
        }
    }

    public void setUsedPosition(IAtom atom, AnnotationPosition position) {
        // TODO : restore to this line
//        getAtomAnnotationPositions(atom).set(position.ordinal());
        BitSet pos = getAtomAnnotationPositions(atom);
//        if (pos == null) {
//            System.out.println("pos null");
//        } else if (position == null) {
//            System.out.println("position null");
//        }
        pos.set(position.ordinal());
    }

    public AnnotationPosition alignmentToAnnotationPosition(int align) {
        switch (align) {
            case 1: return AnnotationPosition.E;
            case -1: return AnnotationPosition.W;
            case -2: return AnnotationPosition.N;
            case 2: return AnnotationPosition.S;
            default: return AnnotationPosition.E;
        }
    }

    public void addBondToAtomAnnotationPositions(IBond bond) {
        IAtom atom0 = bond.getAtom(0);
        IAtom atom1 = bond.getAtom(1);
        AnnotationPosition bondPosition;
        BitSet positions;

        positions = getAtomAnnotationPositions(atom0);
        bondPosition = calculateBondPosition(atom0, atom1);
        positions.set(bondPosition.ordinal());

        positions = getAtomAnnotationPositions(atom1);
        bondPosition = calculateBondPosition(atom1, atom0);
        positions.set(bondPosition.ordinal());

    }

    public AnnotationPosition calculateBondPosition(IAtom atomFrom, IAtom atomTo) {
        AnnotationPosition pos = calculateRelativePosition(atomFrom.getPoint2d(), atomTo.getPoint2d());
//        System.out.println("Relative position for " + atomFrom.getID() + " and " + atomTo.getID() + " is " + pos);
        return pos;
    }

    public AnnotationPosition calculateRelativePosition(Point2d fromPoint, Point2d toPoint) {
        Vector2d bondVector = new Vector2d(toPoint);
        bondVector.sub(fromPoint);
        bondVector.normalize();

        double xAng = Math.toDegrees(bondVector.angle(POS_X));
        double yAng = Math.toDegrees(bondVector.angle(POS_Y));
        if (xAng < 22.5 && (yAng > 67.5 && yAng < 115.5)) {
            return AnnotationPosition.E;
        } else if ((xAng > 22.5 && xAng < 67.5) && (yAng > 115.5 && yAng < 155.5)) {
            return AnnotationPosition.NE;
        } else if ((xAng > 67.5 && xAng < 115.5) && (yAng > 155.5)) {
            return AnnotationPosition.N;
        } else if ((xAng > 115.5 && xAng < 155.5) && (yAng > 115.5 && yAng < 155.5)) {
            return AnnotationPosition.NW;
        } else if (xAng > 155.5 && (yAng > 67.5 && yAng < 115.5)) {
            return AnnotationPosition.W;
        } else if ((xAng > 115.5 && xAng < 155.5) && (yAng > 22.5 && yAng < 67.5)) {
            return AnnotationPosition.SW;
        } else if ((xAng > 67.5 && xAng < 115.5) && yAng < 22.5) {
            return AnnotationPosition.S;
        } else if ((xAng > 22.5 && xAng < 67.5) && (yAng > 22.5 && yAng < 67.5)) {
            return AnnotationPosition.SE;
        }

        return AnnotationPosition.E;    // whatever
    }

    private void blockRingSegment(IAtom atom, List<AnnotationPosition> ringPositions) {
        BitSet positions = getAtomAnnotationPositions(atom);
        // erk
        if (ringPositions.size() != 2) {
            return;
        }
        AnnotationPosition a = ringPositions.get(0);
        AnnotationPosition b = ringPositions.get(1);
        if (positionsEqual(a, b, AnnotationPosition.N, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.NW.ordinal());
            positions.set(AnnotationPosition.W.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.N, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.NE.ordinal());
            positions.set(AnnotationPosition.E.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.S)) {
            positions.set(AnnotationPosition.W.ordinal());
            positions.set(AnnotationPosition.SW.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.S)) {
            positions.set(AnnotationPosition.E.ordinal());
            positions.set(AnnotationPosition.SE.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.W, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.SW.ordinal());
            positions.set(AnnotationPosition.S.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.E, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.SE.ordinal());
            positions.set(AnnotationPosition.S.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.E)) {
            positions.set(AnnotationPosition.N.ordinal());
            positions.set(AnnotationPosition.NE.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.W)) {
            positions.set(AnnotationPosition.NW.ordinal());
            positions.set(AnnotationPosition.N.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.NE)) {
            positions.set(AnnotationPosition.N.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.SW, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.S.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NW, AnnotationPosition.SW)) {
            positions.set(AnnotationPosition.W.ordinal());
        } else if (positionsEqual(a, b, AnnotationPosition.NE, AnnotationPosition.SE)) {
            positions.set(AnnotationPosition.E.ordinal());
        }

    }

    private boolean positionsEqual(AnnotationPosition a, AnnotationPosition b,
            AnnotationPosition c, AnnotationPosition d) {
        return (a == c && b == d) || (a == d && b == c);
    }


    public void addRingCenterToAtomAnnotationPosition(
            IAtom atom, List<IAtom> connectedAtomsInRing) {
        Point2d p1 = atom.getPoint2d();
        List<AnnotationPosition> ringPositions = new ArrayList<AnnotationPosition>();
        for (IAtom connectedAtom : connectedAtomsInRing) {
            Point2d p2 = connectedAtom.getPoint2d();
            ringPositions.add(calculateRelativePosition(p1, p2));
        }
        blockRingSegment(atom, ringPositions);
    }

    public boolean isUsed(IAtom atom, AnnotationPosition suggestedPosition) {
        int index = suggestedPosition.ordinal();
        return getAtomAnnotationPositions(atom).get(index);
    }

    public void reset() {
        atomAnnotationPositions.clear();
    }

    public enum AnnotationPosition {

        N, W, S, E, NW, NE, SW, SE
    }

}
