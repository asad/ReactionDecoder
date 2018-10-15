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
import java.awt.Color;
import static java.awt.Color.BLACK;
import static java.awt.Color.DARK_GRAY;
import static java.awt.Color.WHITE;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import static java.lang.Math.min;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import javax.vecmath.Point2f;
import javax.vecmath.Vector2d;
import org.openscience.cdk.PseudoAtom;
import static org.openscience.cdk.geometry.GeometryTools.getBestAlignmentForLabel;
import static org.openscience.cdk.geometry.GeometryTools.getBestAlignmentForLabelXY;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.renderer.color.CDK2DAtomColors;
import org.openscience.cdk.renderer.color.IAtomColorer;
import uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.E;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.N;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.NE;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.NW;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.S;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.SE;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.SW;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.W;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.E;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.NONE;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.R;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.S;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.Z;

/**
 *
 * @author asad
 */
public class DirectAtomDrawer extends AbstractDirectDrawer {

    private Font atomSymbolFont;
    private Font subscriptFont;
    private Font atomIDFont;
    private Font chiralSymbolFont;
    private final IAtomColorer atomColorer;
    private final Map<IAtom, Rectangle2D> drawnAtomBounds;
    private final LabelManager labelManager;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    /**
     *
     * @param params
     * @param labelManager
     */
    public DirectAtomDrawer(Params params, LabelManager labelManager) {
        setParams(params);
        this.labelManager = labelManager;
        atomColorer = new CDK2DAtomColors();
        drawnAtomBounds = new HashMap<>();
        chiralMap = new HashMap<>();
    }

    /**
     *
     * @param chiralMap
     */
    public void setChirals(Map<IAtom, IStereoAndConformation> chiralMap) {
        this.chiralMap = chiralMap;
    }

    /**
     *
     * @param atomSymbolFont
     */
    public void setAtomSymbolFont(Font atomSymbolFont) {
        this.atomSymbolFont = atomSymbolFont;
    }

    /**
     *
     * @param subscriptFont
     */
    public void setSubscriptFont(Font subscriptFont) {
        this.subscriptFont = subscriptFont;
    }

    /**
     *
     * @param atomIDFont
     */
    public void setAtomIDFont(Font atomIDFont) {
        this.atomIDFont = atomIDFont;
    }

    /**
     *
     * @param chiralSymbolFont
     */
    public void setChiralSymbolFont(Font chiralSymbolFont) {
        this.chiralSymbolFont = chiralSymbolFont;
    }

    /**
     *
     * @param atoms
     * @return
     */
    public Rectangle2D getDrawnBounds(List<IAtom> atoms) {
        Rectangle2D totalBounds = null;
        for (IAtom atom : atoms) {
            Rectangle2D bounds = drawnAtomBounds.get(atom);
            if (bounds == null) {
                continue;
            }
            if (totalBounds == null) {
                totalBounds = (Rectangle2D) bounds.clone();
            }
            totalBounds.add(bounds);
        }
        return totalBounds;

    }

    /**
     *
     * @param molecule
     * @param g
     */
    public void drawAtoms(IAtomContainer molecule, Graphics2D g) {
        Map<IAtom, Integer> lonePairMap = null;
        if (params.drawLonePairs) {
            lonePairMap = getLonePairCounts(molecule);
        }

        for (IAtom atom : molecule.atoms()) {
            int lonePairCount = 0;
            if (params.drawLonePairs) {
                Integer lonePairCountInteger = lonePairMap.get(atom);
                if (lonePairCountInteger == null) {
                    lonePairCount = 0;
                } else {
                    lonePairCount = lonePairCountInteger;
                }
            }
            drawnAtomBounds.put(atom, drawAtom(atom, molecule, lonePairCount, g));
        }
    }

    /**
     *
     * @param atom
     * @param molecule
     * @param lonePairCount
     * @param g
     * @return
     */
    public Rectangle2D drawAtom(
            IAtom atom, IAtomContainer molecule, int lonePairCount, Graphics2D g) {
        Rectangle2D symbolBounds;
        if (shouldDraw(atom, molecule)) {
            symbolBounds = drawAtomSymbol(atom, g);
            if (isCharged(atom)) {
                Rectangle2D chargeBounds = drawCharge(atom, g);
                symbolBounds.add(chargeBounds);
            }

            if (params.drawImplicitHydrogens) {
                Integer implicitHydrogenCount = atom.getImplicitHydrogenCount();
                if (implicitHydrogenCount != null
                        && implicitHydrogenCount > 0) {
                    int align
                            = getBestAlignmentForLabel(molecule, atom);
                    AnnotationPosition suggestedPosition
                            = labelManager.alignmentToAnnotationPosition(align);

                    // special case for H2O
                    if (atom.getSymbol().equals("O")
                            && molecule.getConnectedBondsCount(atom) == 0) {
                        suggestedPosition = W;
                    }

                    if (labelManager.isUsed(atom, suggestedPosition)) {
                        suggestedPosition = labelManager.getNextSparePosition(atom);
                    }
                    labelManager.setUsedPosition(atom, suggestedPosition);
                    Rectangle2D hBounds
                            = drawImplicitHydrogens(atom, implicitHydrogenCount, suggestedPosition, g);
                    if (hBounds != null) {  // TODO - shouldn't be null!
                        symbolBounds.add(hBounds);
                    }
                }
            }
        } else if (params.drawRS && chiralMap.containsKey(atom)) {
            symbolBounds = drawChiralSymbol(atom, chiralMap.get(atom), g);
        } else {
            Point2d p = atom.getPoint2d();
            symbolBounds = new Rectangle2D.Double(p.x, p.y, 0, 0);
        }

        if (params.drawAtomID) {
            Rectangle2D idBounds = drawAtomID(atom, molecule, g);
            if (idBounds != null) {
                symbolBounds.add(idBounds);
            }
        }

        if (params.drawLonePairs) {
            Stroke stroke = g.getStroke();
            g.setStroke(new BasicStroke(0.05f));
            drawElectronPairs(atom, molecule, lonePairCount, g);
            g.setStroke(stroke);
        }
        return symbolBounds;
    }

    private Rectangle2D drawChiralSymbol(IAtom atom, IStereoAndConformation chirality, Graphics2D g) {
        String text = "(-)";
        Point2d p = atom.getPoint2d();
        if (null != chirality) {
            switch (chirality) {
                case NONE:
                    return new Rectangle2D.Double(p.x, p.y, 0, 0);
                case R:
                    text = "(R)";
                    break;
                case S:
                    text = "(S)";
                    break;
                case E:
                    text = "(E)";
                    break;
                case Z:
                    text = "(Z)";
                    break;
                default:
                    text = "(-)";
                    break;
            }
        }
        g.setFont(chiralSymbolFont);
        Color color = DARK_GRAY;
        return drawText(text, p, color, g);
    }

    /**
     *
     * @param atom
     * @param implicitHydrogenCount
     * @param pos
     * @param g
     * @return
     */
    public Rectangle2D drawImplicitHydrogens(
            IAtom atom, int implicitHydrogenCount, AnnotationPosition pos, Graphics2D g) {

        String text = atom.getSymbol();
        Point2d p = atom.getPoint2d();

        g.setFont(atomSymbolFont);
        Point2f pc = getTextPoint(g, text, p.x, p.y);

        // width on screen of the text
        Rectangle2D hBounds = getTextBounds(g, "H");
        double atomSymbolWidth = getTextBounds(g, text).getWidth();
        double hWidth = hBounds.getWidth();
        double hHeight = hBounds.getHeight();
        double subscriptWidth = 0;
        Rectangle2D totalHBounds = null;
        g.setColor(Color.BLACK);
        if (pos == AnnotationPosition.E) {
            double cx = p.x + (atomSymbolWidth / 2) + (hWidth / 2);
            double cy = p.y;
            Point2f hP = getTextPoint(g, "H", cx, cy);

            String hString = "H";
            g.drawString(hString, hP.x, hP.y);
            totalHBounds = new Rectangle2D.Double(
                    cx - (hWidth / 2),
                    cy - (hHeight / 2),
                    hWidth, hHeight);
            if (implicitHydrogenCount > 1) {
                g.setFont(subscriptFont);
                String hCount = String.valueOf(implicitHydrogenCount);
                Rectangle2D subscriptBounds = getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();
                cx += (hWidth / 2) + (subscriptWidth / 2);
                cy += params.subscriptHeight;
                Point2f sP = getTextPoint(g, hCount, cx, cy);
                double subscriptHeight = subscriptBounds.getHeight();
                Rectangle2D finalHBounds
                        = new Rectangle2D.Double(
                                cx - (subscriptWidth / 2),
                                cy - (subscriptHeight / 2),
                                subscriptWidth,
                                subscriptHeight);
                g.setColor(WHITE);
                g.fill(finalHBounds);
                g.setColor(BLACK);
                g.drawString(hCount, sP.x, sP.y);
                g.setFont(atomSymbolFont);
                totalHBounds.add(finalHBounds);
            }
        } else if (pos == W) {

            float x;
            float y;
            if (implicitHydrogenCount > 1) {
                String hCount = String.valueOf(implicitHydrogenCount);
                g.setFont(subscriptFont);
                Rectangle2D subscriptBounds = getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();

                x = (float) (pc.x - subscriptWidth);
                y = pc.y + params.subscriptHeight;
                g.drawString(hCount, x, y);
                g.setFont(atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds = new Rectangle2D.Double(
                        x - (subscriptWidth / 2),
                        y - (subscriptHeight / 2),
                        subscriptWidth, subscriptHeight);
            }

            x = (float) (pc.x - (atomSymbolWidth / 2) - subscriptWidth - (hWidth / 2));
            y = pc.y;
            String hString = "H";
            Rectangle2D hDrawnBounds = new Rectangle2D.Double(
                    p.x - (atomSymbolWidth / 2) - subscriptWidth - hWidth,
                    p.y - (hBounds.getHeight() / 2),
                    hWidth, hHeight);
            g.setColor(WHITE);
            g.fill(hDrawnBounds);
            g.setColor(BLACK);
            g.drawString(hString, x, y);
            if (totalHBounds == null) {
                totalHBounds = hDrawnBounds;
            } else {
                totalHBounds.add(hDrawnBounds);
            }
        }
        return totalHBounds;
    }

    /**
     *
     * @param atom
     * @param g
     * @return
     */
    public Rectangle2D drawAtomSymbol(IAtom atom, Graphics2D g) {
        String text = atom.getSymbol();
        if (atom instanceof PseudoAtom) {
            text = ((IPseudoAtom) atom).getLabel();
        }
        g.setFont(atomSymbolFont);
        Point2d p = atom.getPoint2d();
        return drawText(text, p, colorForAtom(atom), g);
    }

    private Rectangle2D drawText(String text, Point2d p, Color color, Graphics2D g) {
        Point2f pc = getTextPoint(g, text, p.x, p.y);
        Rectangle2D stringBounds = getTextBounds(g, text);
        double sW2 = stringBounds.getWidth() / 2;
        double sH2 = stringBounds.getHeight() / 2;
        double x = p.x - sW2;
        double y = p.y - sH2;
        g.setColor(WHITE);
        Rectangle2D bounds = new Rectangle2D.Double(x, y, sW2 * 2, sH2 * 2);
        g.fill(bounds);
        g.setColor(color);
        g.drawString(text, pc.x, pc.y);
        return bounds;
    }

    /**
     *
     * @param atom
     * @param container
     * @param g
     * @return
     */
    public Rectangle2D drawAtomID(IAtom atom, IAtomContainer container, Graphics2D g) {
        String atomID = atom.getID();

        if (atomID == null) {
            return null;
        }
        g.setFont(atomSymbolFont);
        Rectangle2D atomSymbolBounds;
        Point2d p = atom.getPoint2d();
        if (shouldDraw(atom, container)) {
            atomSymbolBounds = getTextBounds(g, atom.getSymbol());
        } else {
            atomSymbolBounds = new Rectangle2D.Double(p.x, p.y, 1, 1);
        }
        g.setFont(atomIDFont);
        Rectangle2D bounds = getTextBounds(g, atomID);
        Point2d pID = new Point2d(p);
        AnnotationPosition suggestedPosition
                = labelManager.alignmentToAnnotationPosition(getBestAlignmentForLabelXY(container, atom));
        AnnotationPosition pos;
        if (labelManager.isUsed(atom, suggestedPosition)) {
            pos = labelManager.getNextSparePosition(atom);
        } else {
            pos = suggestedPosition;
        }

        //        System.out.println("Alignment for atom " + atomID + " " + pos
        //                + " given annotations at "
        //                + labelManager.getAnnotationPositionsAsString(atom));
        double aW2 = atomSymbolBounds.getWidth() / 2;
        double bW2 = bounds.getWidth() / 2;
        double aH2 = atomSymbolBounds.getHeight() / 2;
        double bH2 = bounds.getHeight() / 2;

        if (null != pos) {
            switch (pos) {
                case N:
                    pID.y -= aH2 + bH2;
                    break;
                case NE:
                    pID.x += aW2 + bW2;
                    pID.y -= aH2 + bH2;
                    break;
                case E:
                    pID.x += aW2 + bW2;
                    break;
                case SE:
                    pID.x += aW2 + bW2;
                    pID.y += aH2 + bH2;
                    break;
                case S:
                    pID.y += aH2 + bH2;
                    break;
                case SW:
                    pID.x -= aW2 + bW2;
                    pID.y += aH2 + bH2;
                    break;
                case W:
                    pID.x -= aW2 + bW2;
                    break;
                case NW:
                    pID.x -= aW2 + bW2;
                    pID.y -= aH2 + bH2;
                    break;
                default:
                    pID.x += aW2 + bW2;
                    break;
            }
        }

        if (pos != null) {
            labelManager.setUsedPosition(atom, pos);
        } else {
//            System.LOGGER.debug("position null for ID " + atomID);
        }

        Point2f tp = getTextPoint(g, atomID, pID.x, pID.y);
        g.setColor(BLACK);
        g.drawString(atomID, tp.x, tp.y);
        g.setFont(atomSymbolFont);

        return new Rectangle2D.Double(
                pID.x - (bounds.getWidth() / 2),
                pID.y - (bounds.getHeight() / 2),
                bounds.getWidth(),
                bounds.getHeight());
    }

    /**
     *
     * @param atom
     * @param container
     * @param lonePairCount
     * @param g
     * @return
     */
    public Rectangle2D drawElectronPairs(
            IAtom atom, IAtomContainer container,
            int lonePairCount, Graphics2D g) {
        if (lonePairCount == 0) {
            return null;
        }

        Point2d atomPoint = atom.getPoint2d();
        Rectangle2D atomSymbolBounds = getTextBounds(g, atom.getSymbol());
        BitSet positions = labelManager.getAtomAnnotationPositions(atom);

        double r = params.electronRadius;
        double d = r * 2;
        for (int i = 0; i < lonePairCount; i++) {
            AnnotationPosition position = labelManager.getNextSparePosition(positions);
            Vector2d v = labelManager.getVectorFromPosition(position);
            Vector2d leftPerp = labelManager.getLeftPerpendicularFromPosition(position);
            Vector2d rightPerp = labelManager.getRightPerpendicularFromPosition(position);

            double dx = ((atomSymbolBounds.getWidth() / 2) + d) * v.x;
            double dy = ((atomSymbolBounds.getHeight() / 2) + d) * v.y;

            Point2d lp = new Point2d(atomPoint.x + dx, atomPoint.y + dy);
            Point2d llp = new Point2d(lp);
            llp.scaleAdd(params.lonePairSeparation / 2, leftPerp, llp);
            Point2d rlp = new Point2d(lp);
            rlp.scaleAdd(params.lonePairSeparation / 2, rightPerp, rlp);

            g.fill(new Ellipse2D.Double(llp.x - r, llp.y - r, d, d));
            g.fill(new Ellipse2D.Double(rlp.x - r, rlp.y - r, d, d));

            positions.set(position.ordinal());
        }
        return null;
    }

    private boolean shouldDraw(IAtom atom, IAtomContainer atomContainer) {
        String symbol = atom.getSymbol();
        if (symbol.equals("C")) {
            if (params.drawCarbons) {
                return true;
            } else if (params.drawTerminalCarbons
                    && isTerminal(atom, atomContainer)) {
                return true;
            } else {
                return getAttachedMultipleBondCount(atom, atomContainer) > 1;
            }
        } else if (symbol.equals("H")) {
            return params.drawExplicitHydrogens;
        }
        return true;
    }

    private int getAttachedMultipleBondCount(
            IAtom atom, IAtomContainer atomContainer) {
        int count = 0;
        for (IBond bond : atomContainer.getConnectedBondsList(atom)) {
            if (bond.getOrder() != SINGLE) {
                count++;
            }
        }
        return count;
    }

    /**
     *
     * @param atom
     * @return
     */
    public boolean isCharged(IAtom atom) {
        Integer formalCharge = atom.getFormalCharge();
        return formalCharge != null && formalCharge != 0;
    }

    private boolean isTerminal(IAtom atom, IAtomContainer atomContainer) {
        int numberOfHeavyAtomsConnected = 0;
        for (IAtom connected : atomContainer.getConnectedAtomsList(atom)) {
            if (!connected.getSymbol().equals("H")) {
                numberOfHeavyAtomsConnected++;
            }
        }
        return numberOfHeavyAtomsConnected < 2;
    }

    private Rectangle2D drawCharge(IAtom atom, Graphics2D g) {
        BitSet annotationPositions = labelManager.getAtomAnnotationPositions(atom);

        Integer formalCharge = atom.getFormalCharge();
        String chargeText = getChargeString(formalCharge);
        Rectangle2D atomBounds = getTextBounds(g, atom.getSymbol());
        Rectangle2D chargeBounds = getTextBounds(g, chargeText);
        g.setColor(BLACK);

        Point2d atomPoint = atom.getPoint2d();
        Point2d chargePoint = new Point2d(atomPoint);
        double chargeDim = min(chargeBounds.getWidth(),
                chargeBounds.getHeight());

        // preferred position for charge is NE (superscript)
        chargePoint.x += (atomBounds.getWidth() / 2) + (chargeDim / 2);
        chargePoint.y -= (atomBounds.getHeight() / 2);
        annotationPositions.set(NE.ordinal());

        Point2f sp = getTextPoint(g, chargeText, chargePoint.x, chargePoint.y);
        Rectangle2D chargeBox = new Rectangle2D.Double(
                chargePoint.x - (chargeBounds.getWidth() / 2),
                chargePoint.y - (chargeBounds.getHeight() / 2),
                chargeBounds.getWidth(),
                chargeBounds.getHeight());
        g.setColor(WHITE);
        g.fill(chargeBox);
        g.setColor(BLACK);
        g.drawString(chargeText, sp.x, sp.y);
        return chargeBox;
    }

    private String getChargeString(Integer formalCharge) {
        if (formalCharge == 1) {
            return "+";
        } else if (formalCharge == -1) {
            return "-";
        } else if (formalCharge > 1) {
            return formalCharge + "+";
        } else if (formalCharge < -1) {
            return formalCharge + "-";
        } else {
            return "";
        }
    }

    private Map<IAtom, Integer> getLonePairCounts(IAtomContainer atomContainer) {
        Map<IAtom, Integer> lonePairMap = new HashMap<>();
        for (ILonePair lonePair : atomContainer.lonePairs()) {
            IAtom atom = lonePair.getAtom();
            int lonePairCount;
            if (lonePairMap.containsKey(atom)) {
                lonePairCount = lonePairMap.get(atom);
            } else {
                lonePairCount = 0;
            }
            lonePairMap.put(atom, lonePairCount + 1);
        }
        return lonePairMap;
    }

    /**
     *
     * @param atom
     * @return
     */
    public Color colorForAtom(IAtom atom) {
        return atomColorer.getAtomColor(atom);
    }
}
