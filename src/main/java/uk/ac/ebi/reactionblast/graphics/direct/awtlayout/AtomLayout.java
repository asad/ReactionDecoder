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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import static java.lang.Math.min;
import static java.lang.String.valueOf;
import java.util.BitSet;
import java.util.HashMap;
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
import uk.ac.ebi.reactionblast.graphics.direct.LabelManager;
import uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.E;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.N;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.NE;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.NW;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.S;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.SE;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.SW;
import static uk.ac.ebi.reactionblast.graphics.direct.LabelManager.AnnotationPosition.W;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.E;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.NONE;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.R;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.S;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.Z;

/**
 * Layout an atom's symbol, and surrounding annotations, such as : charge,
 * implicit hydrogens, etc.
 *
 *
 * @author maclean
 *
 */
public class AtomLayout extends AbstractAWTLayout<IAtom> {

    private Font atomSymbolFont;
    private Font subscriptFont;
    private Font atomIDFont;
    private Font chiralSymbolFont;
    private final IAtomColorer atomColorer;
    private final LabelManager labelManager;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    /**
     *
     * @param parent
     * @param params
     * @param labelManager
     */
    public AtomLayout(AbstractAWTLayout parent, Params params, LabelManager labelManager) {
        super.parent = parent;
        setParams(params);
        this.labelManager = labelManager;
        atomColorer = new CDK2DAtomColors();
        chiralMap = new HashMap<>();
    }

    /**
     *
     * @param atom
     * @param g
     * @return
     */
    @Override
    public BoundsTree layout(IAtom atom, Graphics2D g) {
        currentObject = atom;
        String id = atom.getID();

        // mol reference needed for things like connected bonds
        IAtomContainer molecule = null;
        if (parent != null) {
            molecule = (IAtomContainer) parent.getCurrentObject();
        }

        boundsTree = new BoundsTree(atom.getID());
        if (molecule == null || shouldDraw(atom, molecule)) {
            boundsTree.add(id + ":symbol", layoutAtomSymbol(atom, g));
            if (isCharged(atom)) {
                Rectangle2D chargeBounds = layoutCharge(atom, g);
                boundsTree.add(id + ":charge", chargeBounds);
            }

            if (params.drawImplicitHydrogens) {
                Integer implicitHydrogenCount = atom.getImplicitHydrogenCount();
                if (implicitHydrogenCount != null
                        && implicitHydrogenCount > 0) {
                    int align = 1;
                    if (molecule != null) {
                        getBestAlignmentForLabel(molecule, atom);
                    }
                    AnnotationPosition suggestedPosition
                            = labelManager.alignmentToAnnotationPosition(align);

                    // special case for H2O
                    if (atom.getSymbol().equals("O")
                            && (molecule == null
                            || molecule.getConnectedBondsCount(atom) == 0)) {
                        suggestedPosition = W;
                    }

                    if (labelManager.isUsed(atom, suggestedPosition)) {
                        suggestedPosition = labelManager.getNextSparePosition(atom);
                    }
                    labelManager.setUsedPosition(atom, suggestedPosition);
                    Rectangle2D hBounds
                            = layoutImplicitHydrogens(atom, implicitHydrogenCount, suggestedPosition, g);
                    if (hBounds != null) {  // TODO - shouldn't be null!
                        boundsTree.add(id + ":hs", hBounds);
                    }
                }
            }
        } else if (params.drawRS && chiralMap.containsKey(atom)) {
            boundsTree.add(id + ":chiral", layoutChiralSymbol(atom, chiralMap.get(atom), g));
        } else {
            Point2d p = atom.getPoint2d();
            boundsTree.add(id + ":symbol", new Point2D.Double(p.x, p.y));
        }

        if (params.drawAtomID && molecule != null) {
            Rectangle2D idBounds = layoutAtomID(atom, molecule, g);
            if (idBounds != null) {
                boundsTree.add(id + ":id", idBounds);
            }
        }

        if (params.drawLonePairs && molecule != null) {
            int lonePairCount = 0;
            for (ILonePair lonePair : molecule.lonePairs()) {
                if (lonePair.contains(atom)) {
                    lonePairCount++;
                }
            }
            if (lonePairCount > 0) {
                Stroke stroke = g.getStroke();
                g.setStroke(new BasicStroke(0.05f));
                layoutElectronPairs(atom, molecule, lonePairCount, g);
                g.setStroke(stroke);
            }
        }
        return boundsTree;
    }

    /**
     *
     * @param atom
     * @param g
     * @return
     */
    public Rectangle2D layoutAtomSymbol(IAtom atom, Graphics2D g) {
        String text = atom.getSymbol();
        if (atom instanceof PseudoAtom) {
            text = ((IPseudoAtom) atom).getLabel();
        }
        g.setFont(atomSymbolFont);
        Point2d p = atom.getPoint2d();
        return layoutText(text, p, g);
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

    private Rectangle2D layoutChiralSymbol(IAtom atom, IStereoAndConformation chirality, Graphics2D g) {
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
        return layoutText(text, p, g);
    }

    /**
     *
     * @param atom
     * @param implicitHydrogenCount
     * @param pos
     * @param g
     * @return
     */
    public Rectangle2D layoutImplicitHydrogens(
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
        if (pos == AnnotationPosition.E) {
            double cx = p.x + (atomSymbolWidth / 2) + (hWidth / 2);
            double cy = p.y;

            totalHBounds = new Rectangle2D.Double(
                    cx - (hWidth / 2),
                    cy - (hHeight / 2),
                    hWidth, hHeight);
            if (implicitHydrogenCount > 1) {
                g.setFont(subscriptFont);
                String hCount = valueOf(implicitHydrogenCount);
                Rectangle2D subscriptBounds = getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();
                cx += (hWidth / 2) + (subscriptWidth / 2);
                cy += params.subscriptHeight;
                g.setFont(atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds.add(new Rectangle2D.Double(
                        cx - (subscriptWidth / 2),
                        cy - (subscriptHeight / 2),
                        subscriptWidth,
                        subscriptHeight));
            }
        } else if (pos == W) {

            float x;
            float y;
            if (implicitHydrogenCount > 1) {
                String hCount = valueOf(implicitHydrogenCount);
                g.setFont(subscriptFont);
                Rectangle2D subscriptBounds = getTextBounds(g, hCount);
                subscriptWidth = subscriptBounds.getWidth();

                x = (float) (pc.x - subscriptWidth);
                y = pc.y + params.subscriptHeight;
                g.setFont(atomSymbolFont);
                double subscriptHeight = subscriptBounds.getHeight();
                totalHBounds = new Rectangle2D.Double(
                        x - (subscriptWidth / 2),
                        y - (subscriptHeight / 2),
                        subscriptWidth, subscriptHeight);
            }

            x = (float) (pc.x - (atomSymbolWidth / 2) - subscriptWidth - (hWidth / 2));
            y = pc.y;
            Rectangle2D hDrawnBounds = new Rectangle2D.Double(
                    x - (hWidth / 2),
                    y - (hHeight / 2),
                    hWidth, hHeight);
            if (totalHBounds == null) {
                totalHBounds = hDrawnBounds;
            } else {
                totalHBounds.add(hDrawnBounds);
            }
        }
        return totalHBounds;
    }

    private Rectangle2D layoutText(String text, Point2d p, Graphics2D g) {
        Rectangle2D stringBounds = getTextBounds(g, text);
        double sW2 = stringBounds.getWidth() / 2;
        double sH2 = stringBounds.getHeight() / 2;
        double x = p.x - sW2;
        double y = p.y - sH2;
        return new Rectangle2D.Double(x, y, sW2 * 2, sH2 * 2);
    }

    /**
     *
     * @param atom
     * @param container
     * @param g
     * @return
     */
    public Rectangle2D layoutAtomID(IAtom atom, IAtomContainer container, Graphics2D g) {
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
            //                System.LOGGER.debug("position null for ID " + atomID);
        }

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
    public Rectangle2D layoutElectronPairs(
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

    private Rectangle2D layoutCharge(IAtom atom, Graphics2D g) {
        BitSet annotationPositions = labelManager.getAtomAnnotationPositions(atom);

        Integer formalCharge = atom.getFormalCharge();
        String chargeText = getChargeString(formalCharge);
        Rectangle2D atomBounds = getTextBounds(g, atom.getSymbol());
        Rectangle2D chargeBounds = getTextBounds(g, chargeText);

        Point2d atomPoint = atom.getPoint2d();
        Point2d chargePoint = new Point2d(atomPoint);
        double chargeDim = min(chargeBounds.getWidth(),
                chargeBounds.getHeight());

        // preferred position for charge is NE (superscript)
        chargePoint.x += (atomBounds.getWidth() / 2) + (chargeDim / 2);
        chargePoint.y -= (atomBounds.getHeight() / 2);
        annotationPositions.set(NE.ordinal());

        return new Rectangle2D.Double(
                chargePoint.x - (chargeBounds.getWidth() / 2),
                chargePoint.y - (chargeBounds.getHeight() / 2),
                chargeBounds.getWidth(),
                chargeBounds.getHeight());
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

    /**
     *
     * @param atom
     * @return
     */
    public Color colorForAtom(IAtom atom) {
        return atomColorer.getAtomColor(atom);
    }

    /**
     *
     * @param obj
     * @param rootLabel
     * @param graphics
     * @return
     */
    @Override
    public BoundsTree layout(IAtom obj, String rootLabel, Graphics2D graphics) {
        // TODO Auto-generated method stub
        // XXX not really used.
        return null;
    }

    /**
     *
     */
    public void reset() {
        labelManager.reset();
    }
}
