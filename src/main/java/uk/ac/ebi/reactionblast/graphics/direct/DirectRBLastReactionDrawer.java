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

import java.awt.Color;
import static java.awt.Color.BLACK;
import static java.awt.Color.GRAY;
import java.awt.Font;
import static java.awt.Font.PLAIN;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import static java.lang.Math.max;
import static java.lang.Math.min;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.graphics.direct.ColorRamp.getColors;
import uk.ac.ebi.reactionblast.graphics.direct.awtlayout.AbstractAWTReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.AbstractDirectReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;
import uk.ac.ebi.reactionblast.mapping.blocks.Block;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockPair;
import uk.ac.ebi.reactionblast.mapping.helper.RBlastReaction;

/**
 *
 * @author asad
 */
public class DirectRBLastReactionDrawer extends AbstractDirectDrawer {

    private final static ILoggingTool LOGGER
            = createLoggingTool(DirectRBLastReactionDrawer.class);

    private DirectReactionDrawer reactionDrawer;

    /**
     *
     * @param layout
     */
    public DirectRBLastReactionDrawer(AbstractDirectReactionLayout layout) {
        this(new Params(), layout);
    }

    /**
     *
     * @param params
     * @param layout
     */
    public DirectRBLastReactionDrawer(
            Params params, AbstractDirectReactionLayout layout) {
        setParams(params);
        reactionDrawer = new DirectReactionDrawer(params, layout);
    }

    /**
     *
     * @param params
     * @param layout
     * @param exactLayout
     */
    public DirectRBLastReactionDrawer(Params params,
            AbstractDirectReactionLayout layout,
            AbstractAWTReactionLayout exactLayout) {
        setParams(params);
        reactionDrawer = new DirectReactionDrawer(params, layout, exactLayout);
    }

    /**
     *
     * @return
     */
    public DirectReactionDrawer getReactionDrawer() {
        return reactionDrawer;
    }

    /**
     *
     * @param rBlastReaction
     * @param w
     * @param h
     * @return
     */
    public Image drawRBlastReaction(RBlastReaction rBlastReaction, int w, int h) {
        BufferedImage image = makeBlankImage(w, h);
        Graphics2D g = (Graphics2D) image.getGraphics();
        drawRBlastReaction(rBlastReaction, w, h, g);
        g.dispose();
        return image;
    }

    /**
     *
     * @param rBlastReaction
     * @param w
     * @param h
     * @param g
     * @return
     */
    public Rectangle2D drawRBlastReaction(RBlastReaction rBlastReaction, int w, int h, Graphics2D g) {
        AffineTransform originalTransform = g.getTransform();

        // layout the reaction
        IReaction reaction = rBlastReaction.getReaction();
        reactionDrawer.getLayout().shouldInvert = true;

        // mappings between the unique atom container labels and IDs
        Map<String, String> labelMap = reactionDrawer.makeLabelMap(reaction);

        BoundsTree boundsTree = reactionDrawer.getExactReactionBounds(reaction, g);
        List<String> labels = new ArrayList<>(labelMap.keySet());

        // bug to do with adding zero-dimensional bboxes to root
        Rectangle2D totalBoundingBox = boundsTree.getBounds(labels);
        boundsTree.setRoot(totalBoundingBox);   // ensure root is correct :(

        // calculate zoom, and center
        double zoom = reactionDrawer.calculateZoom(w, h, totalBoundingBox);
        BoundsTree centeredBoundsTree
                = reactionDrawer.centerOn(reaction, 0, 0, boundsTree);

        // bug to do with adding zero-dimensional bboxes to root
        totalBoundingBox = centeredBoundsTree.getBounds(labels);
        centeredBoundsTree.setRoot(totalBoundingBox); // ensure root is correct :(

        if (params.highlightSubgraphs) {
            highlightSubgraphs(rBlastReaction);
        }

        if (params.drawRS) {
            DirectMoleculeDrawer molDrawer = reactionDrawer.getMoleculeDrawer();
            molDrawer.addToChiralMap(rBlastReaction.getAtomStereoProductMap());
            molDrawer.addToChiralMap(rBlastReaction.getAtomStereoReactantMap());
        }

        double rxnWidth = totalBoundingBox.getWidth();
        double rxnHeight = totalBoundingBox.getHeight();

        double finalWidth = (zoom * rxnWidth) + (2 * params.borderX);
        double finalHeight = (zoom * rxnHeight) + (2 * params.borderY);

        // draw the actual reaction
        reactionDrawer.drawReaction(reaction, centeredBoundsTree, w, h, zoom, g);

        // paint the bond change marks
        drawBondChangeMarks(rBlastReaction, g);

        // paint the extra mapping stuff on top
        if (params.drawSubgraphBoxes) {
            drawSubgraphBoxes(rBlastReaction, g);
        }

        if (params.drawReactionID) {
            reactionDrawer.drawReactionID(reaction, (int) (w / zoom), g);
        }

        if (params.drawLabelPanel) {
//            System.out.println("drawing labels");
            g.setTransform(originalTransform);
            g.translate((w / 2), (h / 2));
            g.setFont(new Font(params.labelPanelFont, PLAIN, params.labelPanelFontSize));
            AffineTransform labelTransform = new AffineTransform();
            labelTransform.scale(zoom, zoom);

            double labelGap = params.labelGap;
            double labelHeight = reactionDrawer.getMaxLabelHeight(
                    centeredBoundsTree, labelMap, g);
            double labelShift = (totalBoundingBox.getHeight() / 2) + (labelHeight / 2) + labelGap;
            labelTransform.translate(0, labelShift);
            BoundsTree labelBoundsTree = centeredBoundsTree.transform(labelTransform);
            // TODO : label color selection
            g.setColor(BLACK);

            reactionDrawer.drawLabelPanel(labelMap, labelBoundsTree, g);
            finalHeight += labelHeight + labelGap;
        }

        if (params.shouldCrop) {
            double dX = w - finalWidth;
            double dY = h - finalHeight;
            int cropX = max(0, (int) dX / 2);
            int cropY = max(0, (int) dY / 2);
            int cropW = (int) min(finalWidth, w);
            int cropH = (int) min(finalHeight, w);
//            System.out.println("CROPPING totalBounds " 
//                    + BoundsPrinter.toString(totalBoundingBox));
//            System.out.println("zoom " + zoom + " dX " + dX + " dY " + dY 
//                            + " crop " + cropX + " " + cropY + " " 
//                                       + cropW + " " + cropH);
            if ((cropX + cropW > w) || (cropY + cropH > h)) {
                LOGGER.warn("Not cropping to ["
                        + cropX + ", " + cropY + "] "
                        + cropW + " x " + cropH + " as "
                        + (cropX + cropW) + " > " + w + " or "
                        + (cropY + cropH) + " > " + h);
                return new Rectangle2D.Double(0, 0, finalWidth, finalHeight);
            }
            return new Rectangle2D.Double(cropX, cropY, cropW, cropH);
        } else {
            return new Rectangle2D.Double(0, 0, finalWidth, finalHeight);
        }
    }

    /**
     *
     * @param rBlastReaction
     * @param g
     */
    public void drawSubgraphBoxes(RBlastReaction rBlastReaction, Graphics2D g) {
        List<Color> colors;
        if (params.colorSubgraphBoxes) {
            colors = getColors(rBlastReaction.getMappedSubgraphs().size());
        } else {
            colors = new ArrayList<>();
        }
        int blockIndex = 0;
        for (BlockPair subgraphMapping : rBlastReaction.getMappedSubgraphs()) {
            Block productBlock = subgraphMapping.getProductBlock();
            Block reactantBlock = subgraphMapping.getReactantBlock();
            Color color = getColorForBlock(colors, blockIndex);
            drawBlockBounds(productBlock, color, g);
            drawBlockBounds(reactantBlock, color, g);
            if (params.drawSubgraphMappingLines) {
                drawBlockMapping(productBlock, reactantBlock, g);
            }
            blockIndex++;
        }
    }

    /**
     *
     * @param rBlastReaction
     */
    public void highlightSubgraphs(RBlastReaction rBlastReaction) {
        List<Color> colors = getColors(rBlastReaction.getMappedSubgraphs().size());
        int blockIndex = 0;
        for (BlockPair subgraphMapping : rBlastReaction.getMappedSubgraphs()) {
            Block productBlock = subgraphMapping.getProductBlock();
            Block reactantBlock = subgraphMapping.getReactantBlock();
            Color color = getColorForBlock(colors, blockIndex);
            highlightSubgraph(reactantBlock, color);
            highlightSubgraph(productBlock, color);
            blockIndex++;
        }
    }

    private void highlightSubgraph(Block block, Color color) {
        reactionDrawer.highlightSubgraph(block.getAtoms(), color);
    }

    /**
     *
     * @param colors
     * @param index
     * @return
     */
    public Color getColorForBlock(List<Color> colors, int index) {
        if (colors.isEmpty() || index > colors.size() || !params.colorSubgraphBoxes) {
            return GRAY;
        } else {
            return colors.get(getWheelIndex(index, colors.size()));
        }
    }

    /**
     *
     * @param rBlastReaction
     * @param g
     */
    public void drawBondChangeMarks(RBlastReaction rBlastReaction, Graphics2D g) {

        if (params.drawBondFormedCleavedMarks) {
            List<IBond> bondCleavedReactant = rBlastReaction.getBondsCleavedInReactant();
            drawBondExistentialMarks(bondCleavedReactant, g);
        }

        if (params.drawBondOrderChangedMarks) {
            List<IBond> bondOrderChangedReactant = rBlastReaction.getBondsOrderChangedInReactant();
            drawBondChangeMarks(bondOrderChangedReactant, g);
        }

        if (params.drawBondStereoChanges) {
            List<IBond> bondStereoChangedReactant = rBlastReaction.getBondsStereoChangedInReactant();
            drawBondChangeMarks(bondStereoChangedReactant, g);
        }

        if (params.drawBondFormedCleavedMarks) {
            List<IBond> bondFormedProduct = rBlastReaction.getBondsFormedInProduct();
            drawBondExistentialMarks(bondFormedProduct, g);
        }

        if (params.drawBondOrderChangedMarks) {
            List<IBond> bondOrderChangedProduct = rBlastReaction.getBondsOrderChangedInProduct();
            drawBondChangeMarks(bondOrderChangedProduct, g);
        }

        if (params.drawBondStereoChanges) {
            List<IBond> bondStereoChangedProduct = rBlastReaction.getBondsStereoChangedInProduct();
            drawBondChangeMarks(bondStereoChangedProduct, g);
        }
    }

    private void drawBondExistentialMarks(List<IBond> bondsCleaved, Graphics2D g) {
        double markLength = params.bondMarkLength;
        for (IBond bond : bondsCleaved) {
            Point2d p1 = bond.getAtom(0).getPoint2d();
            Point2d p2 = bond.getAtom(1).getPoint2d();
            Point2d center = new Point2d(p1);
            center.interpolate(p2, 0.5);

            Vector2d bondVector = new Vector2d(p1);
            bondVector.sub(p2);
            bondVector.normalize();
            Vector2d negBondVector = new Vector2d(bondVector);
            negBondVector.negate();

            Point2d pc1 = new Point2d(center);
            pc1.scaleAdd(params.doubleMarkGap, bondVector, pc1);

            Point2d pc2 = new Point2d(center);
            pc2.scaleAdd(params.doubleMarkGap, negBondVector, pc2);

            Vector2d perp = new Vector2d(-bondVector.y, bondVector.x);
            Vector2d negPerp = new Vector2d(perp);
            negPerp.negate();

            Point2d pp11 = new Point2d(pc1);
            pp11.scaleAdd(markLength / 2, perp, pp11);
            Point2d pp12 = new Point2d(pc1);
            pp12.scaleAdd(markLength / 2, negPerp, pp12);

            drawLine(pp11, pp12, g);

            Point2d pp21 = new Point2d(pc2);
            pp21.scaleAdd(markLength / 2, perp, pp21);
            Point2d pp22 = new Point2d(pc2);
            pp22.scaleAdd(markLength / 2, negPerp, pp22);

            drawLine(pp21, pp22, g);
        }
    }

    private void drawBondChangeMarks(List<IBond> bondsChanged, Graphics2D g) {
        double markLength = params.bondMarkLength;
        for (IBond bond : bondsChanged) {
            Point2d p1 = bond.getAtom(0).getPoint2d();
            Point2d p2 = bond.getAtom(1).getPoint2d();
            Point2d center = new Point2d(p1);
            center.interpolate(p2, 0.5);

            Vector2d bondVector = new Vector2d(p1);
            bondVector.sub(p2);
            bondVector.normalize();

            Vector2d perp = new Vector2d(-bondVector.y, bondVector.x);
            Vector2d negPerp = new Vector2d(perp);
            negPerp.negate();

            Point2d pp1 = new Point2d(center);
            pp1.scaleAdd(markLength / 2, perp, pp1);
            Point2d pp2 = new Point2d(center);
            pp2.scaleAdd(markLength / 2, negPerp, pp2);

            drawLine(pp1, pp2, g);
        }
    }

    private void drawBlockMapping(
            Block productBlock, Block reactantBlock, Graphics2D g) {
        Point2d productCenter = productBlock.getCenterPoint();
        Point2d reactantCenter = reactantBlock.getCenterPoint();
        drawLine(productCenter, reactantCenter, g);
    }

    private int getWheelIndex(int index, int n) {
        int sum = 0;
        for (int i = 0; i < index; i++) {
            sum += (n / 2) + ((i % 2) * ((n + 1) % 2));
        }
        return sum % n;
    }

    private void drawBlockBounds(Block block, Color color, Graphics2D g) {
        Rectangle2D bounds = reactionDrawer.getDrawnBounds(block.getAtoms());
        if (bounds == null) {
            bounds = block.getBounds();
        }
        g.setColor(color);
        double centerX = bounds.getCenterX();
        double centerY = bounds.getCenterY();
        int w = ((int) bounds.getWidth()) + (params.subgraphBoxXBorder * 2);
        int h = ((int) bounds.getHeight()) + (params.subgraphBoxYBorder * 2);
        int x = ((int) centerX) - (w / 2);
        int y = ((int) centerY) - (h / 2);
        g.drawRect(x, y, w, h);
    }
}
