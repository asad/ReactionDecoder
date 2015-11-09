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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import javax.vecmath.Point2f;
import javax.vecmath.Vector2d;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.graphics.direct.Params.XAlign;
import uk.ac.ebi.reactionblast.graphics.direct.Params.YAlign;
import uk.ac.ebi.reactionblast.graphics.direct.awtlayout.AbstractAWTReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.AbstractDirectReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsPrinter;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;
import uk.ac.ebi.reactionblast.graphics.direct.layout.LeftToRightReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.TopToBottomReactionLayout;

/**
 * Directly draw, modifying the points in the atom containers.
 * 
 * @author maclean
 *
 */
public class DirectReactionDrawer extends AbstractDirectDrawer {

    private AbstractDirectReactionLayout reactionLayout;
    private AbstractAWTReactionLayout exactReactionLayout;
    private DirectMoleculeDrawer moleculeDrawer;
    private DirectArrowDrawer arrowDrawer;

    public DirectReactionDrawer(AbstractDirectReactionLayout layout) {
        this(new Params(), layout);
    }

    public DirectReactionDrawer(Params params) {
        setParams(params);
        
        // XXX FIXME
        reactionLayout = null;
        moleculeDrawer = new DirectMoleculeDrawer(params);
        arrowDrawer = new DirectArrowDrawer(params);
    }
    
    public DirectReactionDrawer(Params params, AbstractDirectReactionLayout layout) {
        this(params, layout, null);
    }
    
    public DirectReactionDrawer(Params params, 
                                AbstractDirectReactionLayout layout,
                                AbstractAWTReactionLayout exactReactionLayout) {
        setParams(params);
        
        reactionLayout = layout;
        layout.setParams(params);
        this.exactReactionLayout = exactReactionLayout;
        if (exactReactionLayout != null) {
            exactReactionLayout.setParams(params);
        }
        
        moleculeDrawer = new DirectMoleculeDrawer(params);
        arrowDrawer = new DirectArrowDrawer(params);
    }

    public AbstractDirectReactionLayout getLayout() {
        return reactionLayout;
    }

    private void setupLayout() {
        if (reactionLayout == null) {
            if (params.layoutLeftToRight) {
                reactionLayout = new LeftToRightReactionLayout();
            } else {
                reactionLayout = new TopToBottomReactionLayout();
            }
        }
    }

    public Image drawReaction(IReaction reaction) {
        return drawReaction(reaction, true);
    }

    /**
     * Draw a zoomed reaction.
     * 
     * @param reaction
     * @param w
     * @param h
     * @return
     */
    public Image drawReaction(IReaction reaction, int w, int h) {
        return drawReaction(reaction, w, h, true);
    }
    
    public Map<String, String> makeLabelMap(IReaction reaction) {
        Map<String, String> labelMap = new HashMap<String, String>();
        
        // annoying, but necessary for bounds labels
        reaction.getReactants().setID("r");
        reaction.getProducts().setID("p");
        
        String rxnID = reaction.getID();
        
        int counter = 0;
        for (IAtomContainer atomContainer : reaction.getReactants().atomContainers()) {
            String acID = atomContainer.getID(); 
            String boundsLabel = rxnID + "_" + "r" + "_" + acID + ":" + counter;
            labelMap.put(boundsLabel, acID);
            counter++;
        }
        counter = 0;
        for (IAtomContainer atomContainer : reaction.getProducts().atomContainers()) {
            String acID = atomContainer.getID(); 
            String boundsLabel = rxnID + "_" + "p" + "_" + acID + ":" + counter;
            labelMap.put(boundsLabel, acID);
            counter++;
        }
        return labelMap;
    }

    /**
     * Draw a zoomed reaction.
     * 
     * @param reaction
     * @param w
     * @param h
     * @param invert
     * @return
     */
    public Image drawReaction(IReaction reaction, int w, int h, boolean invert) {
        
        setupLayout();
        reactionLayout.shouldInvert = invert;
        
        // mappings between the unique atom container labels and IDs
        Map<String, String> labelMap = makeLabelMap(reaction);

        BoundsTree boundsTree =
                reactionLayout.layout(reaction, reactionLayout.getAxis());
        Vector2d reactionAxis = reactionLayout.getAxis();
        
        BufferedImage image = makeBlankImage(w, h);
        
        Graphics2D g = (Graphics2D) image.getGraphics();
        if (params.useAntialias) {
            g.setRenderingHint(
                    RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }
        
        // store the original transform, for use in drawing label panels
        AffineTransform originalTransform = g.getTransform();
        
        if (exactReactionLayout != null) {
            boundsTree = exactReactionLayout.layout(reaction, g);
        }

//        Rectangle2D totalBoundingBox = boundsTree.getRoot(); // XXX BUG
        Rectangle2D totalBoundingBox = boundsTree.getBounds(new ArrayList<String>(labelMap.keySet()));
        boundsTree.setRoot(totalBoundingBox);   // ensure root is correct :(
        
        double zoom = calculateZoom(w, h, totalBoundingBox);
        BoundsTree centeredBoundsTree = centerOn(reaction, 0, 0, boundsTree);
        
//      totalBoundingBox = centeredBoundsTree.getRoot(); // XXX BUG
        totalBoundingBox = centeredBoundsTree.getBounds(new ArrayList<String>(labelMap.keySet()));
        centeredBoundsTree.setRoot(totalBoundingBox); // ensure root is correct :(
        
        double rxnWidth  = totalBoundingBox.getWidth();
        double rxnHeight = totalBoundingBox.getHeight();

        double finalWidth  = (zoom * rxnWidth)  + (2 * params.borderX);
        double finalHeight = (zoom * rxnHeight) + (2 * params.borderY);
        
        // move the origin to the middle of the canvas, and scale 
        g.translate(w / 2, h / 2);
        g.scale(zoom, zoom);
        
        // if there is a label panel, further transform to make space for it
        if (params.drawLabelPanel) {
            // for now, assumes label panel is always at the bottom
            finalHeight += params.labelPanelHeight;
        }
        
        // actually draw the reaction
        drawReaction(reaction, g, centeredBoundsTree, reactionAxis);

        if (params.drawLabelPanel) {
            // the original, unchanged transform
            g.setTransform(originalTransform);
            g.translate((w/2), (h/2));
            g.setFont(new Font(params.labelPanelFont, Font.PLAIN, params.labelPanelFontSize));
            AffineTransform labelTransform = new AffineTransform();
            labelTransform.scale(zoom, zoom);
            
            double labelGap = params.labelGap;
            double labelHeight = getMaxLabelHeight(centeredBoundsTree, labelMap, g);
            double labelShift = (totalBoundingBox.getHeight() / 2) + (labelHeight / 2) + labelGap;
            labelTransform.translate(0, labelShift);
            BoundsTree labelBoundsTree = centeredBoundsTree.transform(labelTransform);
            // TODO : label color selection
            g.setColor(Color.BLACK);
            
            drawLabelPanel(labelMap, labelBoundsTree, g);
            finalHeight += labelHeight + labelGap;
        }
        
        // alter the image object
        if (params.shouldCrop) {
            double dX = w - finalWidth;
            double dY = h - finalHeight;
            int cropX = Math.max(0, (int) dX / 2);
            int cropY = Math.max(0, (int) dY / 2);
            int cropW = (int) Math.min(finalWidth, w);
            int cropH = (int) Math.min(finalHeight, w);
//            System.out.println("CROPPING totalBounds " 
//                    + BoundsPrinter.toString(totalBoundingBox));
//            System.out.println("zoom " + zoom + " dX " + dX + " dY " + dY 
//                            + " crop " + cropX + " " + cropY + " " 
//                                       + cropW + " " + cropH);
            if ((cropX + cropW > w) || (cropY + cropH > h)) {
                System.out.println("Not cropping " + (cropX + cropW) + " " + w + " "
                                                   + (cropY + cropH) + " " + h);
                return image;
            }
            return image.getSubimage(cropX, cropY, cropW, cropH);
        } else {
            return image;
        }
    }

    /**
     * Draw a natural-scale reaction.
     * 
     * @param reaction
     * @param invert
     * @return
     */
    public Image drawReaction(IReaction reaction, boolean invert) {
        setupLayout();
        int borderX = params.borderX;
        int borderY = params.borderY;

        reactionLayout.shouldInvert = invert;
        
        if (exactReactionLayout != null) {
            reaction.getReactants().setID("r");
            reaction.getProducts().setID("p");
        }

        BoundsTree boundsTree =
                reactionLayout.layout(reaction, reactionLayout.getAxis());
        Vector2d reactionAxis = reactionLayout.getAxis();

        if (exactReactionLayout != null) {
            Image dummyImage = makeBlankImage(1, 1);
            Graphics2D g = (Graphics2D) dummyImage.getGraphics();
            boundsTree = exactReactionLayout.layout(reaction, g);
            Rectangle2D bb = boundsTree.getRoot();
            double dx = (boundsTree.getWidth()/2) -bb.getCenterX();
            double dy = (boundsTree.getHeight()/2) -bb.getCenterY();
            System.err.println(BoundsPrinter.toString(bb) + " " + dx + " " + dy);
            
            // AARGH!
            boundsTree = shift(reaction, boundsTree, dx, 0);
        }
        
        int w = (int)boundsTree.getWidth() + (2 * borderX);
        int h = (int)boundsTree.getHeight() + (2 * borderY);
        Image image = makeBlankImage(w, h);
        Graphics2D g = (Graphics2D) image.getGraphics();
        drawReaction(reaction, g, boundsTree, reactionAxis);
        return image;
    }

    /**
     * Draw a zoomed reaction.
     * 
     * @param reaction
     * @param boundsTree
     * @param w
     * @param h
     * @param zoom
     * @param g
     */
    public void drawReaction(
            IReaction reaction, BoundsTree boundsTree, int w, int h, double zoom, Graphics2D g) {
        if (params.useAntialias) {
            g.setRenderingHint(
                    RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }
        g.translate(w / 2, h / 2);
        g.scale(zoom, zoom);
        Vector2d reactionAxis;
        if (reactionLayout != null) {
            if (exactReactionLayout == null) { 
                reactionAxis = reactionLayout.getAxis();
            } else {
                reaction.getReactants().setID("r");
                reaction.getProducts().setID("p");
                reactionAxis = exactReactionLayout.getAxis();
            }
        } else {
            reactionAxis = new Vector2d(1, 0);
        }
        drawReaction(reaction, g, boundsTree, reactionAxis);
    }
    
    public void drawLabelPanel(Map<String, String> labelMap, 
                               BoundsTree labelBoundsTree, Graphics2D g) {
        MoleculeLabelDrawer molLabelDrawer = 
            new MoleculeLabelDrawer(Axis.X, params);
        molLabelDrawer.draw(labelMap, labelBoundsTree, g);
    }
    
    public double getMaxLabelHeight(BoundsTree tree, Map<String, String> labels, Graphics2D g) {
        double maxHeight = 0;
        Font font = new Font(params.labelPanelFont, Font.PLAIN, params.labelPanelFontSize);
        FontRenderContext frc = g.getFontRenderContext();
        for (String boundsLabel : labels.keySet()) {
            String label = labels.get(boundsLabel); 
            Rectangle2D bounds = tree.get(boundsLabel);
            float boundsWidth = (float)bounds.getWidth();
            TextLayout textLayout = new TextLayout(label, font, frc);
            if (boundsWidth <= 0) continue; // XXX
            TextLayout justifiedLayout = textLayout.getJustifiedLayout(boundsWidth);
            double height = justifiedLayout.getBounds().getHeight();
            if (height > maxHeight) {
                maxHeight = height;
            }
        }
        return maxHeight;
    }
        

    private void drawReaction(IReaction reaction, 
            Graphics2D g, BoundsTree boundsTree, Vector2d reactionAxis) {
        
        IAtomContainerSet reactants = reaction.getReactants();
        //      int reactantAxisPos = reactionLayout.getReactantAxisPos();     // TODO
        double reactantAxisPos = molSetBounds(reactants).getCenterY(); // TODO
        
        String reactionID = reaction.getID();
        drawMoleculeSet(reactants, reactionID, reactantAxisPos, boundsTree, g);
        
        // FIXME
        Point2d arrowCenter;
        int index = reaction.getReactantCount() - 1;
        IAtomContainer lastReactant = 
            reaction.getReactants().getAtomContainer(index);
        String reactantsID = reaction.getReactants().getID(); 
        String boundsID = reactionID + "_" + reactantsID + "_" + lastReactant.getID() + ":" + index;
//        System.out.println("looking up bounds " + boundsID);
        Rectangle2D lastRBound = boundsTree.get(boundsID);
        double xPos;
        double yPos;
        if (reactionLayout.getArrowAxis() == Axis.X) {
            xPos = lastRBound.getMaxX() + params.arrowGap + (params.arrowLength / 2);
            yPos = boundsTree.getRoot().getCenterY();
        } else {
            xPos = boundsTree.getRoot().getCenterX();
            yPos = lastRBound.getMaxY() + params.arrowGap + (params.arrowLength / 2);
        }
        arrowCenter = new Point2d(xPos, yPos); 
        
//        System.out.println("arrow center @ " + arrowCenter);
        if (params.arrowType == Params.ArrowType.FORWARD) {
            arrowDrawer.drawArrow(g, arrowCenter, reactionAxis);
        } else if (params.arrowType == Params.ArrowType.BACKWARD) {
            Vector2d backAxis = new Vector2d(reactionAxis);
            backAxis.negate();
            arrowDrawer.drawArrow(g, arrowCenter, backAxis);
        } else {
            arrowDrawer.drawArrow(g, arrowCenter, reactionAxis);
        }

        IAtomContainerSet products = reaction.getProducts();
        //      int productAxisPos = reactionLayout.getProductAxisPos();     // TODO
        double productAxisPos = molSetBounds(products).getCenterY(); // TODO
        drawMoleculeSet(products, reactionID, productAxisPos, boundsTree, g);

        if (params.drawMappings) {
            drawMappings(reaction, g);
        }
    }
    
    /**
     * DEBUG method to show the bounds tree. 
     *  
     * @param tree
     * @param g
     */
    private void drawBoundsTree(BoundsTree tree, List<String> labels, Color color, Graphics2D g) {
        java.util.Random random = new java.util.Random();
        for (String label : labels) {
            Rectangle2D bounds = tree.get(label);
//            int dx = random.nextInt(5) * ((random.nextBoolean())? 1 : -1);
//            int dy = random.nextInt(5) * ((random.nextBoolean())? 1 : -1);
//            bounds.setRect(bounds.getMinX() + dx, bounds.getMinY() + dy, bounds.getWidth(), bounds.getHeight());
            g.setColor(color);
            g.draw(bounds);
            g.setColor(Color.RED);
            Point2f p = super.getTextPoint(g, label, bounds.getCenterX(), bounds.getCenterY());
            g.drawString(label, p.x, p.y);
        }
    }
    
    public void printBoundsTree(BoundsTree tree, List<String> labels) {
        for (String label : labels) {
            Rectangle2D r = tree.get(label);
            if (r == null) {
                System.out.println(label + ":NULL");
            } else {
                System.out.println(label + ":" + BoundsPrinter.toString(r));
            }
        }
    }

    public BoundsTree centerOnOrigin(IReaction reaction, BoundsTree boundsTree) {
        return centerOn(reaction, 0, 0, boundsTree);
    }
    
    public BoundsTree centerOn(IReaction reaction, double cx, double cy, BoundsTree boundsTree) {
        double boundsX = boundsTree.getRoot().getCenterX();
        double boundsY = boundsTree.getRoot().getCenterY();
        double dx = cx - boundsX;
        double dy = cy - boundsY;
//        System.out.println("shifting (" + boundsX + ", " + boundsY + ") by "+ dx + " " + dy);
        return shift(reaction, boundsTree, dx, dy);
    }

    public BoundsTree shift(
            IReaction reaction, BoundsTree unshiftedTree, double dx, double dy) {
        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet products = reaction.getProducts();
        String rootLabel = reaction.getID();
        
        String reactantID = reactants.getID();
        BoundsTree reactantTree = unshiftedTree.getSubtree(rootLabel + "_" + reactantID); 
        BoundsTree rBoundsTree = shift(reactants, rootLabel, reactantTree, reactantID, dx, dy);

        String productID = products.getID();
        BoundsTree productTree = unshiftedTree.getSubtree(rootLabel + "_" + productID);
        BoundsTree pBoundsTree = shift(products, rootLabel, productTree, products.getID(), dx, dy);
        
        BoundsTree boundsTree = 
            new BoundsTree(reaction.getID(), rBoundsTree, pBoundsTree);
        
        if (exactReactionLayout == null) {
            double pos = reactionLayout.getArrowPos();
            if (reactionLayout.getArrowAxis() == Axis.X) {
                reactionLayout.setArrowPos(pos + dy);
            } else {
                reactionLayout.setArrowPos(pos + dx);
            }
        } else {
            double pos = reactionLayout.getArrowPos();
            if (exactReactionLayout.getArrowAxis() == Axis.X) {
                exactReactionLayout.setArrowPos(pos + dy);
            } else {
                exactReactionLayout.setArrowPos(pos + dx);
            }
        }
        
//        Rectangle2D rBB = boundsTree.get(reaction.getID());
//        System.out.println("reaction center " + rBB.getCenterX() + " " + rBB.getCenterY());
//        System.out.println("arrow center " + reactionLayout.getArrowCenter());
        
        return boundsTree;
    }

    // XXX : seems unnecessary to pass in a subtree, when we could just look it
    // up in the main bounds tree
    private BoundsTree shift(IAtomContainerSet reactants, String rootLabel,
                             BoundsTree unshiftedMolSetTree, 
                             String label, double dx, double dy) {
        BoundsTree boundsTree = new BoundsTree(label);
        int counter = 0;
        for (IAtomContainer atomContainer : reactants.atomContainers()) {
            String fullLabel = rootLabel + "_" + label + "_" + atomContainer.getID() + ":" + counter;
            String subLabel = label + "_" + atomContainer.getID()  + ":" + counter;
//            System.out.println(fullLabel);
            
//            System.out.println("Atoms From" + BoundsPrinter.toString(GeometryTools.getRectangle2D(atomContainer)));
            GeometryTools.translate2D(atomContainer, dx, dy);
            Rectangle2D uBounds = unshiftedMolSetTree.get(fullLabel);
            Rectangle2D sBounds = new Rectangle2D.Double(uBounds.getMinX() + dx, 
                                                         uBounds.getMinY() + dy, 
                                                         uBounds.getWidth(), 
                                                         uBounds.getHeight());
//            System.out.println("From " + BoundsPrinter.toString(uBounds));
//            System.out.println("To   " + BoundsPrinter.toString(sBounds));
//            System.out.println("Atoms To" + BoundsPrinter.toString(GeometryTools.getRectangle2D(atomContainer)));
            boundsTree.add(subLabel, sBounds);
            counter++;
        }
        return boundsTree;
    }

    private Axis getAxisOfExpansion(
            int targetWidth, int targetHeight, double zoomedWidth, double zoomedHeight) {
//        int borderX = params.borderX;
//        int borderY = params.borderY;
//        int border2 = borderX * 2;
        double widthRatio = targetWidth / zoomedWidth;
        double heightRatio = targetHeight / zoomedHeight;
        if (widthRatio < heightRatio) {
            return Axis.X;
        } else {
            return Axis.Y;
        }
    }

    public void align(IReaction reaction,
            int w, int h, double actualWidth, double actualHeight, double zoom) {
        double b2X = params.borderX * 2;
        double b2Y = params.borderY * 2;
        double bXz = params.borderX * zoom;
        double bYz = params.borderY * zoom;
        double ww = (w + b2X) / zoom;
        double hh = (h + b2Y) / zoom;
    
        Axis expansionAxis = getAxisOfExpansion(w, h, actualWidth, actualHeight);
        XAlign xAlign = params.leftRightAlignment;
        YAlign yAlign = params.topBottomAlignment;
        double boundsCenterX = actualWidth / 2;
        double boundsCenterY = actualHeight / 2;
        if (expansionAxis == Axis.Y && xAlign == XAlign.LEFT) {
            return;
        } else if (expansionAxis == Axis.Y && xAlign == XAlign.CENTER) {
            double dx = (ww / 2) - (boundsCenterX + bXz);
            shift(reaction, null, dx, 0);   // FIXME
        } else if (expansionAxis == Axis.X && yAlign == YAlign.CENTER) {
            double dy = (hh / 2) - (boundsCenterY + bYz);
            shift(reaction, null, 0, dy);   // FIXME
        } else {
            // TODO
        }
    }

    public void drawReactionID(IReaction reaction, int w, Graphics2D g) {
        String id = reaction.getID();
        if (id == null) {
            return;
        }
        g.setColor(Color.BLACK);
        g.drawString(id, w / 2, 10);
    }

    public BoundsTree getReactionBounds(IReaction reaction) {
        return reactionLayout.layout(reaction, reactionLayout.getAxis());
    }
    
    public BoundsTree getExactReactionBounds(IReaction reaction, Graphics2D g) {
        BoundsTree tree = reactionLayout.layout(reaction, reactionLayout.getAxis());
        if (exactReactionLayout != null) {
            return exactReactionLayout.layout(reaction, g);
        }
        return tree;
    }

    public double calculateZoom(
            int targetWidth, int targetHeight, Rectangle2D totalBounds) {
        int borderX = params.borderX;
        int borderX2 = borderX * 2;
        int borderY = params.borderY;
        int borderY2 = borderY * 2;
        double preZoomedWidth = totalBounds.getWidth() + borderX2;
        double preZoomedHeight = totalBounds.getHeight() + borderY2;
        return calculateZoom(
                targetWidth, targetHeight, preZoomedWidth, preZoomedHeight);
    }

    private double calculateZoom(int targetWidth, int targetHeight,
            double actualWidth, double actualHeight) {
        return Math.min(targetWidth / actualWidth,
                targetHeight / actualHeight);
    }

    public Rectangle2D getDrawnBounds(List<IAtom> atoms) {
        return moleculeDrawer.getDrawnBounds(atoms);
    }

    public Vector2d getReactionAxis() {
        return reactionLayout.getAxis();
    }

    public void drawMappings(IReaction reaction, Graphics2D g) {
        g.setColor(Color.LIGHT_GRAY);
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            Point2d p0 = a0.getPoint2d();
            Point2d p1 = a1.getPoint2d();
            g.drawLine((int) p0.x, (int) p0.y, (int) p1.x, (int) p1.y);
        }
    }

    private Rectangle2D molSetBounds(IAtomContainerSet molSet) {
        Rectangle2D bounds = null;
        for (IAtomContainer ac : molSet.atomContainers()) {
            Rectangle2D currentBounds = GeometryTools.getRectangle2D(ac);
            if (bounds == null) {
                bounds = (Rectangle2D) currentBounds.clone();
            } else {
                bounds.add(currentBounds);
            }
        }
        return bounds;
    }

    public void drawMoleculeSet(IAtomContainerSet reactants, String reactionID, double yAxis, BoundsTree boundsTree, Graphics2D g) {
        if (params.layoutLeftToRight) {
            params.moleculeLabelFontSize = params.leftToRightMoleculeLabelFontSize;
        } else {
            params.moleculeLabelFontSize = params.topToBottomMoleculeLabelFontSize;
        }
        for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
            IAtomContainer current = reactants.getAtomContainer(i);
            if (i > 0) {
                IAtomContainer previous = reactants.getAtomContainer(i - 1);
                String previousLabel = reactionID + "_" + reactants.getID() + "_" + previous.getID() + ":" + (i - 1);
//                System.out.println("getting " + previousLabel);
                drawPlus(previous, previousLabel, yAxis, boundsTree, g);
                moleculeDrawer.drawMolecule(current, g);
            } else {
                moleculeDrawer.drawMolecule(current, g);
            }
        }
    }

    public void drawPlus(IAtomContainer ac, String acLabel, double yAxis, BoundsTree boundsTree, Graphics2D g) {
        int plusGap = params.plusGap;

        Rectangle2D bounds;
        if (boundsTree == null) {
            bounds = GeometryTools.getRectangle2D(ac);
        } else {
//            System.out.print("looking up " + acLabel);
            bounds = boundsTree.get(acLabel);
//            System.out.println(" found " + BoundsPrinter.toString(bounds));
//            System.out.println(boundsTree);
        }

        Rectangle2D textBounds = getTextBounds(g, "+");
        double tbW = textBounds.getWidth();
        double tbH = textBounds.getHeight();
        double halfWidth = tbW / 2;
        double halfHeight = tbH/ 2;
        
        // FIXME : T2B is bounds.getMaxY...
        double posAlongAxis = bounds.getMaxX() + plusGap + halfWidth;
        
        Point2f p = getTextPoint(g, "+", posAlongAxis, yAxis);
        g.setColor(Color.WHITE);
        g.fill(new Rectangle2D.Double(
                posAlongAxis - halfWidth, yAxis - halfHeight, tbW, tbH));
        g.setColor(Color.BLACK);
        g.setFont(new Font("ARIAL", Font.PLAIN, params.plusFontSize));
//        System.out.println("drawing plus at " + p);
        g.drawString("+", p.x, p.y);
    }

    @Override
    public Params getParams() {
        return params;
    }

    public void highlightSubgraph(List<IAtom> atoms, Color color) {
        IAtomContainer highlightContainer = new AtomContainer();
        for (IAtom atom : atoms) {
            highlightContainer.addAtom(atom);
        }
        moleculeDrawer.addHighlights(highlightContainer, color);
    }

    public DirectMoleculeDrawer getMoleculeDrawer() {
        return moleculeDrawer;
    }
}
