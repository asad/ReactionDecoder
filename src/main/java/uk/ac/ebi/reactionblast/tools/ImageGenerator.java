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
package uk.ac.ebi.reactionblast.tools;

import java.awt.Color;
import static java.awt.Color.BLACK;
import static java.awt.Color.GREEN;
import static java.awt.Color.RED;
import static java.awt.Color.WHITE;
import java.awt.Dimension;
import java.awt.Graphics2D;
import static java.awt.GraphicsEnvironment.isHeadless;
import java.awt.Image;
import static java.awt.RenderingHints.KEY_ANTIALIASING;
import static java.awt.RenderingHints.VALUE_ANTIALIAS_ON;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import static java.awt.image.BufferedImage.TYPE_4BYTE_ABGR;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import static java.lang.String.valueOf;
import static java.lang.System.err;
import static java.lang.System.setProperty;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static javax.imageio.ImageIO.write;
import javax.vecmath.Vector2d;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.FORWARD;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;
import org.openscience.smsd.AtomAtomMapping;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.DirectRBLastReactionDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.Highlighter;
import uk.ac.ebi.reactionblast.graphics.direct.OutlineHighlighter;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.RootSystem;
import static uk.ac.ebi.reactionblast.graphics.direct.SignatureRootFinder.findRootSystems;
import uk.ac.ebi.reactionblast.graphics.direct.SimpleHighlighter;
import uk.ac.ebi.reactionblast.graphics.direct.awtlayout.AbstractAWTReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.awtlayout.LeftToRightAWTReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.awtlayout.TopToBottomAWTReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.AbstractDirectReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;
import uk.ac.ebi.reactionblast.graphics.direct.layout.LeftToRightReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.TopToBottomReactionLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitGridLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import uk.ac.ebi.reactionblast.mapping.helper.RBlastReaction;
import uk.ac.ebi.reactionblast.signature.SignatureMatcher;
import static uk.ac.ebi.reactionblast.tools.LayoutCheck.getMoleculeWithLayoutCheck;

/**
 *
 * @author asad
 */
public class ImageGenerator {

    private static final ILoggingTool LOGGER
            = createLoggingTool(ImageGenerator.class);

    /**
     *
     */
    public final static int SUB_IMAGE_WIDTH = 300;

    /**
     *
     */
    public final static int SUB_IMAGE_HEIGHT = 300;

    static {
        /* works fine! ! */
 /*
         This makes the awt headless
         */

        setProperty("java.awt.headless", "true");
        LOGGER.info("Headless enabled: " + isHeadless());
        /* ---> prints true */
    }

    /**
     *
     * @param width
     * @param height
     * @return
     */
    public synchronized static Image getBlankImage(int width, int height) {
        return new BufferedImage(width, height, TYPE_4BYTE_ABGR);
    }

    /**
     *
     * @param reaction
     * @param layout
     * @param awtLayout
     * @param width
     * @param height
     * @param outFile
     * @throws IOException
     */
    protected static synchronized void makeReactionCenterHighlightedReactionToFile(
            IReaction reaction,
            AbstractDirectReactionLayout layout,
            AbstractAWTReactionLayout awtLayout,
            int width, int height,
            File outFile) throws IOException {
        Params params = new Params();

        params.leftToRightMoleculeLabelFontSize = 10;

        params.drawMappings = false;
        params.drawHighlights = true;
        params.highlightsAbove = true;

        params.drawAtomID = false;

        params.drawMoleculeID = false;
        params.drawLabelPanel = true;
        params.drawAromaticCircles = true;

        params.useCircularHighlight = false;

        params.drawSubgraphBoxes = false;
        params.drawBondStereoChanges = false;
        params.drawBondFormedCleavedMarks = true;
        params.drawBondOrderChangedMarks = true;

        params.arrowGap = 30;
        params.arrowLength = 60;
        params.drawFatArrow = true;
        params.drawArrowFilled = true;

        params.borderY = 40;

        params.drawRS = true;
        params.shouldCrop = true;

        RBlastReaction rblReaction = new RBlastReaction(reaction, true);
        Map<IAtomContainer, List<RootSystem>> rootSystems
                = findRootSystems(rblReaction);

        DirectRBLastReactionDrawer reactionDrawer
                = new DirectRBLastReactionDrawer(params, layout, awtLayout);
        Color rootColor = RED;
        Color neighbourColor = GREEN;
        DirectMoleculeDrawer moleculeDrawer
                = reactionDrawer.getReactionDrawer().getMoleculeDrawer();
        moleculeDrawer.getHighlighters().clear();   // XXX HACK
        for (IAtomContainer atomContainer : rootSystems.keySet()) {
            List<RootSystem> rootSystemList = rootSystems.get(atomContainer);
            for (RootSystem rootSystem : rootSystemList) {
                IAtomContainer rootContainer
                        = reaction.getBuilder().newInstance(IAtomContainer.class);
                rootSystem.getRoots().stream().forEach((root) -> {
                    rootContainer.addAtom(root);
                });
                IAtomContainer neighbourContainer
                        = reaction.getBuilder().newInstance(IAtomContainer.class);
                rootSystem.getLeaves().stream().forEach((leaf) -> {
                    neighbourContainer.addAtom(leaf);
                });
                Highlighter highlighter = new SimpleHighlighter(params);
                highlighter.addHighlights(rootContainer, rootColor);
                highlighter.addHighlights(neighbourContainer, neighbourColor);
                moleculeDrawer.addHighlighter(highlighter);
            }
        }

        BufferedImage image = (BufferedImage) getBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        Rectangle2D finalBounds
                = reactionDrawer.drawRBlastReaction(rblReaction, width, height, g);
        if (params.shouldCrop
                && (finalBounds.getWidth() != width
                || finalBounds.getHeight() != height)) {
            image = image.getSubimage((int) finalBounds.getX(),
                    (int) finalBounds.getY(),
                    (int) finalBounds.getWidth(),
                    (int) finalBounds.getHeight());
        }
        g.dispose();
        write(image, "PNG", outFile);
    }

    /**
     *
     * @param reaction
     * @param layout
     * @param awtLayout
     * @param width
     * @param height
     * @param shouldCrop
     * @param outFile
     * @throws IOException
     */
    protected static synchronized void makeLeftToRighHighlightedReactionToFile(
            IReaction reaction,
            AbstractDirectReactionLayout layout,
            AbstractAWTReactionLayout awtLayout,
            int width, int height,
            boolean shouldCrop,
            File outFile) throws IOException {

        RBlastReaction rblReaction = new RBlastReaction(reaction, true);
        DirectRBLastReactionDrawer drawer
                = new DirectRBLastReactionDrawer(
                        new Params(),
                        layout,
                        awtLayout);

        drawer.getParams().drawMappings = false;
        drawer.getParams().drawAromaticCircles = false;
        /*
        * set ids to false
         */
        drawer.getParams().drawAtomID = false;
        drawer.getParams().drawLonePairs = false;
        drawer.getParams().drawMoleculeID = true;
        //Make this false
        drawer.getParams().drawSubgraphBoxes = false;
        drawer.getParams().highlightSubgraphs = true;
        drawer.getParams().drawSubgraphMappingLines = false;
        drawer.getParams().highlightsBelow = false;
        drawer.getParams().highlightsAbove = true;
        drawer.getParams().drawAromaticCircles = true;
        drawer.getParams().highlightAlpha = 0.25f;
        drawer.getParams().drawRS = true;
        drawer.getParams().labelYGap = 25;
        drawer.getParams().borderY = 40;
        drawer.getParams().borderX = 40;
        drawer.getParams().arrowGap = 30;
        drawer.getParams().arrowLength = 60;
        drawer.getParams().drawArrowFilled = true;
        drawer.getParams().drawFatArrow = true;
        drawer.getParams().shouldCrop = shouldCrop;
        drawer.getParams().leftToRightMoleculeLabelFontSize = 10;

        Image drawRBlastReaction = drawer.drawRBlastReaction(rblReaction, width, height);
        write((RenderedImage) drawRBlastReaction, "PNG", outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param width
     * @param height
     * @param shouldCrop
     * @param outFile
     * @throws IOException
     */
    protected static synchronized void makeLeftToRighHighlightedReactionToFile(
            IReaction cdkReaction,
            int width, int height,
            boolean shouldCrop,
            File outFile) throws IOException {

        RBlastReaction rbReaction = new RBlastReaction(cdkReaction, true);

        DirectRBLastReactionDrawer drawer
                = new DirectRBLastReactionDrawer(new Params(),
                        new LeftToRightReactionLayout(),
                        new LeftToRightAWTReactionLayout());

        drawer.getParams().drawMappings = false;
        drawer.getParams().drawAromaticCircles = false;
        /*
        * set ids to false
         */
        drawer.getParams().drawAtomID = false;
        drawer.getParams().drawLonePairs = false;
        drawer.getParams().drawMoleculeID = true;
        //Make this false
        drawer.getParams().drawSubgraphBoxes = false;
        drawer.getParams().highlightSubgraphs = true;
        drawer.getParams().drawSubgraphMappingLines = false;
        drawer.getParams().highlightsBelow = false;
        drawer.getParams().highlightsAbove = true;
        drawer.getParams().drawAromaticCircles = true;
        drawer.getParams().highlightAlpha = 0.25f;
        drawer.getParams().drawRS = true;
        drawer.getParams().labelYGap = 25;
        drawer.getParams().borderY = 40;
        drawer.getParams().arrowGap = 30;
        drawer.getParams().arrowLength = 60;
        drawer.getParams().drawFatArrow = true;
        drawer.getParams().shouldCrop = shouldCrop;
        drawer.getParams().leftToRightMoleculeLabelFontSize = 10;

        /*
        * Hack the code to crop by Asad else use //java.awt.Image image =
        * drawer.drawRBlastReaction(rbReaction, width, height); for usual image
         */
        BufferedImage image = (BufferedImage) getBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        Rectangle2D finalBounds
                = drawer.drawRBlastReaction(rbReaction, width, height, g);
        if (shouldCrop
                && (finalBounds.getWidth() != width
                || finalBounds.getHeight() != height)) {
            image = image.getSubimage((int) finalBounds.getX(),
                    (int) finalBounds.getY(),
                    (int) finalBounds.getWidth(),
                    (int) finalBounds.getHeight());
        }
        g.dispose();

        write(image, "PNG", outFile);

    }

    /**
     *
     * @param cdkReaction
     * @param width
     * @param height
     * @param outFile
     * @throws IOException
     */
    protected static synchronized void makeTopToBottomRHighlightedReactionToFile(
            IReaction cdkReaction,
            int width, int height,
            File outFile) throws IOException {

        RBlastReaction rbReaction = new RBlastReaction(cdkReaction, true);

        DirectRBLastReactionDrawer drawer
                = new DirectRBLastReactionDrawer(
                        new Params(),
                        new TopToBottomReactionLayout(),
                        new TopToBottomAWTReactionLayout());
        drawer.getParams().drawMappings = false;
        drawer.getParams().drawAromaticCircles = false;
        drawer.getParams().drawAtomID = true;
        drawer.getParams().drawLonePairs = false;
        //Make this false
        drawer.getParams().drawSubgraphBoxes = false;
        drawer.getParams().highlightSubgraphs = true;
        drawer.getParams().drawSubgraphMappingLines = false;
        drawer.getParams().highlightsBelow = false;
        drawer.getParams().highlightsAbove = true;
        drawer.getParams().drawAromaticCircles = true;
        drawer.getParams().highlightAlpha = 0.25f;
        drawer.getParams().drawRS = true;
        drawer.getParams().labelYGap = 25;
        drawer.getParams().borderY = 40;
        drawer.getParams().arrowGap = 30;
        drawer.getParams().arrowLength = 60;
        drawer.getParams().drawFatArrow = true;
        drawer.getParams().drawArrowFilled = true;
        drawer.getParams().drawLabelPanel = false;
        drawer.getParams().drawMoleculeID = true;
        drawer.getParams().topToBottomMoleculeLabelFontSize = 10;

        java.awt.Image image = drawer.drawRBlastReaction(rbReaction, width, height);
        write((RenderedImage) image, "PNG", outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void LeftToRightReactionLayoutImageSmall(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {
        int width = 600;
        int height = 400;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeLeftToRighHighlightedReactionToFile(cdkReaction, width, height, true, outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void LeftToRightReactionCenterImageSmall(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {
        int width = 600;
        int height = 400;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeReactionCenterHighlightedReactionToFile(cdkReaction,
                new LeftToRightReactionLayout(),
                new LeftToRightAWTReactionLayout(), width, height, outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void TopToBottomReactionLayoutImageSmall(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {

        int height = 400;
        int width = 600;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeTopToBottomRHighlightedReactionToFile(cdkReaction, width, height, outFile);

    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void LeftToRightReactionLayoutImage(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {
        int height = 800;
        int width = 1200;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeLeftToRighHighlightedReactionToFile(cdkReaction, width, height, false, outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void LeftToRightReactionCenterImage(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {
        int height = 800;
        int width = 1200;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeReactionCenterHighlightedReactionToFile(cdkReaction,
                new LeftToRightReactionLayout(),
                new LeftToRightAWTReactionLayout(), width, height, outFile);
    }

    /**
     *
     * @param cdkReaction
     * @param rmrID
     * @param outputDir
     * @throws Exception
     */
    public synchronized static void TopToBottomReactionLayoutImage(
            IReaction cdkReaction, String rmrID, String outputDir) throws Exception {
        int height = 800;
        int width = 1200;
        File outFile = new File(getDir(outputDir), rmrID + ".png");
        makeTopToBottomRHighlightedReactionToFile(cdkReaction, width, height, outFile);
    }

    private synchronized static File getDir(String outputDir) {
        File file = new File(outputDir);
        if (!file.exists()) {
            boolean success = file.mkdirs();
            if (!success) {
                LOGGER.debug("Could not make dir " + file);
            }
        }
        return file;
    }
    private final List<QueryTargetPair> queryTargetPairs;
    private final Params params;

    /**
     *
     */
    public ImageGenerator() {
        queryTargetPairs = new ArrayList<>();
        params = new Params();
    }

    private IReaction layoutReaction(
            IReaction mappedReaction, String reactionID) {
        IReaction reactionWithLayout = new Reaction();
        reactionWithLayout.setDirection(FORWARD);
        reactionWithLayout.setID(reactionID);

        for (IAtomContainer ac : mappedReaction.getReactants().atomContainers()) {
            IAtomContainer moleculeWithLayoutCheck = getMoleculeWithLayoutCheck(ac);
            moleculeWithLayoutCheck.setID(ac.getID());
            reactionWithLayout.addReactant(ac, mappedReaction.getReactantCoefficient(ac));
        }
        for (IAtomContainer ac : mappedReaction.getProducts().atomContainers()) {
            IAtomContainer moleculeWithLayoutCheck = getMoleculeWithLayoutCheck(ac);
            moleculeWithLayoutCheck.setID(ac.getID());
            reactionWithLayout.addProduct(ac, mappedReaction.getProductCoefficient(ac));
        }

        for (IMapping m : mappedReaction.mappings()) {
            reactionWithLayout.addMapping(m);
        }

        reactionWithLayout.setFlags(mappedReaction.getFlags());
        reactionWithLayout.setProperties(mappedReaction.getProperties());
        return reactionWithLayout;

    }

    /**
     *
     * @param query
     * @param target
     * @param label
     * @param maxac
     * @throws IOException
     * @throws Exception
     */
    public void addImages(
            IAtomContainer query, IAtomContainer target, String label, Map<Integer, Integer> maxac) throws IOException, Exception {

        SingleMoleculeLayout msl = new SingleMoleculeLayout(params);
        msl.layout(query, new Vector2d(0.0, 0.0));
        msl.layout(target, new Vector2d(0.0, 0.0));

        IAtomContainer cloneOfQuery = new org.openscience.cdk.AtomContainer(query).clone();
        IAtomContainer cloneOfTarget = new org.openscience.cdk.AtomContainer(target).clone();

        IAtomContainer querySubgraph = query.getBuilder().newInstance(IAtomContainer.class, cloneOfQuery);
        IAtomContainer targetSubgraph = target.getBuilder().newInstance(IAtomContainer.class, cloneOfTarget);
        List<IAtom> n1 = new ArrayList<>(query.getAtomCount());
        List<IAtom> n2 = new ArrayList<>(target.getAtomCount());

        for (Map.Entry<Integer, Integer> aMaps : maxac.entrySet()) {
            IAtom qAtom = cloneOfQuery.getAtom(aMaps.getKey());
            IAtom tAtom = cloneOfTarget.getAtom(aMaps.getValue());
            qAtom.setID(aMaps.getKey().toString());
            tAtom.setID(aMaps.getValue().toString());
            n1.add(qAtom);
            n2.add(tAtom);
        }

        for (IAtom atom : cloneOfQuery.atoms()) {
            if (!n1.contains(atom)) {
                querySubgraph.removeAtom(atom);
            }
        }

        for (IAtom atom : cloneOfTarget.atoms()) {
            if (!n2.contains(atom)) {
                targetSubgraph.removeAtom(atom);
            }
        }

        queryTargetPairs.add(
                new QueryTargetPair(
                        cloneOfQuery, cloneOfTarget, querySubgraph, targetSubgraph, label));
    }

    /**
     *
     * @param query
     * @param target
     * @param label
     * @param maxac
     * @throws IOException
     * @throws Exception
     */
    public void addImages(IAtomContainer query, IAtomContainer target, String label, AtomAtomMapping maxac) throws IOException, Exception {

        SingleMoleculeLayout msl = new SingleMoleculeLayout(params);
        msl.layout(query, new Vector2d(0.0, 0.0));
        msl.layout(target, new Vector2d(0.0, 0.0));

        IAtomContainer cloneOfQuery = new org.openscience.cdk.AtomContainer(query).clone();
        IAtomContainer cloneOfTarget = new org.openscience.cdk.AtomContainer(target).clone();

        IAtomContainer querySubgraph = query.getBuilder().newInstance(IAtomContainer.class, cloneOfQuery);
        IAtomContainer targetSubgraph = target.getBuilder().newInstance(IAtomContainer.class, cloneOfTarget);
        List<IAtom> n1 = new ArrayList<>(query.getAtomCount());
        List<IAtom> n2 = new ArrayList<>(target.getAtomCount());

        for (Map.Entry<IAtom, IAtom> aMaps : maxac.getMappingsByAtoms().entrySet()) {
            IAtom qAtom = aMaps.getKey();
            IAtom tAtom = aMaps.getValue();
            qAtom.setID(valueOf(maxac.getQueryIndex(qAtom)));
            tAtom.setID(valueOf(maxac.getQueryIndex(tAtom)));
            n1.add(qAtom);
            n2.add(tAtom);
        }

        for (IAtom atom : cloneOfQuery.atoms()) {
            if (!n1.contains(atom)) {
                querySubgraph.removeAtom(atom);
            }
        }

        for (IAtom atom : cloneOfTarget.atoms()) {
            if (!n2.contains(atom)) {
                targetSubgraph.removeAtom(atom);
            }
        }

        queryTargetPairs.add(
                new QueryTargetPair(
                        cloneOfQuery, cloneOfTarget, querySubgraph, targetSubgraph, label));
    }

    /**
     *
     * @param outImageFileName
     * @param qName
     * @param tName
     */
    public void createImage(String outImageFileName, String qName, String tName) {

        // layout, and set the highlight subgraphs
        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        IChemObjectBuilder builder = getInstance();
        IAtomContainerSet leftHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        IAtomContainerSet rightHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        queryTargetPairs.stream().map((pair) -> {
            moleculeDrawer.addHighlights(pair.querySubgraph);
            return pair;
        }).map((pair) -> {
            moleculeDrawer.addHighlights(pair.targetSubgraph);
            return pair;
        }).map((pair) -> {
            leftHandMoleculeSet.addAtomContainer(pair.query);
            return pair;
        }).forEach((pair) -> {
            rightHandMoleculeSet.addAtomContainer(pair.target);
        });

        // calculate the total dimensions of the final image
        int width = SUB_IMAGE_WIDTH * 2;
        int height = SUB_IMAGE_HEIGHT * queryTargetPairs.size();

        // make the image, and draw the molecules on it
        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        List<IAtomContainer> mols = new ArrayList<>();
        queryTargetPairs.stream().map((pair) -> {
            mols.add(pair.query);
            return pair;
        }).forEach((pair) -> {
            mols.add(pair.target);
        });
        ZoomToFitGridLayout layoutDrawer = new ZoomToFitGridLayout(moleculeDrawer, queryTargetPairs.size(), 2);
        layoutDrawer.layout(mols, new Dimension(SUB_IMAGE_WIDTH, SUB_IMAGE_HEIGHT), g);

        float labelX = SUB_IMAGE_WIDTH / 2;
        float labelY = 15;
        g.setColor(BLACK);
        for (QueryTargetPair pair : queryTargetPairs) {
            g.drawString(pair.label, labelX, labelY);
            labelY += SUB_IMAGE_HEIGHT;
        }
        g.dispose();

        try {
            write((RenderedImage) image, "PNG", new File(outImageFileName + ".png"));
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    /**
     *
     * @return
     */
    public RenderedImage createImage() {

        // layout, and set the highlight subgraphs
        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        IChemObjectBuilder builder = getInstance();
        IAtomContainerSet leftHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        IAtomContainerSet rightHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        queryTargetPairs.stream().map((pair) -> {
            moleculeDrawer.addHighlights(pair.querySubgraph);
            return pair;
        }).map((pair) -> {
            moleculeDrawer.addHighlights(pair.targetSubgraph);
            return pair;
        }).map((pair) -> {
            leftHandMoleculeSet.addAtomContainer(pair.query);
            return pair;
        }).forEach((pair) -> {
            rightHandMoleculeSet.addAtomContainer(pair.target);
        });

        // calculate the total dimensions of the final image
        int width = SUB_IMAGE_WIDTH * 2;
        int height = SUB_IMAGE_HEIGHT * queryTargetPairs.size();

        // make the image, and draw the molecules on it
        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        List<IAtomContainer> mols = new ArrayList<>();
        queryTargetPairs.stream().map((pair) -> {
            mols.add(pair.query);
            return pair;
        }).forEach((pair) -> {
            mols.add(pair.target);
        });
        ZoomToFitGridLayout layoutDrawer = new ZoomToFitGridLayout(moleculeDrawer, queryTargetPairs.size(), 2);
        layoutDrawer.layout(mols, new Dimension(SUB_IMAGE_WIDTH, SUB_IMAGE_HEIGHT), g);

        float labelX = SUB_IMAGE_WIDTH / 2;
        float labelY = 15;
        g.setColor(BLACK);
        for (QueryTargetPair pair : queryTargetPairs) {
            g.drawString(pair.label, labelX, labelY);
            labelY += SUB_IMAGE_HEIGHT;
        }
        g.dispose();

        return (RenderedImage) image;

    }

    /**
     *
     * @param outputDirName
     * @param molecule
     * @param molID
     * @throws IOException
     */
    public synchronized void directMoleculeImageNaturalScale(File outputDirName, IAtomContainer molecule, String molID) throws IOException {

        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        Params p1 = moleculeDrawer.getParams();
        p1.drawAtomID = false;
        p1.atomSymbolFontSize = 14;
        p1.moleculeLabelFontSize = 14;
        p1.bondLength = 50;
        p1.filledWedgeWidth = 10;
        p1.dashedWedgeStroke = 1.5f;
        p1.dashedGapFactor = 0.1;
        p1.drawMoleculeID = true;

        SingleMoleculeLayout layout = new SingleMoleculeLayout(p1, false);
        Vector2d center = new Vector2d(0, 0);
        layout.invert(molecule);
        BoundsTree boundsTree = layout.layout(molecule, center);
        // TODO : double-valued width and height
        Dimension canvasSize = new Dimension((int) boundsTree.getWidth(), (int) boundsTree.getHeight());
        int width = canvasSize.width + (2 * p1.borderX);
        int height = canvasSize.height + (2 * p1.borderY);
        int centerX = width / 2;
        int centerY;
        if (p1.drawMoleculeID) {
            centerY = (int) ((height / 2) - p1.labelYGap);
        } else {
            centerY = height / 2;
        }
//        layout.translateTo(molecule, centerX, centerY);   // TODO check and FIXME
        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        moleculeDrawer.drawMolecule(molecule, g);
        File outFile = new File(outputDirName, molID + ".png");
        write((RenderedImage) image, "PNG", outFile);
    }

    /**
     *
     * @param outputDirName
     * @param molecule
     * @param molID
     * @throws IOException
     */
    public synchronized void directMoleculeImageZoomedToFit(
            File outputDirName, IAtomContainer molecule, String molID) throws IOException {
        int width = 800;
        int height = 600;
        directMoleculeImageZoomedToFit(outputDirName, molecule, molID, width, height);
    }

    /**
     *
     * @param outputDirName
     * @param molecule
     * @param molID
     * @param width
     * @param height
     * @throws IOException
     */
    public synchronized void directMoleculeImageZoomedToFit(File outputDirName, IAtomContainer molecule, String molID, int width, int height) throws IOException {

        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        Params par = moleculeDrawer.getParams();
        par.drawAtomID = false;
        par.bondLength = 50;
        par.moleculeLabelFontSize = 20;
        par.atomSymbolFontSize = 20;
        par.labelYGap = 30.0;
        par.filledWedgeWidth = 10;
        par.dashedWedgeStroke = 1.0f;
        par.dashedGapFactor = 0.1;
        par.borderX = 40;
        par.borderY = 40;
        par.drawImplicitHydrogens = true;

        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        ZoomToFitLayout layout = new ZoomToFitLayout(moleculeDrawer);
        layout.invert(molecule);
        layout.layout(molecule, new Dimension(width, height), g);
        File outFile = new File(outputDirName, molID + ".png");
        write((RenderedImage) image, "PNG", outFile);
    }

    /**
     *
     * @param outputDir
     * @param cdkReaction
     * @param rmrID
     * @throws Exception
     */
    public synchronized void drawTopToBottomReactionLayout(String outputDir, IReaction cdkReaction, String rmrID) throws Exception {
        drawTopToBottomReactionLayout(new CreateDirectory().createDirectory(outputDir, false), cdkReaction, rmrID);
    }

    /**
     *
     * @param outputDirName
     * @param cdkReaction
     * @param rmrID
     * @throws Exception
     */
    public synchronized void drawTopToBottomReactionLayout(File outputDirName, IReaction cdkReaction, String rmrID) throws Exception {
        int width = 800;
        int height = 1000;

        RBlastReaction rbReaction = new RBlastReaction(cdkReaction, true);

        DirectRBLastReactionDrawer drawer
                = new DirectRBLastReactionDrawer(
                        new Params(),
                        new TopToBottomReactionLayout(),
                        new TopToBottomAWTReactionLayout());
        drawer.getParams().drawMappings = false;
        drawer.getParams().drawAromaticCircles = true;
        drawer.getParams().drawAtomID = true;
        drawer.getParams().drawLonePairs = false;
        //Make this false
        drawer.getParams().drawSubgraphBoxes = false;
        drawer.getParams().highlightSubgraphs = true;
        drawer.getParams().drawSubgraphMappingLines = false;
        drawer.getParams().highlightsBelow = false;
        drawer.getParams().highlightsAbove = true;
        drawer.getParams().drawAromaticCircles = true;
        drawer.getParams().highlightAlpha = 0.25f;
        drawer.getParams().drawRS = true;
        drawer.getParams().leftToRightMoleculeLabelFontSize = 10;
        drawer.getParams().labelYGap = 25;
        drawer.getParams().borderY = 40;
        drawer.getParams().arrowGap = 30;
        drawer.getParams().arrowLength = 60;
        drawer.getParams().drawFatArrow = true;
        drawer.getParams().drawArrowFilled = true;

        Image image = drawer.drawRBlastReaction(rbReaction, width, height);
        File outFile = new File(outputDirName, rmrID + ".png");
        write((RenderedImage) image, "PNG", outFile);
    }

    /**
     *
     * @param outputDirName
     * @param cdkReaction
     * @param rmrID
     * @throws Exception
     */
    public synchronized void drawLeftToRightReactionLayout(String outputDirName, IReaction cdkReaction, String rmrID) throws Exception {
        drawLeftToRightReactionLayout(new CreateDirectory().createDirectory(outputDirName, false), cdkReaction, rmrID);
    }

    /**
     *
     * @param outputDirName
     * @param mappedReaction
     * @param reactionID
     * @throws Exception
     */
    public synchronized void drawLeftToRightReactionLayout(
            File outputDirName, IReaction mappedReaction, String reactionID) throws Exception {
        int width = 2048;
        int height = 600;

        /*
         Layout reaction to avoid image errors
         */
        IReaction reactionWithLayout = layoutReaction(mappedReaction, reactionID);
        RBlastReaction rbReaction = new RBlastReaction(reactionWithLayout, true);

        DirectRBLastReactionDrawer drawer
                = new DirectRBLastReactionDrawer(
                        new Params(),
                        new LeftToRightReactionLayout(),
                        new LeftToRightAWTReactionLayout());
        drawer.getParams().drawMappings = false;
        drawer.getParams().drawAtomID = false;
        drawer.getParams().drawLonePairs = false;
        //Make this false
        drawer.getParams().drawSubgraphBoxes = false;
        drawer.getParams().highlightSubgraphs = true;
        drawer.getParams().drawSubgraphMappingLines = false;
        drawer.getParams().highlightsBelow = false;
        drawer.getParams().highlightsAbove = true;
        drawer.getParams().drawAromaticCircles = true;

        drawer.getParams().drawRS = true;
        drawer.getParams().leftToRightMoleculeLabelFontSize = 10;
        drawer.getParams().labelYGap = 25;
        drawer.getParams().borderY = 40;
        drawer.getParams().arrowGap = 30;
        drawer.getParams().arrowLength = 60;
        drawer.getParams().drawFatArrow = true;
        drawer.getParams().drawArrowFilled = true;

        /*
         * For Lighter images
         *   drawer.getParams().highlightAlpha = 0.25f;
         *   drawer.getParams().bondStrokeWidth = default;
         */
 /* for darker presentation images
         * drawer.getParams().highlightAlpha = 0.30f;
         * drawer.getParams().bondStrokeWidth=2.0f;
         */
        drawer.getParams().highlightAlpha = 0.30f;
        drawer.getParams().bondStrokeWidth = 2.0f;

        Image image = drawer.drawRBlastReaction(rbReaction, width, height);
        File outFile = new File(outputDirName, reactionID + ".png");
        write((RenderedImage) image, "PNG", outFile);
    }

    /**
     * Draw a left-to-right reaction with atoms highlighted as reaction centers.
     * Files are drawn as "$ouputDir/$reactionID.png".
     *
     * @param mappedReaction mapped reaction
     * @param reactionCenterSigs the reaction center signatures
     * @param reactionID the id of the reaction
     * @param outputDir the directory to put the files in
     * @throws IOException if there is a problem writing the file
     */
    public void drawLeftToRightReactionCenterMoleculeImage(File outputDir,
            IReaction mappedReaction, List<String> reactionCenterSigs, String reactionID) throws IOException {
        int width = 1000;
        int height = 800;

        Params par = new Params();
        par.drawMappings = false;
        par.drawHighlights = true;
        par.highlightsAbove = true;

        par.drawAtomID = false;

        par.drawMoleculeID = true;
        par.labelYGap *= 2;
        par.drawAromaticCircles = true;

        par.useCircularHighlight = true;
        par.circularHighlightIsConcentric = true;
        par.circularHighlightTransparentFilled = false;
        par.circularHighlightShowAtoms = false;

        par.drawSubgraphBoxes = false;
        par.drawBondStereoChanges = false;
        par.drawBondFormedCleavedMarks = false;
        par.drawBondOrderChangedMarks = false;

        par.drawRS = true;

        /*
         Layout reaction to avoid image errors
         */
        IReaction reactionWithLayout = layoutReaction(mappedReaction, reactionID);
        DirectRBLastReactionDrawer reactionDrawer = new DirectRBLastReactionDrawer(
                par, new LeftToRightReactionLayout());

        RBlastReaction rblReaction = new RBlastReaction(reactionWithLayout, true);
        setHighightsFromSignatures(reactionDrawer, rblReaction, reactionCenterSigs);

        // TODO : better
        Image image = reactionDrawer.getReactionDrawer().getMoleculeDrawer().makeBlankImage(width, height);

        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        reactionDrawer.drawRBlastReaction(rblReaction, width, height, g);
        g.dispose();
        File file = new File(outputDir, reactionID + ".png");
        write((RenderedImage) image, "PNG", file);
    }

    private void setHighightsFromSignatures(DirectRBLastReactionDrawer drawer, RBlastReaction rblReaction, List<String> signatures) {
        /*
         Layout reaction to avoid image errors
         */
        IReaction reaction = layoutReaction(rblReaction.getReaction(), "Signature_" + rblReaction.getReaction().getID());

        /*
         Layout reaction to avoid image errors
         */
        // TODO : do this in one step
        SignatureMatcher matcher = new SignatureMatcher(1, 3);
        List<IAtom> roots = matcher.getMatchingRootAtoms(signatures, reaction);
        List<IAtom> filteredRoots = new ArrayList<>();
        filterByBonds(roots, filteredRoots, rblReaction.getBondsCleavedInReactant());
        filterByBonds(roots, filteredRoots, rblReaction.getBondsFormedInProduct());
        filterByBonds(roots, filteredRoots, rblReaction.getBondsOrderChangedInReactant());
        filterByBonds(roots, filteredRoots, rblReaction.getBondsOrderChangedInProduct());
        filterByAtoms(roots, filteredRoots, rblReaction.getAtomStereoProductMap().keySet());
        filterByAtoms(roots, filteredRoots, rblReaction.getAtomStereoReactantMap().keySet());

        Map<IAtom, IAtomContainer> atomToAtomContainerMap = getAtomToAtomContainerMap(reaction);
        Color rootColor = RED;
        Color neighbourColor = GREEN;
        DirectMoleculeDrawer moleculeDrawer = drawer.getReactionDrawer().getMoleculeDrawer();
        moleculeDrawer.getHighlighters().clear();   // XXX HACK
        for (IAtom atom : filteredRoots) {
            IAtomContainer rootContainer = reaction.getBuilder().newInstance(IAtomContainer.class);
            IAtomContainer atomContainer = atomToAtomContainerMap.get(atom);
            rootContainer.addAtom(atom);
//            System.out.println("root " + atom.getID() + " " + atomContainer.getID());

            IAtomContainer neighbourContainer = reaction.getBuilder().newInstance(IAtomContainer.class);
            atomContainer.getConnectedAtomsList(atom).stream().forEach((neighbour) -> {
                neighbourContainer.addAtom(neighbour);
//                System.out.println("neighbour " + neighbour.getID());
            });
            Highlighter highlighter = new OutlineHighlighter(moleculeDrawer.getParams());
            highlighter.addHighlights(rootContainer, rootColor);
            highlighter.addHighlights(neighbourContainer, neighbourColor);
            moleculeDrawer.addHighlighter(highlighter);
        }
    }

    private Map<IAtom, IAtomContainer> getAtomToAtomContainerMap(IReaction reaction) {
        Map<IAtom, IAtomContainer> map = new HashMap<>();
        getAllAtomContainers(reaction).stream().forEach((atomContainer) -> {
            for (IAtom atom : atomContainer.atoms()) {
                map.put(atom, atomContainer);
            }
        });
        return map;
    }

    private void filterByBonds(List<IAtom> atomList, List<IAtom> filteredList, List<IBond> bonds) {
        atomList.stream().forEach((IAtom atom) -> {
            bonds.stream().filter((bond) -> (bond.contains(atom))).forEach((_item) -> {
                filteredList.add(atom);
            });
        });
    }

    private void filterByAtoms(List<IAtom> atomList, List<IAtom> filteredList, Set<IAtom> validAtoms) {
        Set setDiff = new HashSet(atomList);
        setDiff.retainAll(validAtoms);
        filteredList.addAll(setDiff);
    }

    private class QueryTargetPair {

        public final IAtomContainer query;
        public final IAtomContainer target;
        public final IAtomContainer querySubgraph;
        public final IAtomContainer targetSubgraph;
        public final String label;

        QueryTargetPair(IAtomContainer query, IAtomContainer target, IAtomContainer querySubgraph, IAtomContainer targetSubgraph, String label) {
            this.query = query;
            this.target = target;
            this.querySubgraph = querySubgraph;
            this.targetSubgraph = targetSubgraph;
            this.label = label;
        }
    }

}
