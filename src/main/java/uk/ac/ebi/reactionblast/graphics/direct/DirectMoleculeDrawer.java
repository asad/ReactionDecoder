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
import java.awt.Font;
import static java.awt.Font.PLAIN;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import javax.vecmath.Point2f;
import org.openscience.cdk.exception.Intractable;
import static org.openscience.cdk.geometry.GeometryTools.getRectangle2D;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;

/**
 *
 * @author asad
 */
public class DirectMoleculeDrawer extends AbstractDirectDrawer {

    private final static ILoggingTool LOGGER
            = createLoggingTool(DirectMoleculeDrawer.class);

    private Font moleculeIDFont;
    private List<Highlighter> highlightDrawers;
    private LabelManager labelManager;
    private DirectAtomDrawer atomDrawer;
    private DirectBondDrawer bondDrawer;
    private Map<IAtom, IStereoAndConformation> chiralMap;

    /**
     *
     * @param params
     */
    public DirectMoleculeDrawer(Params params) {
        setParams(params);
        params.bondLength = 20;

        // make an initial highlight drawer and add to a list
        highlightDrawers = new ArrayList<>();
        Highlighter highlightDrawer;
        if (params.useCircularHighlight) {
            highlightDrawer = new OutlineHighlighter(params);
        } else {
            highlightDrawer = new SimpleHighlighter(params);
        }
        highlightDrawers.add(highlightDrawer);

        labelManager = new LabelManager();
        atomDrawer = new DirectAtomDrawer(params, labelManager);
        bondDrawer = new DirectBondDrawer(params, labelManager);
        chiralMap = new HashMap<>();
    }

    /**
     *
     */
    public DirectMoleculeDrawer() {
        this(new Params());
    }

    /**
     *
     * @param chirals
     */
    public void addToChiralMap(Map<IAtom, IStereoAndConformation> chirals) {
        chiralMap.putAll(chirals);
    }

    /**
     *
     * @param atoms
     * @return
     */
    public Rectangle2D getDrawnBounds(List<IAtom> atoms) {
        return atomDrawer.getDrawnBounds(atoms);
    }

    /**
     * Removes all the highlights from the drawer.
     */
    public void clearHighlights() {
        for (Highlighter highlightDrawer : highlightDrawers) {
            highlightDrawer.clearHighlights();
        }
    }

    /**
     * Get the first highlighter in the list, or create one if none exists.
     *
     * @return a highlighter
     */
    public Highlighter getFirstHighlighter() {
        Highlighter highlightDrawer;
        if (highlightDrawers.isEmpty()) {
            if (params.useCircularHighlight) {
                highlightDrawer = new OutlineHighlighter(params);
            } else {
                highlightDrawer = new SimpleHighlighter(params);
            }
            highlightDrawers.add(highlightDrawer);
        } else {
            highlightDrawer = highlightDrawers.get(0);
        }
        return highlightDrawer;
    }

    /**
     * Get the list of highlighters.
     *
     * @return a reference to the list of highlight drawers
     */
    public List<Highlighter> getHighlighters() {
        return highlightDrawers;
    }

    /**
     * Add a highlighter to the list.
     *
     * @param highlighter a class implementing the highlighter interface
     */
    public void addHighlighter(Highlighter highlighter) {
        this.highlightDrawers.add(highlighter);
    }

    /**
     * Set the highlights for all atoms and bonds in the highlight container to
     * this color.
     *
     * @param highlightContainer
     * @param color
     */
    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        Highlighter highlightDrawer = getFirstHighlighter();
        highlightDrawer.addHighlights(highlightContainer, color);
    }

    /**
     *
     * @param atoms
     * @param color
     */
    public void addHighlights(List<IAtom> atoms, Color color) {
        Map<IAtom, Color> atomColorMap = new HashMap<>();
        atoms.forEach((atom) -> {
            atomColorMap.put(atom, color);
        });
        Highlighter highlightDrawer = getFirstHighlighter();
        highlightDrawer.addToHighlights(atomColorMap);
    }

    /**
     * Set the highlights for all atoms and bonds in the container to the color
     * set in Params.highlightColor.
     *
     * @param highlightContainer
     */
    public void addHighlights(IAtomContainer highlightContainer) {
        addHighlights(highlightContainer, params.highlightColor);
    }

    /**
     * Set the highlights for all the atoms and bonds to the color in
     * Params.highlightColor.
     *
     * @param atoms
     * @param bonds
     */
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        Highlighter highlightDrawer = getFirstHighlighter();
        highlightDrawer.addHighlights(atoms, bonds);
    }

    /**
     * Set the highlights for all the atoms in the list to the color in
     * Params.highlightColor.
     *
     * @param atoms
     */
    public void addHighlights(List<IAtom> atoms) {
        addHighlights(atoms, new ArrayList<>());
    }

    /**
     *
     * @param colorMap
     */
    public void addToHighlights(Map<IAtom, Color> colorMap) {
        Highlighter highlightDrawer = getFirstHighlighter();
        highlightDrawer.addToHighlights(colorMap);
    }

    /**
     *
     * @param molecule
     * @param g
     */
    public void drawMolecule(IAtomContainer molecule, Graphics2D g) {
        // reset label manager
        labelManager.reset();

        // setup fonts
        atomDrawer.setAtomSymbolFont(new Font("ROMAN", PLAIN, params.atomSymbolFontSize));
        atomDrawer.setSubscriptFont(new Font("ROMAN", PLAIN, params.subscriptTextSize));
        atomDrawer.setAtomIDFont(new Font("ROMAN", PLAIN, params.atomIDFontSize));
        atomDrawer.setChiralSymbolFont(new Font("ROMAN", PLAIN, params.chiralSymbolFontSize));

        moleculeIDFont = new Font("ROMAN", PLAIN, params.moleculeLabelFontSize);

        Color savedColor = g.getColor();
        if (params.drawBounds) {
            Rectangle2D bounds = getRectangle2D(molecule);
            g.draw(bounds);
        }

        if (params.drawHighlights && params.highlightsBelow) {
            drawHighlights(molecule, g);
        }

        atomDrawer.setChirals(chiralMap);
        try {
            bondDrawer.drawBonds(molecule, g);
        } catch (Intractable ex) {
            LOGGER.error(Level.SEVERE, null, ex);
        }
        atomDrawer.drawAtoms(molecule, g);

        if (params.drawHighlights && params.highlightsAbove) {
            drawHighlights(molecule, g);
        }

        if (params.drawMoleculeID) {
            drawMoleculeID(molecule, g);
        }

        g.setColor(savedColor);
    }

    private void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        highlightDrawers.forEach((highlightDrawer) -> {
            highlightDrawer.drawHighlights(molecule, g);
        });
    }

    /**
     *
     * @param atomContainer
     * @param g
     * @return
     */
    public Rectangle2D drawMoleculeID(IAtomContainer atomContainer, Graphics2D g) {
        String id = atomContainer.getID();
        if (id == null) {
            return null;
        }
        Rectangle2D moleculeBounds = getRectangle2D(atomContainer);
        double labelCenterX = moleculeBounds.getCenterX();
        double labelCenterY = moleculeBounds.getMaxY() + params.labelYGap;
        Point2f textPoint = getTextPoint(g, id, labelCenterX, labelCenterY);
        g.setFont(moleculeIDFont);
        g.setColor(BLACK);
        g.drawString(id, textPoint.x, textPoint.y);
        Rectangle2D textBounds = getTextBounds(g, id);
        return new Rectangle2D.Double(
                labelCenterX - (textBounds.getWidth() / 2),
                labelCenterY - (textBounds.getHeight() / 2),
                textBounds.getWidth(),
                textBounds.getHeight());
    }
}
