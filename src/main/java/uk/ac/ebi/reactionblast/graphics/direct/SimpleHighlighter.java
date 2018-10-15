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
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author asad
 */
public class SimpleHighlighter extends AbstractHighlightDrawer implements Highlighter {

    private Map<IAtom, Color> atomColorMap;
    private final Map<IBond, Color> bondColorMap;

    /**
     *
     * @param params
     */
    public SimpleHighlighter(Params params) {
        super(params);
        atomColorMap = new HashMap<>();
        bondColorMap = new HashMap<>();
    }

    /**
     *
     * @param molecule
     * @param g
     */
    @Override
    public void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        for (IAtom atom : atomColorMap.keySet()) {
            if (molecule.contains(atom)) {
                Color color = atomColorMap.get(atom);
                drawHighlight(atom, color, g);
            }
        }
        for (IBond bond : bondColorMap.keySet()) {
            if (molecule.contains(bond)) {
                Color color = bondColorMap.get(bond);
                drawHighlight(bond, color, g);
            }
        }
    }

    /**
     * Set the highlights for all atoms and bonds in the highlight container to
     * this color.
     *
     * @param highlightContainer
     * @param color
     */
    @Override
    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        registerColor(color);
        for (IAtom atom : highlightContainer.atoms()) {
            atomColorMap.put(atom, color);
        }
        for (IBond bond : highlightContainer.bonds()) {
            bondColorMap.put(bond, color);
        }
    }

    /**
     * Set the highlights for all the atoms and bonds to the color in
     * Params.highlightColor.
     *
     * @param atoms
     * @param bonds
     */
    @Override
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        for (IAtom atom : atoms) {
            atomColorMap.put(atom, params.highlightColor);
        }
        for (IBond bond : bonds) {
            bondColorMap.put(bond, params.highlightColor);
        }
    }

    /**
     * Add the set of atom-to-color mappings to the highlights.
     *
     * @param atomColorMap
     */
    @Override
    public void addToHighlights(Map<IAtom, Color> atomColorMap) {
        this.atomColorMap.putAll(atomColorMap);
    }

    /**
     * Reset all highlights to this map of atoms to colors.
     *
     * @param atomColorMap
     */
    public void setHighlights(Map<IAtom, Color> atomColorMap) {
        this.atomColorMap = atomColorMap;
    }

    /**
     *
     * @param atom
     * @param g
     */
    public void drawHighlight(IAtom atom, Graphics2D g) {
        if (params.highlightsAbove) {
            drawHighlight(atom, translucentHighlightColor, g);
        } else {
            drawHighlight(atom, opaqueHighlightColor, g);
        }
    }

    /**
     *
     * @param atom
     * @param color
     * @param g
     */
    public void drawHighlight(IAtom atom, Color color, Graphics2D g) {
        Color actualColor;
        if (params.highlightsAbove) {
            actualColor = getTranslucentColor(color);
        } else {
            actualColor = color;
        }
        g.setColor(actualColor);
        double r = params.highlightRadius;
        double d = r * 2;
        Point2d p = atom.getPoint2d();
        g.fill(new Ellipse2D.Double(p.x - r, p.y - r, d, d));
    }

    /**
     *
     * @param bond
     * @param color
     * @param g
     */
    public void drawHighlight(IBond bond, Color color, Graphics2D g) {
        Stroke stroke = g.getStroke();
        g.setStroke(new BasicStroke(params.highlightBondStroke));
        Point2d p0 = bond.getAtom(0).getPoint2d();
        Point2d p1 = bond.getAtom(1).getPoint2d();
        drawLine(p0, p1, g);
        g.setStroke(stroke);
    }

    /**
     *
     * @param highlightContainer
     * @param g
     */
    public void drawHighlightContainer(IAtomContainer highlightContainer, Graphics2D g) {
        if (params.highlightsAbove) {
            drawHighlightContainer(highlightContainer, translucentHighlightColor, g);
        } else {
            drawHighlightContainer(highlightContainer, opaqueHighlightColor, g);
        }
    }

    /**
     *
     * @param highlightContainer
     * @param color
     * @param g
     */
    public void drawHighlightContainer(
            IAtomContainer highlightContainer, Color color, Graphics2D g) {
        Color actualColor;
        if (params.highlightsAbove) {
            actualColor = getTranslucentColor(color);
        } else {
            actualColor = color;
        }
//        System.out.println(color + " " + color.getAlpha() + " " + actualColor + actualColor.getAlpha());
        g.setColor(actualColor);
        double r = params.highlightRadius;
        double d = r * 2;
        for (IAtom atom : highlightContainer.atoms()) {
            Point2d p = atom.getPoint2d();
            g.fill(new Ellipse2D.Double(p.x - r, p.y - r, d, d));
        }

        Stroke stroke = g.getStroke();
        g.setStroke(new BasicStroke(params.highlightBondStroke));
        for (IBond bond : highlightContainer.bonds()) {
            Point2d p0 = bond.getAtom(0).getPoint2d();
            Point2d p1 = bond.getAtom(1).getPoint2d();
            drawLine(p0, p1, g);
        }
        g.setStroke(stroke);
    }

    /**
     *
     */
    @Override
    public void clearHighlights() {
        atomColorMap.clear();
        bondColorMap.clear();
    }

}
