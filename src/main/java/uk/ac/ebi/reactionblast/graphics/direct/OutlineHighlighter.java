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
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import static org.openscience.cdk.geometry.GeometryTools.get2DCenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author asad
 */
public class OutlineHighlighter extends AbstractHighlightDrawer implements Highlighter {

    private final Map<IAtomContainer, Color> colorMap;

    /**
     *
     * @param params
     */
    public OutlineHighlighter(Params params) {
        super(params);
        colorMap = new HashMap<>();
    }

    /**
     *
     * @param highlightContainer
     * @param color
     */
    @Override
    public void addHighlights(IAtomContainer highlightContainer, Color color) {
        colorMap.put(highlightContainer, color);
    }

    /**
     *
     * @param atoms
     * @param bonds
     */
    @Override
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds) {
        IAtomContainer highlightContainer = null;
        if (atoms.size() > 0) {
            highlightContainer
                    = atoms.get(0).getBuilder().newInstance(IAtomContainer.class);
        } else if (bonds.size() > 0) {
            highlightContainer
                    = bonds.get(0).getBuilder().newInstance(IAtomContainer.class);
        } else {
            return;
        }

        for (IAtom atom : atoms) {
            highlightContainer.addAtom(atom);
        }

        for (IBond bond : bonds) {
            highlightContainer.addBond(bond);
        }
        addHighlights(highlightContainer, params.highlightColor);
    }

    /**
     *
     * @param colorMap
     */
    @Override
    public void addToHighlights(Map<IAtom, Color> colorMap) {
        // TODO Auto-generated method stub
        // ? Problem is that this highlighter intends to outline all
        // the atoms in one atom container with a single color, not
        // color each atom separately
    }

    /**
     *
     * @param molecule
     * @param g
     */
    @Override
    public void drawHighlights(IAtomContainer molecule, Graphics2D g) {
        Point2d center = null;
        List<IAtomContainer> highlightContainers;
        if (params.circularHighlightIsConcentric) {
            highlightContainers = new ArrayList<>(colorMap.keySet());
            sort(highlightContainers, new Comparator<IAtomContainer>() {

                @Override
                public int compare(IAtomContainer ac0, IAtomContainer ac1) {
                    if (ac0.getAtomCount() < ac1.getAtomCount()) {
                        return 1;
                    } else if (ac0.getAtomCount() > ac1.getAtomCount()) {
                        return -1;
                    } else {
                        return 0;
                    }
                }

            });
            center = get2DCenter(
                    highlightContainers.get(highlightContainers.size() - 1));
        } else {
            highlightContainers = new ArrayList<>(colorMap.keySet());
        }

        for (int containerIndex = 0; containerIndex < highlightContainers.size(); containerIndex++) {
            IAtomContainer highlightContainer = highlightContainers.get(containerIndex);
            Color savedColor = g.getColor();
            if (params.circularHighlightTransparentFilled) {
                g.setColor(getTranslucentColor(colorMap.get(highlightContainer)));
            } else {
                g.setColor(colorMap.get(highlightContainer));
            }

            if (!params.circularHighlightIsConcentric || center == null) {
                center = get2DCenter(highlightContainer);
            }
            double maxDist = 0.0;
            for (IAtom highlightAtom : highlightContainer.atoms()) {
                if (molecule.contains(highlightAtom)) {
                    Point2d point = highlightAtom.getPoint2d();
                    if (point != null) {
                        double d = center.distance(point);
                        if (d > maxDist) {
                            maxDist = d;
                        }
                    }
                    if (params.circularHighlightShowAtoms) {
                        double r = params.highlightRadius;
                        g.fill(new Ellipse2D.Double(
                                point.x - r, point.y - r, r * 2, r * 2));
                    }
                }
            }

            double x;
            double y;
            double dim;
            if (highlightContainer.getAtomCount() == 1 && containerIndex == highlightContainers.size() - 1) {
                x = center.x - params.circularHighlightMinRadius;
                y = center.y - params.circularHighlightMinRadius;
                dim = 2 * params.circularHighlightMinRadius;
            } else {
                x = center.x - maxDist;
                y = center.y - maxDist;
                dim = 2 * maxDist;
            }

            if (params.circularHighlightTransparentFilled) {
                g.fill(new Ellipse2D.Double(x, y, dim, dim));
            } else {
                g.draw(new Ellipse2D.Double(x, y, dim, dim));
            }
            g.setColor(savedColor);
        }
    }

    /**
     *
     */
    @Override
    public void clearHighlights() {
        colorMap.clear();
    }

}
