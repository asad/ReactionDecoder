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
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author asad
 */
public interface Highlighter {

    /**
     *
     * @param highlightContainer
     * @param color
     */
    public void addHighlights(IAtomContainer highlightContainer, Color color);

    /**
     *
     * @param atoms
     * @param bonds
     */
    public void addHighlights(List<IAtom> atoms, List<IBond> bonds);

    /**
     *
     * @param molecule
     * @param g
     */
    public void drawHighlights(IAtomContainer molecule, Graphics2D g);

    /**
     *
     * @param colorMap
     */
    public void addToHighlights(Map<IAtom, Color> colorMap);

    /**
     *
     */
    public void clearHighlights();

}
