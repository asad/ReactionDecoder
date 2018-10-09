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
package uk.ac.ebi.reactionblast.graphics.direct.layout;

import java.awt.geom.Rectangle2D;
import static java.lang.Boolean.TRUE;
import javax.vecmath.Vector2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import static uk.ac.ebi.reactionblast.graphics.direct.layout.MoleculeAligner.alignToMaxWidth;
import static uk.ac.ebi.reactionblast.graphics.direct.layout.MoleculeAligner.alignToMinAreaBox;

/**
 *
 * @author asad
 * @param <T>
 */
public abstract class AbstractDirectLayout<T> {
    
    /**
     *
     */
    public static final String INVERTED = "Coordinates Inverted";
    
    /**
     *
     */
    protected Params params;
    
    /**
     *
     */
    protected BoundsTree boundsTree;
    
    /**
     *
     */
    public boolean shouldInvert;
    
    /**
     *
     */
    public AbstractDirectLayout() {
        this(true);
    }
    
    /**
     *
     * @param shouldInvert
     */
    public AbstractDirectLayout(boolean shouldInvert) {
        this.shouldInvert = shouldInvert;
    }
    
    /**
     *
     * @param obj
     * @param axis
     * @return
     */
    public abstract BoundsTree layout(T obj, Vector2d axis);

    /**
     *
     * @return
     */
    public abstract Vector2d getAxis();
    
    /**
     *
     * @return
     */
    public abstract double getAxisPosition();
    
    /**
     *
     * @return
     */
    public Params getParams() {
        return params;
    }

    /**
     *
     * @param params
     */
    public void setParams(Params params) {
        this.params = params;
    }
    
    /**
     *
     * @param ac
     * @param x
     * @param y
     * @param bounds
     */
    public void translateTo(IAtomContainer ac, double x, double y, Rectangle2D bounds) {
        double dx = x - bounds.getCenterX();
        double dy = y - bounds.getCenterY();
        for (IAtom atom : ac.atoms()) {
            atom.getPoint2d().x += dx;
            atom.getPoint2d().y += dy;
        }
//        bounds.setFrameFromCenter(x, y, bounds.getMinX() + dx, bounds.getMinY() + dy);
//        System.out.print(ac.getID() + " ADL Before : " + BoundsPrinter.toString(bounds));
        bounds.setRect(bounds.getMinX() + dx, bounds.getMinY() + dy, bounds.getWidth(), bounds.getHeight());
//        System.out.println(" After: " + BoundsPrinter.toString(bounds) + " " + dx + " " + dy);
    }
    
    /**
     *
     * @param ac
     */
    public void invert(IAtomContainer ac) {
        if (shouldInvert && 
                ac.getProperty(INVERTED) == null ||
                !((Boolean)ac.getProperty(INVERTED))) {
            for (IAtom atom : ac.atoms()) {
                atom.getPoint2d().y *= - 1;
            }
            ac.setProperty(INVERTED, TRUE);
        }
        shouldInvert = false;
    }
    
    /**
     *
     * @param atomContainer
     * @param molAxis
     */
    public void align(IAtomContainer atomContainer, Vector2d molAxis) {
        switch (params.moleculeAlignMethod) {
            case MAX_AXIS: alignToMaxWidth(atomContainer, molAxis);
            case MIN_AREA: alignToMinAreaBox(atomContainer, molAxis);
            default: alignToMaxWidth(atomContainer, molAxis);
        }
    }

}
