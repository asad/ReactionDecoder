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

import java.awt.Graphics2D;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.graphics.direct.LabelManager;
import uk.ac.ebi.reactionblast.graphics.direct.Params;
import uk.ac.ebi.reactionblast.graphics.direct.layout.BoundsTree;

/**
 * 'Layout' a molecule in the context of a Graphics object; NOTE - does not
 * re-position the atoms...
 *
 * @author maclean
 *
 */
public class MoleculeLayout extends AbstractAWTLayout<IAtomContainer> {

    private AtomLayout atomLayout;

    /**
     *
     * @param params
     */
    public MoleculeLayout(Params params) {
        atomLayout = new AtomLayout(this, params, new LabelManager());
    }

    /**
     *
     * @param parent
     * @param params
     */
    public MoleculeLayout(AbstractAWTLayout parent, Params params) {
        this(params);
        super.parent = parent;
    }

    /**
     *
     * @param atomContainer
     * @param graphics
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainer atomContainer, Graphics2D graphics) {
        return layout(atomContainer, atomContainer.getID(), graphics);
    }

    /**
     *
     * @param atomContainer
     * @param rootLabel
     * @param graphics
     * @return
     */
    @Override
    public BoundsTree layout(IAtomContainer atomContainer, String rootLabel, Graphics2D graphics) {
        atomLayout.reset();
        setGraphics(graphics);
        currentObject = atomContainer;

        boundsTree = new BoundsTree(rootLabel);
        for (IAtom atom : atomContainer.atoms()) {
            // add all the atom bounds to the tree, with prefix of molID
            boundsTree.add(rootLabel, atomLayout.layout(atom, graphics));
        }

        return boundsTree;
    }

}
