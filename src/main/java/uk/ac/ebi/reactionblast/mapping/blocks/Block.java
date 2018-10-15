/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.blocks;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point2d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author asad
 */
public class Block implements Comparable<Block> {

    private final IAtomContainer atomContainer;

//    private List<IAtom> atoms;
    private final Map<IAtom, IAtom> atomMap;

    private SubgraphMoleculeSignature subgraphSignature;

    private String signatureString;

    private Point2d centerPoint;

    private Rectangle2D bounds;

    private Block partner;

    /**
     *
     * @param atomContainer
     */
    public Block(IAtomContainer atomContainer) {
        this.atomContainer = atomContainer;
        this.atomMap = new HashMap<>();
    }

    /**
     *
     * @return
     */
    public int[] getMappingPermutation() {
        int[] mappingPermutation = new int[atomMap.size()];

        Map<Integer, Integer> indexMap = new HashMap<>();
        atomMap.keySet().stream().forEach((atom) -> {
            int atomIndex = atomContainer.indexOf(atom);
            IAtom partnerAtom = atomMap.get(atom);
            int partnerIndex = partner.atomContainer.indexOf(partnerAtom);
            indexMap.put(atomIndex, partnerIndex);
        });//        System.out.println("indexMap " + indexMap);

        Map<Integer, Integer> compactMap
                = getCompactMap(new ArrayList<>(indexMap.keySet()));
        Map<Integer, Integer> compactPartnerMap
                = getCompactMap(new ArrayList<>(indexMap.values()));

        indexMap.keySet().stream().forEach((index) -> {
            mappingPermutation[compactMap.get(index)]
                    = compactPartnerMap.get(indexMap.get(index));
        });

        return mappingPermutation;
    }

    private Map<Integer, Integer> getCompactMap(List<Integer> indices) {
        Map<Integer, Integer> compactMap = new HashMap<>();
        sort(indices);
        int compactIndex = 0;
        for (Integer index : indices) {
            compactMap.put(index, compactIndex);
            compactIndex++;
        }
        return compactMap;
    }

    /**
     *
     * @param partner
     */
    public void setPartner(Block partner) {
        this.partner = partner;
    }

    /**
     *
     * @return
     */
    public Block getPartner() {
        return this.partner;
    }

    /**
     *
     * @param atom
     * @param partner
     */
    public void addMapping(IAtom atom, IAtom partner) {
        atomMap.put(atom, partner);
    }

    /**
     *
     * @param centerPoint
     */
    public void setCenterPoint(Point2d centerPoint) {
        this.centerPoint = centerPoint;
    }

    /**
     *
     * @return
     */
    public Point2d getCenterPoint() {
        if (centerPoint == null) {
            Rectangle2D b = getBounds();
            centerPoint = new Point2d(b.getCenterX(), b.getCenterY());
        }
        return centerPoint;
    }

    /**
     *
     * @return
     */
    public Rectangle2D getBounds() {
        if (bounds == null) {
            IAtomContainer tmp
                    = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            tmp.setAtoms(getAtoms().toArray(new IAtom[0]));
            bounds = calculateBounds(tmp);
        }
        return bounds;
    }

    private Rectangle2D calculateBounds(IAtomContainer atomContainer) {
        Rectangle2D b = null;
        for (IAtom atom : atomContainer.atoms()) {
            Point2d p = atom.getPoint2d();
            if (b == null) {
                b = new Rectangle2D.Double(p.x, p.y, 0, 0);
            } else {
                b.add(p.x, p.y);
            }
        }
        return b;
    }

    @Override
    public int compareTo(Block o) {
        return getSignatureString().compareTo(o.getSignatureString());
    }

    /**
     *
     * @return
     */
    public int getAtomCount() {
        return atomMap.size();
    }

    /**
     *
     * @return
     */
    public int[] getLabels() {
        return getSubgraphSignature().getCanonicalLabels();
    }

    /**
     *
     * @return
     */
    public SubgraphMoleculeSignature getSubgraphSignature() {
        if (subgraphSignature == null) {
            subgraphSignature
                    = new SubgraphMoleculeSignature(atomContainer, getAtoms(), -1);
        }
        return subgraphSignature;
    }

    /**
     *
     * @return
     */
    public String getSignatureString() {
        if (signatureString == null) {
            signatureString = getSubgraphSignature().toCanonicalString();
        }
        return signatureString;
    }

    /**
     *
     * @return
     */
    public IAtomContainer getAtomContainer() {
        return atomContainer;
    }

    /**
     *
     * @return
     */
    public List<IAtom> getAtoms() {
        return new ArrayList<>(atomMap.keySet());
    }

    private List<Integer> getIndices(List<IAtom> atoms, IAtomContainer container) {
        List<Integer> indices = new ArrayList<>();
        atoms.stream().forEach((atom) -> {
            indices.add(container.indexOf(atom));
        });
        return indices;
    }

    @Override
    public String toString() {
        return atomContainer.getID() + " " + getIndices(getAtoms(), atomContainer) + " " + getSignatureString();
    }

}
