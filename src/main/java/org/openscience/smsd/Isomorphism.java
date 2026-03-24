/*
 * Copyright (C) 2009-2020  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd;

import java.io.Serializable;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.Algorithm;
import com.bioinception.smsd.core.SMSD;
import com.bioinception.smsd.core.ChemOptions;

/**
 * Isomorphism adapter that delegates MCS computation to SMSD 3.0.0.
 *
 * Maintains the same public API as the original Isomorphism class
 * but uses the new SMSD 3.0.0 engine for all algorithm work.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public final class Isomorphism extends BaseMapping implements Serializable {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Isomorphism.class);
    static final long serialVersionUID = 0x24845e5c5ae877L;
    private final Algorithm algorithmType;

    /**
     * Initialize query and target molecules and perform MCS search.
     *
     * @param query query mol
     * @param target target mol
     * @param algorithmType algorithm selection (all delegated to SMSD 3.0.0)
     * @param am atom matcher
     * @param bm bond matcher
     * @throws CDKException
     */
    public Isomorphism(
            IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher am, BondMatcher bm) throws CDKException {
        super(query, target, am, bm);
        this.algorithmType = algorithmType;
        mcsBuilder(super.getQuery(), super.getTarget());
        super.setSubgraph(isSubgraph());
    }

    /**
     * Initialize with IQueryAtomContainer.
     *
     * @param query query mol
     * @param target target mol
     * @param algorithmType algorithm selection
     * @throws CDKException
     */
    public Isomorphism(
            IQueryAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType) throws CDKException {
        super(query, target, AtomMatcher.forQuery(), BondMatcher.forQuery());
        this.algorithmType = algorithmType;
        mcsBuilder(super.getQuery(), super.getTarget());
        super.setSubgraph(isSubgraph());
    }

    private void mcsBuilder(IAtomContainer mol1, IAtomContainer mol2) throws CDKException {
        int rAtomCount = mol1.getAtomCount();
        int pAtomCount = mol2.getAtomCount();

        if (rAtomCount == 1 || pAtomCount == 1) {
            singleAtomMapping(mol1, mol2);
        } else {
            smsdMCSAlgorithm(mol1, mol2);
        }
    }

    /**
     * Delegates MCS computation to SMSD 3.0.0.
     */
    private void smsdMCSAlgorithm(IAtomContainer mol1, IAtomContainer mol2) throws CDKException {
        try {
            // First try substructure check
            if (mol1.getAtomCount() > 1 && mol2.getAtomCount() > 1) {
                Substructure sub;
                if (mol1 instanceof IQueryAtomContainer) {
                    sub = new Substructure((IQueryAtomContainer) mol1, mol2, atomMatcher, bondMatcher, true);
                } else {
                    sub = new Substructure(mol1, mol2, atomMatcher, bondMatcher, true);
                }
                if (sub.isSubgraph()) {
                    clearMaps();
                    getMCSList().addAll(sub.getAllAtomMapping());
                    return;
                }
            }

            // Fall back to MCS via SMSD 3.4.0
            ChemOptions chemOptions = new ChemOptions();
            SMSD smsd = new SMSD(mol1, mol2, chemOptions);
            Map<Integer, Integer> mcsResult = smsd.findMCS();

            clearMaps();
            if (mcsResult != null && !mcsResult.isEmpty()) {
                AtomAtomMapping aam = new AtomAtomMapping(mol1, mol2);
                for (Map.Entry<Integer, Integer> entry : mcsResult.entrySet()) {
                    int qIdx = entry.getKey();
                    int tIdx = entry.getValue();
                    if (qIdx >= 0 && qIdx < mol1.getAtomCount()
                            && tIdx >= 0 && tIdx < mol2.getAtomCount()) {
                        IAtom qAtom = mol1.getAtom(qIdx);
                        IAtom tAtom = mol2.getAtom(tIdx);
                        if (qAtom != null && tAtom != null) {
                            aam.put(qAtom, tAtom);
                        }
                    }
                }
                if (!aam.isEmpty()) {
                    getMCSList().add(aam);
                }
            }
        } catch (Exception e) {
            LOGGER.error(Level.SEVERE, "Error in SMSD MCS computation", e);
            throw new CDKException("MCS computation failed: " + e.getMessage(), e);
        }
    }

    /**
     * Handle single atom mapping directly.
     */
    private void singleAtomMapping(IAtomContainer mol1, IAtomContainer mol2) {
        clearMaps();
        if (mol1.getAtomCount() >= 1 && mol2.getAtomCount() >= 1) {
            for (IAtom qAtom : mol1.atoms()) {
                for (IAtom tAtom : mol2.atoms()) {
                    if (AtomMatcher.matchSymbol(qAtom, tAtom)) {
                        AtomAtomMapping aam = new AtomAtomMapping(mol1, mol2);
                        aam.put(qAtom, tAtom);
                        getMCSList().add(aam);
                    }
                }
            }
        }
    }

    /**
     * @return true if query is a subgraph of the target
     */
    @Override
    public boolean isSubgraph() {
        float mappingSize;
        if (getMappingCount() > 0) {
            mappingSize = getAllAtomMapping().iterator().next().getCount();
        } else {
            return false;
        }
        int sourceAtomCount = getQuery().getAtomCount();
        int targetAtomCount = getTarget().getAtomCount();

        if (mappingSize == sourceAtomCount && mappingSize <= targetAtomCount) {
            if (mappingSize == 1) {
                return true;
            } else if (!getAllBondMaps().isEmpty()
                    && getAllBondMaps().iterator().next().size() == getQuery().getBondCount()) {
                return true;
            }
        }
        return false;
    }

}
