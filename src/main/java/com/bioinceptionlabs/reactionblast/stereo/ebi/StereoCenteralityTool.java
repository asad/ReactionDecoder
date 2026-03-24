/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.stereo.ebi;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import com.bioinceptionlabs.reactionblast.mechanism.helper.Utility;
import com.bioinceptionlabs.reactionblast.stereo.IStereoAndConformation;

/**
 * Tool for detecting stereochemistry using CDK's built-in stereo perception.
 * Replaces the old centres-based CIP priority analysis with CDK 2.12's
 * native IStereoElement support (ITetrahedralChirality, IDoubleBondStereochemistry).
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public abstract class StereoCenteralityTool extends Utility {

    private static final ILoggingTool LOGGER
            = createLoggingTool(StereoCenteralityTool.class);
    private static final long serialVersionUID = 17867606807697859L;

    private static IAtom getAtomByID(String id, IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            if (a.getID() != null && a.getID().equals(id)) {
                return a;
            }
        }
        return null;
    }

    /**
     * Detect chirality for all molecules in a reaction using CDK 2D stereo.
     * Uses CDK's built-in IStereoElement perception.
     *
     * @param reaction the reaction to analyze
     * @return map of atoms to their stereo configuration
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static Map<IAtom, IStereoAndConformation> getChirality2D(IReaction reaction)
            throws CDKException, CloneNotSupportedException {
        Map<IAtom, IStereoAndConformation> chiralityMap = new HashMap<>();

        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
            Map<IAtom, IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH);
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
            Map<IAtom, IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH);
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        return chiralityMap;
    }

    /**
     * Detect chirality for a single molecule using CDK's IStereoElement.
     * Maps ITetrahedralChirality → R/S and IDoubleBondStereochemistry → E/Z.
     *
     * @param ac molecule to analyze (should have 2D coordinates)
     * @return map of atoms to their stereo configuration
     */
    @SuppressWarnings("unchecked")
    public static Map<IAtom, IStereoAndConformation> getChirality2D(IAtomContainer ac) {
        Map<IAtom, IStereoAndConformation> chiralityMap = new HashMap<>();

        // Default all atoms to NONE
        for (IAtom atom : ac.atoms()) {
            chiralityMap.put(atom, IStereoAndConformation.NONE);
        }

        try {
            // Extract stereo from CDK's stereo elements (set during SMILES parsing or 2D perception)
            for (IStereoElement<?, ?> element : ac.stereoElements()) {
                if (element instanceof ITetrahedralChirality) {
                    ITetrahedralChirality tc = (ITetrahedralChirality) element;
                    IAtom focus = tc.getChiralAtom();
                    ITetrahedralChirality.Stereo stereo = tc.getStereo();
                    if (stereo == ITetrahedralChirality.Stereo.CLOCKWISE) {
                        chiralityMap.put(focus, IStereoAndConformation.R);
                    } else if (stereo == ITetrahedralChirality.Stereo.ANTI_CLOCKWISE) {
                        chiralityMap.put(focus, IStereoAndConformation.S);
                    }
                    if (focus != null) {
                        focus.setProperty("Stereo", chiralityMap.get(focus));
                    }
                } else if (element instanceof IDoubleBondStereochemistry) {
                    IDoubleBondStereochemistry dbs = (IDoubleBondStereochemistry) element;
                    IDoubleBondStereochemistry.Conformation conf = dbs.getStereo();
                    // Mark both atoms of the double bond
                    IAtom a0 = dbs.getStereoBond().getBegin();
                    IAtom a1 = dbs.getStereoBond().getEnd();
                    IStereoAndConformation sc;
                    if (conf == IDoubleBondStereochemistry.Conformation.OPPOSITE) {
                        sc = IStereoAndConformation.E;
                    } else if (conf == IDoubleBondStereochemistry.Conformation.TOGETHER) {
                        sc = IStereoAndConformation.Z;
                    } else {
                        continue;
                    }
                    chiralityMap.put(a0, sc);
                    chiralityMap.put(a1, sc);
                    if (a0 != null) a0.setProperty("Stereo", sc);
                    if (a1 != null) a1.setProperty("Stereo", sc);
                }
            }
        } catch (Exception e) {
            LOGGER.debug("Stereo perception error: " + e.getMessage());
        }

        // Also set as atom properties for downstream consumers
        chiralityMap.forEach((atom, sc) -> {
            if (atom != null) atom.setProperty("Stereo", sc);
        });

        return chiralityMap;
    }
}
