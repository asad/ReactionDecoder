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
package com.bioinceptionlabs.reactionblast.tools;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.ChemModel;
import static org.openscience.cdk.geometry.GeometryUtil.has2DCoordinates;
import static org.openscience.cdk.graph.ConnectivityChecker.partitionIntoMolecules;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.MoleculeSetManipulator.getAllAtomContainers;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class LayoutCheck {

    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(LayoutCheck.class);

    /**
     *
     * @param mol
     * @return
     */
    public static IAtomContainer getMoleculeWithLayoutCheck(IAtomContainer mol) {
        if (!has2DCoordinates(mol)) {
            try {
                StructureDiagramGenerator sdg = new StructureDiagramGenerator(new AtomContainer(mol));
                sdg.generateCoordinates();
                mol = sdg.getMolecule();
            } catch (Exception e) {
                LOGGER.error(e);
            }
        }

        return mol;
    }

    /**
     *
     * @param mol
     * @return
     */
    public static IChemModel getChemModelWithMoleculeWithLayoutCheck(IAtomContainer mol) {
        IChemModel chemModel = new ChemModel();
        chemModel.setMoleculeSet(partitionIntoMolecules(mol));
        for (IAtomContainer molecule : getAllAtomContainers(chemModel.getMoleculeSet())) {
            if (has2DCoordinates(molecule)) {
                try {
                    StructureDiagramGenerator sdg = new StructureDiagramGenerator(new AtomContainer(molecule));
                    sdg.generateCoordinates();
                    chemModel = (IChemModel) sdg.getMolecule();
                } catch (Exception e) {
                    LOGGER.error(e);
                }
            }
        }
        return chemModel;
    }

    private LayoutCheck() {
    }
}
