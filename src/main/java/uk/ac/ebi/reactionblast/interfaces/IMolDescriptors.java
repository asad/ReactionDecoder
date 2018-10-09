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
package uk.ac.ebi.reactionblast.interfaces;

import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @contact e-mail: asad@ebi.ac.uk
 */
public interface IMolDescriptors {

    /**
     *
     * @return Sum of the atomic polarizabilities (including implicit
     * hydrogens). This class need explicit hydrogens.
     */
    double getAPolDescriptor();

    /**
     *
     * @return Charged Partial Surface Area (CPSA) descriptors
     * @throws org.openscience.cdk.exception.CDKException
     */
    double getCPSADescriptor() throws CDKException;

    /**
     *
     * @return
     */
    int getCovalentCount();

    /**
     *
     * @return A topological descriptor combining distance and adjacency
     * information.
     */
    int getEccentricConnectivityIndexDescriptor();

    /**
     *
     * @return the complexity of a system
     */
    double getFragmentComplexityDescriptor();

    /**
     *
     * @param checkAromaticity
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    int getHBondAcceptors(boolean checkAromaticity) throws CDKException;

    /**
     *
     * @param checkAromaticity
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    int getHBondDoners(boolean checkAromaticity) throws CDKException;

    /**
     *
     * @return Heavy atom c
     */
    int getHeavyAtomCount();

    /**
     *
     * @param checkAromaticity
     * @return he number of atoms in the largest pi system.
     * @throws org.openscience.cdk.exception.CDKException
     */
    int getLargestPiSystemDescriptor(boolean checkAromaticity) throws CDKException;

    /**
     *
     * @return Molecular weight of the molecule
     * @throws org.openscience.cdk.exception.CDKException
     */
    double getMolecularWeight() throws CDKException;

    /**
     *
     * @return Petitjean Number of a molecule.
     */
    double getPetitjeanNumberDescriptor();

    /**
     *
     * @param includeTerminals
     * @param excludeAmides
     * @return The number of rotatable bonds is given by the SMARTS specified by
     * Daylight on SMARTS tutorial
     * @throws org.openscience.cdk.exception.CDKException
     */
    int getRotatableBondsCountDescriptor(boolean includeTerminals, boolean excludeAmides) throws CDKException;

    /**
     *
     * @param checkAromaticity
     * @return Molecular Polar Surface Area
     * @throws org.openscience.cdk.exception.CDKException
     */
    double getTPSADescriptor(boolean checkAromaticity) throws CDKException;

    /**
     *
     * @return returns total charge on the molecule
     */
    double getTotalCharge();

    /**
     *
     * @return
     */
    double getVAdjMaDescriptor();

    /**
     *
     * @param checkAromaticity
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    double getXlogP(boolean checkAromaticity) throws CDKException;

    /**
     *
     * @return Zagreb index: the sum of the squares of atom degree over all
     * heavy atoms i. This can be used instead of heavy atom count
     */
    double getZagrebIndexDescriptor();
}
