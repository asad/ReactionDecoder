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
package uk.ac.ebi.reactionblast.tools.descriptors;

import java.io.FileInputStream;
import java.io.InputStreamReader;
import static java.lang.Boolean.TRUE;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getHeavyAtoms;
import org.openscience.smsd.helper.MoleculeInitializer;
import uk.ac.ebi.reactionblast.interfaces.IMolDescriptors;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.checkAndCleanMolecule;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CDKMolecularDescriptor extends MoleculeInitializer implements IMolDescriptors {

     private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CDKMolecularDescriptor.class);

    private final IAtomContainer molecule;

    /**
     *
     * @param molecule
     * @throws CDKException
     */
    public CDKMolecularDescriptor(IAtomContainer molecule) throws CDKException {
        super();
        this.molecule = checkAndCleanMolecule(molecule);
        initializeMolecule(molecule);
    }

    /**
     *
     * @param molFile
     * @throws java.lang.Exception
     */
    public CDKMolecularDescriptor(String molFile) throws Exception {
        super();
        FileInputStream ReadMolecule = new FileInputStream(molFile);
        MDLV2000Reader MolRead
                = new MDLV2000Reader(new InputStreamReader(ReadMolecule));
        MolRead.close();
        IAtomContainer newMol = MolRead.read(new AtomContainer());
        this.molecule = checkAndCleanMolecule(newMol);
        initializeMolecule(molecule);
    }

    /**
     *
     * @param checkAromaticity
     * @return
     */
    @Override
    public int getHBondAcceptors(boolean checkAromaticity) {

        IMolecularDescriptor acc = new HBondAcceptorCountDescriptor();
        Object[] hBondparams = {checkAromaticity};
        try {
            acc.setParameters(hBondparams);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        int acceptors = ((IntegerResult) acc.calculate(molecule).getValue()).intValue();
        return acceptors;

    }

    /**
     *
     * @param checkAromaticity
     * @return
     */
    @Override
    public int getHBondDoners(boolean checkAromaticity) {
        IMolecularDescriptor don = new HBondDonorCountDescriptor();
        Object[] hBondparams = {checkAromaticity};
        try {
            don.setParameters(hBondparams);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        int donors = ((IntegerResult) don.calculate(molecule).getValue()).intValue();
        return donors;

    }

    /**
     *
     * @param checkAromaticity
     * @return
     */
    @Override
    public double getXlogP(boolean checkAromaticity) {
        IMolecularDescriptor xlogP = new XLogPDescriptor();
        Object[] xlogPparams = {
            checkAromaticity, TRUE};
        try {
            xlogP.setParameters(xlogPparams);
        } catch (CDKException ex) {
             LOGGER.error(SEVERE, null, ex);
        }
        double xlogPvalue = ((DoubleResult) xlogP.calculate(molecule).getValue()).doubleValue();
        return xlogPvalue;

    }

    /**
     *
     * @param checkAromaticity
     * @return Molecular Polar Surface Area
     */
    @Override
    public double getTPSADescriptor(boolean checkAromaticity) {
        TPSADescriptor tpsa = new TPSADescriptor();
        Object[] tpsaParameter = {checkAromaticity};
        try {
            tpsa.setParameters(tpsaParameter);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        double tpsaValue = ((DoubleResult) tpsa.calculate(molecule).getValue()).doubleValue();
        return tpsaValue;
    }

    /**
     *
     * @return Charged Partial Surface Area (CPSA) descriptors
     */
    @Override
    public double getCPSADescriptor() {

        CPSADescriptor cpsa = new CPSADescriptor();
        double cpsaValue = ((DoubleResult) cpsa.calculate(molecule).getValue()).doubleValue();
        return cpsaValue;

    }

    /**
     *
     * @return Zagreb index: the sum of the squares of atom degree over all
     * heavy atoms i. This can be used instead of heavy atom count
     */
    @Override
    public double getZagrebIndexDescriptor() {

        ZagrebIndexDescriptor zag = new ZagrebIndexDescriptor();
        double value = ((DoubleResult) zag.calculate(molecule).getValue()).doubleValue();
        return value;

    }

    /**
     *
     * @return Sum of the atomic polarizabilities (including implicit
     * hydrogens). This class need explicit hydrogens.
     */
    @Override
    public double getAPolDescriptor() {
        APolDescriptor apol = new APolDescriptor();
        double value = ((DoubleResult) apol.calculate(molecule).getValue()).doubleValue();
        return value;

    }

    /**
     *
     * @return the complexity of a system
     */
    @Override
    public double getFragmentComplexityDescriptor() {

        FragmentComplexityDescriptor fcd = new FragmentComplexityDescriptor();
        double value = ((DoubleResult) fcd.calculate(molecule).getValue()).doubleValue();
        return value;
    }

    /**
     *
     * @return Petitjean Number of a molecule.
     */
    @Override
    public double getPetitjeanNumberDescriptor() {

        PetitjeanNumberDescriptor pnd = new PetitjeanNumberDescriptor();
        double value = ((DoubleResult) pnd.calculate(molecule).getValue()).doubleValue();
        return value;
    }

    /**
     *
     * @return
     */
    @Override
    public double getVAdjMaDescriptor() {
        VAdjMaDescriptor vmd = new VAdjMaDescriptor();
        double value = ((DoubleResult) vmd.calculate(molecule).getValue()).doubleValue();
        return value;
    }

    /**
     *
     * @return Molecular weight of the molecule
     */
    @Override
    public double getMolecularWeight() {

        WeightDescriptor mw = new WeightDescriptor();
        Object[] mwparams = {"*"};
        try {
            mw.setParameters(mwparams);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        double value = ((DoubleResult) mw.calculate(molecule).getValue()).doubleValue();
        return value;
    }

    /**
     *
     * @param includeTerminals
     * @param excludeAmides
     * @return The number of rotatable bonds is given by the SMARTS specified by
     * Daylight on SMARTS tutorial
     */
    @Override
    public int getRotatableBondsCountDescriptor(boolean includeTerminals, boolean excludeAmides) {
        RotatableBondsCountDescriptor rbcd = new RotatableBondsCountDescriptor();
        Object[] params = {includeTerminals, excludeAmides};
        int value = -1;
        try {
            rbcd.setParameters(params);
            value = ((IntegerResult) rbcd.calculate(molecule).getValue()).intValue();
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return value;
    }

    /**
     *
     * @param checkAromaticity
     * @return he number of atoms in the largest pi system.
     */
    @Override
    public int getLargestPiSystemDescriptor(boolean checkAromaticity) {
        LargestPiSystemDescriptor lps = new LargestPiSystemDescriptor();
        Object[] mwparams = {checkAromaticity};
        try {
            lps.setParameters(mwparams);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        int value = ((IntegerResult) lps.calculate(molecule).getValue()).intValue();
        return value;
    }

    /**
     *
     * @return A topological descriptor combining distance and adjacency
     * information.
     */
    @Override
    public int getEccentricConnectivityIndexDescriptor() {

        EccentricConnectivityIndexDescriptor eci = new EccentricConnectivityIndexDescriptor();
        int value = ((IntegerResult) eci.calculate(molecule).getValue()).intValue();
        return value;
    }

    /**
     *
     * @return returns total charge on the molecule
     */
    @Override
    public double getTotalCharge() {

        double count = ExtAtomContainerManipulator.getTotalCharge(molecule);
        return count;

    }

    @Override
    public int getHeavyAtomCount() {

        int count = getHeavyAtoms(molecule).size();
        return count;

    }

    /**
     *
     * @return
     */
    @Override
    public int getCovalentCount() {

        int count = getBondArray(molecule).length;
        return count;

    }
}
