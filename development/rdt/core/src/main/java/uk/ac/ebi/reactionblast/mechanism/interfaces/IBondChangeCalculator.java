/*
 * Copyright (C) 2016 asad.
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
package uk.ac.ebi.reactionblast.mechanism.interfaces;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.helper.MoleculeMoleculePair;

/**
 *
 * @author asad
 */
public interface IBondChangeCalculator {

    /**
     *
     * @return
     */
    Map<IAtom, IAtom> getAtomAtomMappings();

    Map<IBond, String> getBondCleavedReactant();

    Map<IBond, String> getBondFormedProduct();

    Map<IBond, String> getBondOrderProduct();

    Map<IBond, String> getBondOrderReactant();

    /**
     * @return the energyDelta
     */
    int getEnergyDelta();

    /**
     *
     * @return @throws CDKException
     */
    IPatternFingerprinter getFormedCleavedWFingerprint() throws CDKException;

    /**
     *
     * @return
     */
    Map<String, Collection<String>> getMoleculeMoleculeTransformationPairs();

    /**
     *
     * @return @throws CDKException
     */
    IPatternFingerprinter getOrderChangesWFingerprint() throws CDKException;

    /**
     *
     * @return (removed the unchanged H atoms)
     * @throws Exception
     */
    IReaction getReaction() throws Exception;

    /**
     *
     * @return
     */
    Collection<IAtom> getReactionCenterSet();

    /**
     * Return Reaction center Fingerprint
     *
     * @return
     * @throws CDKException
     */
    IPatternFingerprinter getReactionCenterWFingerprint() throws CDKException;

    Collection<MoleculeMoleculePair> getReactionCentreTransformationPairs();

    /**
     *
     * @return (removed the unchanged H atoms)
     */
    IReaction getReactionWithCompressUnChangedHydrogens();

    Map<IAtom, String> getStereoCenterAtomsProduct();

    Map<IAtom, String> getStereoCenterAtomsReactant();

    IPatternFingerprinter getStereoChangesWFingerprint() throws CDKException;

    /**
     *
     * @return @throws CDKException
     */
    int getTotalBondBreakingEnergy() throws CDKException;

    /**
     * @return the totalSmallestFragmentSize
     */
    int getTotalSmallestFragmentSize();

    /**
     * Print the atom changes
     *
     * @return
     */
    @Override
    String toString();

    /**
     * Print the atom changes
     *
     * @param bondChangeInfoFile
     * @throws IOException
     */
    void writeBondChanges(File bondChangeInfoFile) throws IOException;

}
