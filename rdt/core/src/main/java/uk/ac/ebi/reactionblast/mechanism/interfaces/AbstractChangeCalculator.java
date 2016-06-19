/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class AbstractChangeCalculator extends ReactionMappingUtility {

    private static final long serialVersionUID = 1523522323512L;

    /*
     * Return KEGG like RPAIRS
     */
    /**
     *
     * @return
     */
    public abstract Map<String, Collection<String>> getMoleculeMoleculeTransformationPairs();

    /**
     * @return bond(s) cleaved at reactant side and with reactant ID
     */
    public abstract Map<IBond, String> getBondCleavedReactant();

    /**
     * @return bond(s) Formed at reactant side and with reactant ID
     */
    public abstract Map<IBond, String> getBondFormedProduct();

    /**
     * @return bond(s) order changed at product side and with product ID
     */
    public abstract Map<IBond, String> getBondOrderProduct();

    /**
     * @return bond(s) order changed at reactant side and with reactant ID
     */
    public abstract Map<IBond, String> getBondOrderReactant();

    /**
     * @return the bond effect by Stereo changes at product side and with
     * product ID
     */
    public abstract Map<IAtom, String> getStereoCenterAtomsProduct();

    /**
     * @return the bond effect by Stereo changes at reactant side and with
     * reactant ID
     */
    public abstract Map<IAtom, String> getStereoCenterAtomsReactant();

    /**
     *
     * @return @throws CDKException
     */
    public abstract IPatternFingerprinter getFormedCleavedWFingerprint() throws CDKException;

    /**
     *
     * @return @throws CDKException
     */
    public abstract IPatternFingerprinter getOrderChangesWFingerprint() throws CDKException;

    /**
     *
     * @return @throws CDKException
     */
    public abstract IPatternFingerprinter getReactionCenterWFingerprint() throws CDKException;

    /**
     *
     * @return (removed the unchanged H atoms)
     * @throws Exception
     */
    public abstract IReaction getReaction() throws Exception;

    /**
     *
     * @return (removed the unchanged H atoms)
     * @throws Exception
     */
    public abstract IReaction getReactionWithCompressUnChangedHydrogens() throws Exception;

    /**
     *
     * @return @throws CDKException
     */
    public abstract IPatternFingerprinter getStereoChangesWFingerprint() throws CDKException;

    /*
     * 
     * @return atom-atom mapping
     */
    public abstract Map<IAtom, IAtom> getAtomAtomMappings();

    /**
     * Print the bond changes
     *
     * @param bondChangeInfoFile
     * @throws IOException
     */
    public abstract void writeBondChanges(File bondChangeInfoFile) throws IOException;

    /**
     * Return the list of Molecule Molecule Transformation List
     *
     * @return
     */
    public abstract Collection<MoleculeMoleculePair> getReactionCentreTransformationPairs();

}
