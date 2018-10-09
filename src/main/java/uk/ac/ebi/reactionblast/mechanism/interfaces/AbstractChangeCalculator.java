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
package uk.ac.ebi.reactionblast.mechanism.interfaces;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.helper.MoleculeMoleculePair;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactantProductPair;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionCenterFragment;
import uk.ac.ebi.reactionblast.mechanism.helper.Utility;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCircularFragment;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getSMILES;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class AbstractChangeCalculator extends Utility {

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
     * @return
     */
    public abstract Collection<ReactionCenterFragment> getReactionCenterFragmentList();

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
    /**
     *
     * @return
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

    /**
     *
     * @param reactantAtom
     * @param productAtom
     * @param atomContainerR
     * @param atomContainerP
     * @return
     * @throws Exception
     */
    protected static MoleculeMoleculePair getMolMolPair(
            IAtom reactantAtom,
            IAtom productAtom,
            IAtomContainer atomContainerR,
            IAtomContainer atomContainerP) throws Exception {

        int atomIndexR = getAtomIndexByID(atomContainerR, reactantAtom);

        String signatureR1 = getSignature(atomContainerR, reactantAtom, 1);
        String signatureR2 = getSignature(atomContainerR, reactantAtom, 2);
        String signatureR3 = getSignature(atomContainerR, reactantAtom, 3);
        String signatureR = getSignature(atomContainerR, reactantAtom, -1);

        IAtomContainer fragR1 = getCircularFragment(atomContainerR, atomIndexR, 1);
        IAtomContainer fragR2 = getCircularFragment(atomContainerR, atomIndexR, 2);
        IAtomContainer fragR3 = getCircularFragment(atomContainerR, atomIndexR, 3);
        IAtomContainer fragR = getCircularFragment(atomContainerR, atomIndexR, -1);

        String signatureP1 = getSignature(atomContainerP, productAtom, 1);
        String signatureP2 = getSignature(atomContainerP, productAtom, 2);
        String signatureP3 = getSignature(atomContainerP, productAtom, 3);
        String signatureP = getSignature(atomContainerP, productAtom, -1);

        int atomIndexP = getAtomIndexByID(atomContainerP, productAtom);

        IAtomContainer fragP1 = getCircularFragment(atomContainerP, atomIndexP, 1);
        IAtomContainer fragP2 = getCircularFragment(atomContainerP, atomIndexP, 2);
        IAtomContainer fragP3 = getCircularFragment(atomContainerP, atomIndexP, 3);
        IAtomContainer fragP = getCircularFragment(atomContainerP, atomIndexP, -1);

        IReaction reaction1 = new Reaction();
        reaction1.addReactant(fragR1, 1.0);
        reaction1.addProduct(fragP1, 1.0);

        IReaction reaction2 = new Reaction();
        reaction2.addReactant(fragR2, 1.0);
        reaction2.addProduct(fragP2, 1.0);

        IReaction reaction3 = new Reaction();
        reaction3.addReactant(fragR3, 1.0);
        reaction3.addProduct(fragP3, 1.0);

        IReaction reaction = new Reaction();
        reaction.addReactant(fragR, 1.0);
        reaction.addProduct(fragP, 1.0);

        String smirks;
        String smirks1;
        String smirks2;
        String smirks3;

        smirks = getSMILES(reaction, true);
        smirks1 = getSMILES(reaction1, true);
        smirks2 = getSMILES(reaction2, true);
        smirks3 = getSMILES(reaction3, true);

        String smartsR;
        String smartsR1;
        String smartsR2;
        String smartsR3;

        smartsR = getSMILES(fragR, true);
        smartsR1 = getSMILES(fragR1, true);
        smartsR2 = getSMILES(fragR2, true);
        smartsR3 = getSMILES(fragR3, true);

        String smartsP;
        String smartsP1;
        String smartsP2;
        String smartsP3;

        smartsP = getSMILES(fragP, true);
        smartsP1 = getSMILES(fragP1, true);
        smartsP2 = getSMILES(fragP2, true);
        smartsP3 = getSMILES(fragP3, true);

        ReactantProductPair rrpName = new ReactantProductPair(atomContainerR.getID(), atomContainerP.getID());
        ReactantProductPair rrpSMARTS = new ReactantProductPair(smartsR, smartsP);
        ReactantProductPair rrpSignature = new ReactantProductPair(signatureR, signatureP);

        MoleculeMoleculePair mmp = new MoleculeMoleculePair(rrpName, rrpSMARTS, rrpSignature, smirks);
        mmp.setSignature1(new ReactantProductPair(signatureR1, signatureP1));
        mmp.setSignature2(new ReactantProductPair(signatureR2, signatureP2));
        mmp.setSignature3(new ReactantProductPair(signatureR3, signatureP3));
        mmp.setSmarts1(new ReactantProductPair(smartsR1, smartsP1));
        mmp.setSmarts2(new ReactantProductPair(smartsR2, smartsP2));
        mmp.setSmarts3(new ReactantProductPair(smartsR3, smartsP3));
        mmp.setSmirks1(smirks1);
        mmp.setSmirks2(smirks2);
        mmp.setSmirks3(smirks3);

        return mmp;
    }
}
