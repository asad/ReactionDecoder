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
package uk.ac.ebi.reactionblast.signature;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import static signature.AbstractVertexSignature.parse;
import signature.ColoredTree;
import static uk.ac.ebi.reactionblast.tools.labelling.AtomContainerAtomPermutor.permute;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;

/**
 * Signature implementation specific to rBLAST.
 *
 * @author maclean
 * @author modified by Asad to use rBLAST SMILES
 */
public class RBlastMoleculeSignature extends BaseMoleculeSignature {

    private boolean useAromatics = true;
    private boolean useCharge = true;
    private boolean isBondSensitive = true;
    private RBlastMoleculeFromSignatureBuilder builder;

    /**
     * Make an object that acts as a factory for atom signatures and can also
     * produce molecule signatures.
     *
     * @param atomContainer
     */
    public RBlastMoleculeSignature(IAtomContainer atomContainer) {
        super(atomContainer);
        builder = new RBlastMoleculeFromSignatureBuilder(atomContainer.getBuilder());
    }

    /**
     *
     * @return
     */
    public int getAtomCount() {
        return this.getVertexCount();
    }

    /**
     * @return true if aromatic symbols will be used on aromatic bonds
     */
    public boolean isUseAromatics() {
        return useAromatics;
    }

    /**
     * @param useAromatics if true, will use aromatic symbols
     */
    public void setUseAromatics(boolean useAromatics) {
        this.useAromatics = useAromatics;
    }

    /**
     * @return true if charge symbols will be output
     */
    public boolean isUseCharge() {
        return useCharge;
    }

    /**
     * @param useCharge if true, charge symbols will be output
     */
    public void setUseCharge(boolean useCharge) {
        this.useCharge = useCharge;
    }

    /**
     * @return true if bond order symbols are used
     */
    public boolean isBondSensitive() {
        return isBondSensitive;
    }

    /**
     * @param isBondSensitive if true, bond order symbols will be used
     */
    public void setBondSensitive(boolean isBondSensitive) {
        this.isBondSensitive = isBondSensitive;
    }

    /**
     * Get the canonical signature string for the entire molecule. To do this,
     * signatures are made for each atom, and the lexicographically minimal one
     * is returned.
     *
     * @return
     */
    public String getMoleculeCanonicalSignatureString() {
        return super.toCanonicalString();
    }

    /**
     * Get the canonical signature string for a particular atom. It is canonical
     * in the sense that
     *
     * @param atomIndex
     * @return
     */
    public String getSignatureStringForAtom(int atomIndex) {
        return getAtomSignature(atomIndex).toCanonicalString();
    }

    /**
     *
     * @param atomIndex
     * @param height
     * @return
     */
    public String getSignatureStringForAtom(int atomIndex, int height) {
        return getAtomSignature(atomIndex, height).toCanonicalString();
    }

    /**
     *
     * @param atomIndex
     * @return
     */
    public RBlastAtomSignature getAtomSignature(int atomIndex) {
        RBlastAtomSignature atomSignature
                = (RBlastAtomSignature) signatureForVertex(atomIndex);
        setFlags(atomSignature);
        return atomSignature;
    }

    /**
     *
     * @param atomIndex
     * @param height
     * @return
     */
    public RBlastAtomSignature getAtomSignature(int atomIndex, int height) {
        RBlastAtomSignature atomSignature = new RBlastAtomSignature(atomIndex, atomContainer, height);
        setFlags(atomSignature);
        return atomSignature;
    }

    private void setFlags(RBlastAtomSignature atomSignature) {
        atomSignature.setUseAromatics(useAromatics);
        atomSignature.setUseCharge(useCharge);
        atomSignature.setBondSensitive(isBondSensitive);
    }

    private void setFlags(RBlastMoleculeFromSignatureBuilder builder) {
        builder.setUseAromatics(useAromatics);
        builder.setUseCharge(useCharge);
        builder.setBondSensitive(isBondSensitive);
    }

    /**
     *
     * @param atomIndex
     * @param height
     * @return
     */
    public String getSmilesForAtomSignature(int atomIndex, int height) {
        return getSmilesForAtomSignature(atomIndex, height, new RBlastMoleculeSignatureLabellingAdaptor());
    }

    /**
     * Get a fragment as a smiles with the atoms ordered by canonical signature.
     *
     * @param atomIndex
     * @param height
     * @param labeller
     * @return
     */
    public String getSmilesForAtomSignature(int atomIndex, int height, ICanonicalMoleculeLabeller labeller) {
        String atomSignatureString = getSignatureStringForAtom(atomIndex, height);
        ColoredTree tree = parse(atomSignatureString);

        builder = new RBlastMoleculeFromSignatureBuilder(atomContainer.getBuilder());
        setFlags(builder);
        builder.makeFromColoredTree(tree);
        IAtomContainer fragment = builder.getAtomContainer();
        fragment = permute(labeller.getCanonicalPermutation(fragment), fragment);
//        System.out.println(new AtomContainerPrinter().toString(fragment));
//        System.out.println(Arrays.toString(labeller.getCanonicalPermutation(fragment)));
//        return "";
        RBlastSmilesGenerator smilesGenerator
                = new RBlastSmilesGenerator(false, labeller);
        try {
            return smilesGenerator.createSMILESWithoutCheckForMultipleMolecules(
                    fragment, false, new boolean[fragment.getBondCount()]);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            return "";
        }
    }

    /**
     * Get a fragment of the underlying molecule.
     *
     * @param atomIndex
     * @param height
     * @return
     */
    public IAtomContainer getFragment(int atomIndex, int height) {
        String atomSignatureString = getSignatureStringForAtom(atomIndex, height);
        return makeMoleculeFromSignature(atomSignatureString);
    }

    /**
     * Convert a signature string into a molecule.
     *
     * @param signatureString
     * @return
     */
    public IAtomContainer makeMoleculeFromSignature(String signatureString) {
        ColoredTree tree = parse(signatureString);
        builder = new RBlastMoleculeFromSignatureBuilder(atomContainer.getBuilder());
        setFlags(builder);
        builder.makeFromColoredTree(tree);
        return builder.getAtomContainer();
    }
}
