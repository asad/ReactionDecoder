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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.Serializable;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;
import org.openscience.smsd.helper.MoleculeInitializer;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
/**
 * Return Substructure match
 */
class CountSubstructures extends MoleculeInitializer implements Serializable {

    private static final long serialVersionUID = 12343289751445148L;
    private SmilesParser sp;
    private IAtomContainer mol;

    public CountSubstructures(IAtomContainer atomContainer) throws CloneNotSupportedException {
        sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        try {
            this.mol = null;
            mol = ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID(atomContainer);
            initializeMolecule(mol);
        } catch (CDKException ex) {
            Logger.getLogger(CountSubstructures.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public int substructureSize(String smiles) throws CDKException {
        try {
            IAtomContainer parseSmiles = sp.parseSmiles(smiles);
            VF2 vf = new VF2(parseSmiles, mol, true, true, true);
            return vf.isSubgraph() ? vf.getFirstAtomMapping().getCount() : 0;
        } catch (InvalidSmilesException ex) {
            VF2 vf = new VF2(SMARTSParser.parse(smiles, mol.getBuilder()), mol);
            return vf.isSubgraph() ? vf.getFirstAtomMapping().getCount() : 0;
        }
    }
}
