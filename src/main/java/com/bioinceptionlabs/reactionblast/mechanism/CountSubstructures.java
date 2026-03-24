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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.Serializable;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.openscience.cdk.smiles.smarts.parser.SMARTSParser.parse;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.AtomBondMatcher.AtomMatcher;
import org.openscience.smsd.AtomBondMatcher.BondMatcher;
import org.openscience.smsd.MoleculeInitializer;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
class CountSubstructures extends MoleculeInitializer implements Serializable {

    private static final ILoggingTool LOGGER
            = createLoggingTool(CountSubstructures.class);
    private static final long serialVersionUID = 12343289751445148L;
    private final SmilesParser sp;
    private IAtomContainer mol;

    CountSubstructures(IAtomContainer atomContainer) throws CloneNotSupportedException {
        sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        try {
            this.mol = null;
            mol = removeHydrogensExceptSingleAndPreserveAtomID(atomContainer);
            initializeMolecule(mol);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    public int substructureSize(String smiles) throws CDKException {
        AtomMatcher atomMatcher = AtomBondMatcher.atomMatcher(false, false);
        BondMatcher bondMatcher = AtomBondMatcher.bondMatcher(false, false);

        try {
            IAtomContainer parseSmiles = sp.parseSmiles(smiles);
            Substructure sub = new Substructure(parseSmiles, mol, atomMatcher, bondMatcher, false);
            return sub.isSubgraph() ? sub.getFirstAtomMapping().getCount() : 0;
        } catch (InvalidSmilesException ex) {
            Substructure sub = new Substructure(parse(smiles, mol.getBuilder()), mol, false);
            return sub.isSubgraph() ? sub.getFirstAtomMapping().getCount() : 0;
        }
    }
}
