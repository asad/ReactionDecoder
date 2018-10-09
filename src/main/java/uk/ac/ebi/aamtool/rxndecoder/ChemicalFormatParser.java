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
package uk.ac.ebi.aamtool.rxndecoder;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.System.err;
import static java.lang.System.exit;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.AtomContainer;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.io.CMLReader;
import static org.openscience.cdk.io.IChemObjectReader.Mode.RELAXED;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class ChemicalFormatParser {

    final static String NEW_LINE = getProperty("line.separator");
    private static final Logger LOG = getLogger(ChemicalFormatParser.class.getName());

    protected IReaction parseCML(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        if (!f.isFile()) {
            err.println("CML file not found! " + f.getName());
            exit(1);
        }
        String[] split = f.getName().split(".cml");
        CMLReader cmlReader = new CMLReader(new FileInputStream(input));
        AtomContainer ac = cmlReader.read(new AtomContainer());
        IReaction r = new Reaction();
        r.addReactant(ac, 1.0);
        r.addProduct(ac, 1.0);
        r.setID(split[0]);
        return r;
    }

    protected IReaction parseRXN(String fileName) {
        File filepath = new File(fileName);
        if (!filepath.isFile()) {
            err.println("RXN file not found! " + filepath.getName());
            exit(1);
        }
        String[] split = filepath.getName().split(".rxn");
        try {
            out.println(NEW_LINE + "Annotating Reaction " + split[0] + NEW_LINE);
            IReaction rxnReactions;
            try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath));) {
                try {
                    rxnReactions = reader.read(new Reaction());
                    reader.close();
                    rxnReactions.setID(split[0]);
                    return rxnReactions;
                } catch (IOException | CDKException ex) {
                    err.println("ERROR in Reading Reaction file " + filepath + NEW_LINE + ex);
                }
            }
        } catch (IOException ex) {
            err.println("Failed to Read and Annotate RXN File ");
            getLogger(ChemicalFormatParser.class.getName()).log(SEVERE, null, ex);
        }
        return null;
    }

    protected IReaction parseReactionSMILES(String reactionSmiles) {
        SmilesParser sp = new SmilesParser(getInstance());
        try {
            IReaction parseReactionSmiles = sp.parseReactionSmiles(reactionSmiles);
            try {
                out.println(NEW_LINE + "Annotating Reaction " + "smiles" + NEW_LINE);
                parseReactionSmiles.setID("smiles");
                return parseReactionSmiles;
            } catch (Exception ex) {
                getLogger(ChemicalFormatParser.class.getName()).log(SEVERE, null, ex);
            }
        } catch (InvalidSmilesException ex) {
            getLogger(ChemicalFormatParser.class.getName()).log(SEVERE, null, ex);
        }
        return null;
    }

    protected IReaction parseSMILES(String smiles) {
        SmilesParser sp = new SmilesParser(getInstance());
        try {
            IAtomContainer mol = sp.parseSmiles(smiles);
            try {
                IReaction parseReactionSmiles = getInstance().newInstance(IReaction.class);
                parseReactionSmiles.addReactant(mol, 1.0);
//                parseReactionSmiles.addProduct(mol, 1.0);
                out.println(NEW_LINE + "Annotating Reaction " + "smiles" + NEW_LINE);
                parseReactionSmiles.setID("smiles");
                return parseReactionSmiles;
            } catch (IllegalArgumentException ex) {
                getLogger(ChemicalFormatParser.class.getName()).log(SEVERE, null, ex);
            }
        } catch (InvalidSmilesException ex) {
            getLogger(ChemicalFormatParser.class.getName()).log(SEVERE, null, ex);
        }
        return null;
    }

    protected IReaction parseMOL2(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        if (!f.isFile()) {
            err.println("Mol2 file not found! " + f.getName());
            exit(1);
        }

        String[] split = f.getName().split(".mol");
        MDLV2000Reader mdlV2000Reader = new MDLV2000Reader(
                new FileReader(input), RELAXED);
        AtomContainer ac = mdlV2000Reader.read(new AtomContainer());
        IReaction r = new Reaction();
        r.addReactant(ac, 1.0);
        r.addProduct(ac, 1.0);
        r.setID(split[0]);
        return r;
    }

    protected IReaction parseSDF(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        if (!f.isFile()) {
            err.println("SDF file not found! " + f.getName());
            exit(1);
        }
        String[] split = f.getName().split(".sdf");
        Mol2Reader mol2Reader = new Mol2Reader(new FileReader(input));
        AtomContainer ac = mol2Reader.read(new AtomContainer());
        IReaction r = new Reaction();
        r.addReactant(ac, 1.0);
        r.addProduct(ac, 1.0);
        r.setID(split[0]);
        return r;
    }
}
