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
package aamtool;

import static aamtool.Annotator.NEW_LINE;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class ChemicalFormatParser {
    private static final Logger LOG = Logger.getLogger(ChemicalFormatParser.class.getName());

    protected IReaction parseCML(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        if (!f.isFile()) {
            Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.WARNING, String.format("CML file not found! " + f.getName()));
            System.exit(1);
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

    protected List<IReaction> parseRXN(String fileNames) {
        /*
         split of file extension
         */
        String[] files = fileNames.split(";");
        List<IReaction> reactions = new ArrayList<>();
        for (String file : files) {
            String[] f = file.split("\\.(?=[^\\.]+$)");
            if (f[0].equals("rxn")) {
                continue;
            }
            String fileName = f[0].trim() + ".rxn";
            File filepath = new File(fileName);
            if (!filepath.isFile()) {
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.WARNING, String.format("RXN file not found! %s", filepath.getName()));
                System.exit(1);
            }
            try {
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.INFO, "Annotating Reaction {0}", filepath.getName());
                IReaction rxnReactions;
                try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath));) {
                    try {
                        rxnReactions = reader.read(new Reaction());
                        reader.close();
                        rxnReactions.setID(filepath.getName().split(".rxn")[0]);
                        reactions.add(rxnReactions);
                    } catch (IOException | CDKException ex) {
                        System.err.println("ERROR in Reading Reaction file " + filepath + NEW_LINE + ex);
                    }
                }
            } catch (IOException ex) {
                System.err.println("Failed to Read and Annotate RXN File ");
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return reactions;
    }

    protected List<IReaction> parseReactionSMILES(String reactionSmiles) {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        String[] smiles = reactionSmiles.split("\\s+");
        List<IReaction> reactions = new ArrayList<>();
        int smilesIndex = 1;
        for (String s : smiles) {
            try {
                IReaction parseReactionSmile = sp.parseReactionSmiles(s);
                try {
                    Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.INFO, "Annotating Reaction " + "smiles");
                    if (smiles.length > 1) {
                        parseReactionSmile.setID("smiles_" + smilesIndex);
                    } else {
                        parseReactionSmile.setID("smiles");
                    }
                    reactions.add(parseReactionSmile);
                } catch (Exception ex) {
                    Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.SEVERE, null, ex);
                }
            } catch (InvalidSmilesException ex) {
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.SEVERE, null, ex);
            }
            smilesIndex++;
        }
        return reactions;
    }

    protected IReaction parseSMILES(String smiles) {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        try {
            IAtomContainer mol = sp.parseSmiles(smiles);
            try {
                IReaction parseReactionSmiles = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
                parseReactionSmiles.addReactant(mol, 1.0);
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.INFO, "Annotating Reaction " + "smiles");
                parseReactionSmiles.setID("smiles");
                return parseReactionSmiles;
            } catch (IllegalArgumentException ex) {
                Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.SEVERE, null, ex);
            }
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    protected IReaction parseMOL2(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        if (!f.isFile()) {
            Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.WARNING, String.format("Mol2 file not found! " + f.getName()));
            System.exit(1);
        }

        String[] split = f.getName().split(".mol");
        MDLV2000Reader mdlV2000Reader = new MDLV2000Reader(
                new FileReader(input), IChemObjectReader.Mode.RELAXED);
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
            Logger.getLogger(ChemicalFormatParser.class.getName()).log(Level.WARNING, String.format("SDF file not found! " + f.getName()));
            System.exit(1);
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
