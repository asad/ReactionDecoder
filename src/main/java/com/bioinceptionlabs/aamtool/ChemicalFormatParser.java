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
package com.bioinceptionlabs.aamtool;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.String.format;
import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.List;
import static java.util.logging.Level.INFO;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.ExtAtomContainerManipulator;
import com.bioinceptionlabs.reactionblast.tools.ChemicalFileIO.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
class ChemicalFormatParser {

    static final String NEW_LINE = getProperty("line.separator");
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(ChemicalFormatParser.class);

    protected static IReaction parseCML(String input) throws FileNotFoundException, CDKException {
        File f = new File(input);
        try {
            f = f.getCanonicalFile();
        } catch (IOException e) {
            throw new FileNotFoundException("Invalid file path: " + input);
        }
        if (!f.isFile()) {
            throw new FileNotFoundException("CML file not found: " + f.getName());
        }
        String[] split = f.getName().split("\\.cml");
        try (FileInputStream fis = new FileInputStream(input);
             CMLReader cmlReader = new CMLReader(fis)) {
            AtomContainer ac = cmlReader.read(new AtomContainer());
            IReaction r = new Reaction();
            r.addReactant(ac, 1.0);
            r.addProduct(ac, 1.0);
            r.setID(split[0]);
            return r;
        } catch (IOException ex) {
            throw new CDKException("Error reading CML file: " + input, ex);
        }
    }

    protected static List<IReaction> parseRXN(String fileNames) {
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
            File filepath;
            try {
                filepath = new File(fileName).getCanonicalFile();
            } catch (IOException e) {
                LOGGER.error(WARNING, format("Invalid file path! %s", fileName));
                continue;
            }
            if (!filepath.isFile()) {
                LOGGER.error(WARNING, format("RXN file not found! %s", filepath.getName()));
                continue;
            }
            LOGGER.info(INFO, "Annotating Reaction {0}", filepath.getName());
            IReaction rxnReactions;
            try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath))) {
                rxnReactions = reader.read(new Reaction());
                rxnReactions.setID(filepath.getName().split("\\.rxn")[0]);
                rxnReactions = convertRoundTripRXNSMILES(rxnReactions);
                reactions.add(rxnReactions);
            } catch (IOException | CDKException ex) {
                LOGGER.error(SEVERE, "ERROR in Reading Reaction file " + filepath, ex);
            }
        }
        return reactions;
    }

    protected static IReaction convertRoundTripRXNSMILES(IReaction ref_reaction) throws CDKException {
        final SmilesGenerator sg = new SmilesGenerator(
                SmiFlavor.AtomAtomMap
                | SmiFlavor.UseAromaticSymbols
                | SmiFlavor.Stereo);
        String createSmilesFromReaction = sg.create(ref_reaction);
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(createSmilesFromReaction);
        parseReactionSmiles.setID(ref_reaction.getID());
        for (int i = 0; i < ref_reaction.getReactantCount(); i++) {
            IAtomContainer atomContainer = parseReactionSmiles.getReactants().getAtomContainer(i);
            String id = ref_reaction.getReactants().getAtomContainer(i).getID();
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
            atomContainer.setID(id);
        }
        for (int i = 0; i < ref_reaction.getProductCount(); i++) {
            IAtomContainer atomContainer = parseReactionSmiles.getProducts().getAtomContainer(i);
            String id = ref_reaction.getProducts().getAtomContainer(i).getID();
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
            atomContainer.setID(id);
        }
        return parseReactionSmiles;
    }

    protected static List<IReaction> parseReactionSMILES(String reactionSmiles) {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String[] smiles = reactionSmiles.split("\\s+");
        List<IReaction> reactions = new ArrayList<>();
        int smilesIndex = 1;
        for (String s : smiles) {
            try {
                IReaction parseReactionSmile = sp.parseReactionSmiles(s);
                try {
                    parseReactionSmile = convertRoundTripRXNSMILES(parseReactionSmile);
                } catch (CDKException e) {
                    LOGGER.error(SEVERE, NEW_LINE, " Sorry - error in Configuring reaction smiles: ", e.getMessage());
                }
                try {
                    LOGGER.info(INFO, "Annotating Reaction " + "smiles");
                    if (smiles.length > 1) {
                        parseReactionSmile.setID("smiles_" + smilesIndex);
                    } else {
                        parseReactionSmile.setID("smiles");
                    }
                    reactions.add(parseReactionSmile);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, NEW_LINE, ex);
                }
            } catch (InvalidSmilesException ex) {
                LOGGER.error(SEVERE, NEW_LINE, ex);
            }
            smilesIndex++;
        }
        return reactions;
    }
}
