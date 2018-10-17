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
package uk.ac.ebi.aamtool;

import java.io.File;
import static java.io.File.separator;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.util.Map;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.aamtool.Annotator.NEW_LINE;
import uk.ac.ebi.reactionblast.tools.ImageGenerator;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000RXNWriter;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class Helper extends ChemicalFormatParser {

    static final String NEW_LINE = getProperty("line.separator");
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Helper.class);

    protected static void getHeader() {
        StringBuilder sb = new StringBuilder();

        sb.append("!--------------------------------------------------------");
        sb.append(NEW_LINE);
        sb.append("Reaction Decoder Tool (RDT)");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Author: Syed Asad Rahman");
        sb.append(NEW_LINE);
        sb.append("e-mail: asad@ebi.ac.uk|s9asad@gmail.com");
        sb.append(NEW_LINE);
        sb.append("c/o EMBL-European BioInformatics Institute (EBI)");
        sb.append(NEW_LINE);
        sb.append("WTGC, CB10 1SD Hinxton");
        sb.append(NEW_LINE);
        sb.append("UK");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Note: The copyright of this software belongs to the author");
        sb.append(NEW_LINE);
        sb.append("and EMBL-European BioInformatics Institute (EBI).");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);

        sb.append("Reference");
        sb.append(NEW_LINE);
        sb.append("Rahman, S.A. et.al.(2016) Reaction Decoder Tool (RDT):");
        sb.append(NEW_LINE);
        sb.append("Extracting Features from Chemical Reactions, Bioinformatics (2016)");
        sb.append(NEW_LINE);
        sb.append("doi: 10.1093/bioinformatics/btw096");
        sb.append(NEW_LINE);
        sb.append("!--------------------------------------------------------");
        sb.append(NEW_LINE);
        out.println(sb.toString());

    }

    /**
     * WreactionWithLayoutite the preactionWithLayoutovided
     * numbereactionWithLayout of blank lineheaderString to the
     * preactionWithLayoutovided OutputStreactionWithLayouteam.
     *
     * @param numberBlankLines NumbereactionWithLayout of blank lineheaderString
     * to wreactionWithLayoutite.
     * @param out OutputStreactionWithLayouteam to which to
     * wreactionWithLayoutite the blank lineheaderString.
     */
    protected static void displayBlankLines(final int numberBlankLines, final OutputStream out) {
        try {
            for (int i = 0; i < numberBlankLines; ++i) {
                out.write(NEW_LINE.getBytes());
            }
        } catch (IOException ioEx) {
            for (int i = 0; i < numberBlankLines; ++i) {
                System.out.println();
            }
        }
    }

    /*
     System.out.println("-- USAGE --");
     printUsage(applicationName + " (Posix)", constructPosixOptions(), System.out);
     displayBlankLines(1, System.out);
     printUsage(applicationName + " (Gnu)", constructGnuOptions(), System.out);
     displayBlankLines(4, System.out);
     System.out.println("-- HELP --");
     printHelp(
     constructPosixOptions(), 80, "POSIX HELP", "End of POSIX Help",
     3, 5, true, System.out);
     displayBlankLines(1, System.out);
     printHelp(
     constructGnuOptions(), 80, "GNU HELP", "End of GNU Help",
     5, 3, true, System.out);
     */
    protected static void printHelp(final OutputStream out, final Options options) {
        final String commandLineSyntax = "java -jar RXNDecoder.jar";
        try (PrintWriter writer = new PrintWriter(out)) {
            final HelpFormatter formatter = new HelpFormatter();
            displayBlankLines(2, out);
            formatter.printHelp(writer, 80, commandLineSyntax, "HELP",
                    options, 5, 3, "End of Helper Help", true);
            writer.flush();
            writer.close();
        }
    }

    protected static void printHelp(final Map<String, Options> optionsMap, final int printedRowWidth,
            final String header, final String footer, final int spacesBeforeOption,
            final int spacesBeforeOptionDescription, final boolean displayUsage, final OutputStream out) {
        final String commandLineSyntax = "java -jar ReactionDecoder.jar";
        try (PrintWriter writer = new PrintWriter(out)) {
            final HelpFormatter helpFormatter = new HelpFormatter();
            optionsMap.keySet().stream().map((headerString) -> {
                helpFormatter.printHelp(
                        writer,
                        printedRowWidth,
                        commandLineSyntax,
                        headerString,
                        optionsMap.get(headerString),
                        spacesBeforeOption,
                        spacesBeforeOptionDescription,
                        "End of Helper " + headerString + " Help",
                        displayUsage);
                return headerString;
            }).map((_item) -> {
                displayBlankLines(2, out);
                return _item;
            }).forEach((_item) -> {
                writer.flush();
            });
            writer.close();
        }
    }

    protected File generateImage(String canonicalRootPath, IReaction mappedReaction, String reactionID) throws Exception {
        File file = new File(canonicalRootPath);
        new ImageGenerator().drawLeftToRightReactionLayout(file, mappedReaction, reactionID);
        return new File(file.getCanonicalFile(), reactionID + ".png");
    }

    protected File generateAAMImage(String canonicalRootPath, IReaction mappedReaction, String reactionID) throws Exception {
        File file = new File(canonicalRootPath);
        new ImageGenerator().drawTopToBottomReactionLayout(file, mappedReaction, reactionID);
        return new File(file.getCanonicalFile(), reactionID + ".png");
    }

    protected File writeRXNMappedFile(String canonicalRootPath, IReaction mappedReaction, String name) throws IOException, CDKException {
//        printReaction(mappedReaction);
        File f = new File(canonicalRootPath + separator + name + ".rxn");
        try (MDLV2000RXNWriter writer = new MDLV2000RXNWriter(new FileWriter(f))) {
            writer.write(mappedReaction);
            writer.close();
        }
        return f;
    }

    protected String printUsageExamples() {
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("Option examples: ");
        sb.append(NEW_LINE);
        sb.append("-Q RXN -q reaction.rxn");
        sb.append(NEW_LINE);
        sb.append("-Q SMI -q \"C=C1CCCC(=O)C1>>O[C@H]1CCCC(=C)C1\"");
        sb.append(NEW_LINE);
        return sb.toString();
    }

}
