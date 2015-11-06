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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.openscience.cdk.interfaces.IReaction;
import org.w3c.dom.Document;
import static aamtool.Helper.getHeader;
import java.util.List;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class AAMTool extends Annotator {

    private final static boolean DEBUG = false;

    private void FormatXMLToFile(Document doc, String fileName) throws TransformerConfigurationException, TransformerException {

        // write xml to file
        TransformerFactory transformerFactory = TransformerFactory.newInstance();

        Transformer transformer = transformerFactory.newTransformer();
        transformer.setOutputProperty(OutputKeys.METHOD, "xml");
        transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");

        doc.setXmlStandalone(true);
        DOMSource source = new DOMSource(doc);

        /*
         Write to a File
         */
        File file = new File(fileName + ".xml");
        StreamResult result = new StreamResult(file);
        transformer.transform(source, result);

        System.out.println("Output is presented in xml format: " + file.getAbsolutePath());

        if (DEBUG) {
            // Show output on console during development
            result = new StreamResult(System.out);
            transformer.transform(source, result);
        }
    }

    private void FormatTextToFile(StringBuilder doc, String fileName) throws UnsupportedEncodingException, FileNotFoundException, IOException {
        File file = new File(fileName + ".txt");
        try (Writer writer = new OutputStreamWriter(new FileOutputStream(file), "UTF-8")) {
            writer.write(doc.toString());
        }

        System.out.println("Output is presented in text format: " + file.getAbsolutePath());

        if (DEBUG) {
            // Show output on console during development
            System.out.println(doc.toString());
        }
    }

    public AAMTool() {
        super();

    }

    /**
     * @param args the command line areactionWithLayoutgumentheaderString
     */
    public static void main(String[] args) {
        try {
            CommandLineOptions cmd = new CommandLineOptions();
            Options createAAMOptions = cmd.createAAMOptions();

            PosixParser parser1 = new PosixParser();
            CommandLine aamLine = parser1.parse(createAAMOptions, args, true);

            /*
             * Print the Header
             */
            getHeader();

            if (aamLine.hasOption('j') && aamLine.getOptionValue("j").equalsIgnoreCase("AAM")
                    && aamLine.hasOption('Q') && aamLine.hasOption('q') && aamLine.hasOption('f')) {

                System.out.println("-- AAM --");
                AAMTool rxn = new AAMTool();
                rxn.AAMTask(aamLine, createAAMOptions);

            } else if (aamLine.hasOption('j') && aamLine.getOptionValue("j").equalsIgnoreCase("AAM")) {
                System.out.println("-- AAM USAGE --");
                printHelp(System.out, createAAMOptions);
            } else {
                System.out.println("-- AAMTool HELP --");
                Map<String, Options> options = new TreeMap<>();
                options.put("Atom-Atom Mapping (AAM Tool)", createAAMOptions);
                printHelp(options, 80, "EC-BLAST", "End of Help",
                        5, 3, true, System.out);
            }
        } catch (ParseException ex) {
            Logger.getLogger(AAMTool.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(AAMTool.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private void AAMTask(CommandLine aamLine, Options createAAMOptions) throws Exception {

        // TODO code application logic here
        String optionValue = aamLine.getOptionValue("q");

        if (aamLine.hasOption('m')) {
            REPORT_ALL_MAPPINGS = true;
        }

        if (aamLine.hasOption('g')) {
            GENERATE_IMAGE = false;
            GENERATE_AAMIMAGE = true;
        }

        if (aamLine.hasOption('p')) {
            PREFIX = aamLine.getOptionValue("p");
        }

        if (aamLine.hasOption('u')) {
            REMAP = false;
        }

        IReaction reaction = null;

        switch (aamLine.getOptionValue("Q")) {
            case "SMI":
                if (optionValue.contains(">>")) {
                    List<IReaction> parseReactions = parseReactionSMILES(optionValue);
                    if (parseReactions.iterator().hasNext()) {
                        reaction = parseReactions.iterator().next();
                    }
                } else {
                    System.err.println("Not a valid reaction SMILES");
                }
                break;
            case "RXN":
                List<IReaction> parseReactions = parseRXN(optionValue);
                if (parseReactions.iterator().hasNext()) {
                    reaction = parseReactions.iterator().next();
                }
                break;
            default:
                displayBlankLines(2, System.out);
                System.out.println("-- USAGE --");
                printHelp(System.out, createAAMOptions);
                break;
        }
        if (reaction == null) {
            return;
        }

        String jobFileName;

        if (!PREFIX.isEmpty()) {
            jobFileName = PREFIX + "_ECBLAST_" + reaction.getID() + "_AAM";
        } else {
            jobFileName = "ECBLAST_" + reaction.getID() + "_AAM";
        }

        ReactionMechanismTool annotateReaction = getReactionMechanismTool(reaction, REMAP);
        boolean writeFiles = writeFiles(jobFileName, annotateReaction);

        if (writeFiles && aamLine.getOptionValue("f").equalsIgnoreCase("XML")) {
            DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
            // root element
            org.w3c.dom.Document doc = docBuilder.newDocument();

            org.w3c.dom.Element rootElement = doc.createElement("EC_BLAST");
            doc.appendChild(rootElement);
            annotateReactionAsXML(annotateReaction, jobFileName, doc, rootElement);
            FormatXMLToFile(doc, jobFileName);
            System.out.println("XML File saved!");

        } else if (writeFiles && aamLine.getOptionValue("f").equalsIgnoreCase("TEXT")) {
            StringBuilder sb = new StringBuilder();
            annotateReactionAsText(annotateReaction, reaction.getID() + "_AAM", sb);
            FormatTextToFile(sb, jobFileName);
        } else if (writeFiles && aamLine.getOptionValue("f").equalsIgnoreCase("BOTH")) {

            DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
            // root element
            org.w3c.dom.Document doc = docBuilder.newDocument();

            org.w3c.dom.Element rootElement = doc.createElement("EC_BLAST");
            doc.appendChild(rootElement);
            annotateReactionAsXML(annotateReaction, jobFileName, doc, rootElement);

            StringBuilder sb = new StringBuilder();
            annotateReactionAsText(annotateReaction, jobFileName, sb);

            /*
             Write XML and TEXT file
             */
            FormatTextToFile(sb, jobFileName);
            FormatXMLToFile(doc, jobFileName);
            System.out.println("XML File saved!");

        } else {
            displayBlankLines(2, System.out);
            System.out.println("-- USAGE --");
            printHelp(System.out, createAAMOptions);
        }
    }
}
