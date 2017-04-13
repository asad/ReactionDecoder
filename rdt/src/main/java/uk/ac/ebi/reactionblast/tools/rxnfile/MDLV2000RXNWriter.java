/* $RCSfile$
 * $Author$ 
 * $Date$
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.text.NumberFormat;
import static java.text.NumberFormat.getNumberInstance;
import java.util.Iterator;
import static java.util.Locale.ENGLISH;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.REMARK;
import static org.openscience.cdk.CDKConstants.TITLE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.io.DefaultChemObjectWriter;
import org.openscience.cdk.io.formats.IResourceFormat;
import static org.openscience.cdk.io.formats.MDLFormat.getInstance;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Writes a reaction to a MDL rxn or SDF file. Attention: Stoichiometric
 * coefficients have to be natural numbers.
 *
 * <pre>
 * MDLRXNWriter writer = new MDLRXNWriter(new FileWriter(new File("output.mol")));
 * writer.write((AtomContainer)molecule);
 * writer.close();
 * </pre>
 *
 * See {
 *
 * @cdk.cite DAL92}.
 *
 * @cdk.module io
 * 
 *
 * @cdk.keyword file format, MDL RXN file
 */

public class MDLV2000RXNWriter extends DefaultChemObjectWriter {

    private static ILoggingTool logger
            = createLoggingTool(MDLV2000RXNWriter.class);
    private static final Logger LOG = getLogger(MDLV2000RXNWriter.class.getName());
    private BufferedWriter writer;
    private int reactionNumber;

    /**
     *
     */
    public Map rdFields = null;

    /**
     * Constructs a new MDLWriter that can write an array of Molecules to a
     * Writer.
     *
     * @param out The Writer to write to
     */
    public MDLV2000RXNWriter(Writer out) {
        try {
            if (out instanceof BufferedWriter) {
                writer = (BufferedWriter) out;
            } else {
                writer = new BufferedWriter(out);
            }
        } catch (Exception exc) {
        }
        this.reactionNumber = 1;
    }

    /**
     * Constructs a new MDLWriter that can write an array of Molecules to a
     * given OutputStream.
     *
     * @param output The OutputStream to write to
     */
    public MDLV2000RXNWriter(OutputStream output) {
        this(new OutputStreamWriter(output));
    }

    /**
     *
     */
    public MDLV2000RXNWriter() {
        this(new StringWriter());
    }

    /**
     *
     * @return
     */
    @Override
    public IResourceFormat getFormat() {
        return getInstance();
    }

    /**
     *
     * @param out
     * @throws CDKException
     */
    @Override
    public void setWriter(Writer out) throws CDKException {
        if (out instanceof BufferedWriter) {
            writer = (BufferedWriter) out;
        } else {
            writer = new BufferedWriter(out);
        }
    }

    /**
     *
     * @param output
     * @throws CDKException
     */
    @Override
    public void setWriter(OutputStream output) throws CDKException {
        setWriter(new OutputStreamWriter(output));
    }

    /**
     * Here you can set a map which will be used to build rd fields in the file.
     * The entries will be translated to rd fields like this:<br>
     * &gt; &lt;key&gt;<br>
     * &gt; value<br>
     * empty line<br>
     *
     * @param map The map to be used, map of String-String pairs
     */
    public void setRdFields(Map map) {
        rdFields = map;
    }

    /**
     * Flushes the output and closes this object.
     *
     * @throws java.io.IOException
     */
    @Override
    public void close() throws IOException {
        writer.close();
    }

    /**
     *
     * @param classObject
     * @return
     */
    @Override
    public boolean accepts(Class classObject) {
        Class[] interfaces = classObject.getInterfaces();
        for (Class intf : interfaces) {
            if (IReaction.class.equals(intf)) {
                return true;
            }
            if (IReactionSet.class.equals(intf)) {
                return true;
            }
        }
        Class superClass = classObject.getSuperclass();
        if (superClass != null) {
            return this.accepts(superClass);
        }
        return false;
    }

    /**
     * Writes a IChemObject to the MDL RXN file formated output. It can only
     * output ChemObjects of type Reaction
     *
     * @param object class must be of type AtomContainer or MoleculeSet.
     * @throws org.openscience.cdk.exception.CDKException
     *
     * @see org.openscience.cdk.ChemFile
     */
    @Override
    public void write(IChemObject object) throws CDKException {
        if (object instanceof IReactionSet) {
            writeReactionSet((IReactionSet) object);
        } else if (object instanceof IReaction) {
            writeReaction((IReaction) object);
        } else {
            throw new CDKException("Only supported is writing ReactionSet, Reaction objects.");
        }
    }

    /**
     * Writes an array of Reaction to an OutputStream in MDL rdf format.
     *
     * @param som Array of Reactions that is written to an OutputStream
     */
    private void writeReactionSet(IReactionSet reactions) throws CDKException {

        for (Iterator<IReaction> it = reactions.reactions().iterator(); it.hasNext();) {
            writeReaction(it.next());
        }
    }

    /**
     * Writes a Reaction to an OutputStream in MDL sdf format.
     *
     * @param reaction A Reaction that is written to an OutputStream
     */
    private void writeReaction(IReaction reaction) throws CDKException {

        /*Fixed correct reactant product count*/
        int reactantCount = 0;
        for (IAtomContainer e : reaction.getReactants().atomContainers()) {
            reactantCount += reaction.getReactantCoefficient(e).intValue();
        }
        int productCount = 0;
        for (IAtomContainer p : reaction.getProducts().atomContainers()) {
            productCount += reaction.getProductCoefficient(p).intValue();
        }

        if (reactantCount <= 0 || productCount <= 0) {
            throw new CDKException("Either no reactants or no products present.");
        }

        try {
            // taking care of the $$$$ signs:
            // we do not write such a sign at the end of the first reaction, thus we have to write on BEFORE the second reaction
            if (reactionNumber == 2) {
                writer.write("$$$$");
                writer.newLine();
            }
            writer.write("$RXN");
            writer.newLine();

            // reaction name
            String line = (String) reaction.getProperty(TITLE);
            if (line == null) {
                String rid = reaction.getID() == null ? "" : reaction.getID();
                line = "  " + "EC-BLAST" + "     " + rid;
            }
            if (line.length() > 80) {
                line = line.substring(0, 80);
            }
            writer.newLine();
            writer.write(line);
            // user/program/date&time/reaction registry no. line
            writer.newLine();
            // comment line
            line = (String) reaction.getProperty(REMARK);
            if (line == null) {
                line = "";
            }
            if (line.length() > 80) {
                line = line.substring(0, 80);
            }
            writer.write(line);
            writer.newLine();

            line = "";
            line += formatMDLInt(reactantCount, 3);
            line += formatMDLInt(productCount, 3);
            writer.write(line);
            writer.newLine();

            int i = 0;
            for (IMapping mapping : reaction.mappings()) {
                Iterator<IChemObject> it = mapping.relatedChemObjects().iterator();
                /*
                 Do not overwrite the existing labels
                 */

                if (it.next().getProperty(ATOM_ATOM_MAPPING) == null) {
                    it.next().setProperty(ATOM_ATOM_MAPPING, i + 1);
                    it.next().setProperty(ATOM_ATOM_MAPPING, i + 1);
                    i++;
                }
            }
            writeMoleculeSet(reaction.getReactants());
            writeMoleculeSet(reaction.getProducts());

            //write sdfields, if any
            if (rdFields != null) {
                Set set = rdFields.keySet();
                Iterator iterator = set.iterator();
                while (iterator.hasNext()) {
                    Object element = iterator.next();
                    writer.write("> <" + element + ">");
                    writer.newLine();
                    writer.write(rdFields.get(element).toString());
                    writer.newLine();
                    writer.newLine();
                }
            }
            // taking care of the $$$$ signs:
            // we write such a sign at the end of all except the first molecule
            if (reactionNumber != 1) {
                writer.write("$$$$");
                writer.newLine();
            }
            reactionNumber++;

        } catch (IOException ex) {
            logger.error(ex.getMessage());
            logger.debug(ex);
            throw new CDKException("Exception while writing MDL file: " + ex.getMessage(), ex);
        }
    }

    /**
     * Writes a MoleculeSet to an OutputStream for the reaction.
     *
     * @param som The MoleculeSet that is written to an OutputStream
     */
    private void writeMoleculeSet(IAtomContainerSet som) throws IOException, CDKException {

        for (int i = 0; i < som.getAtomContainerCount(); i++) {
            IAtomContainer mol = som.getAtomContainer(i);
            for (int j = 0; j < som.getMultiplier(i); j++) {
                StringWriter sw = new StringWriter();
                writer.write("$MOL");
                writer.newLine();
                MDLV2000Writer mdlwriter = null;
                try {
                    mdlwriter = new MDLV2000Writer(sw);
                    // GMT : added to force writing of aromatic bond types
                    // Asad: Swith off aromatic bond writing
                    // mdlwriter.getIOSettings()[1].setSetting("true");
                } catch (Exception ex) {
                    logger.error(ex.getMessage());
                    logger.debug(ex);
                    throw new CDKException("Exception while creating MDLWriter: " + ex.getMessage(), ex);
                }
                mdlwriter.write(mol);
                writer.write(sw.toString());
            }
        }
    }

    /**
     * Formats an int to fit into the connection table and changes it to a
     * String.
     *
     * @param i The int to be formated
     * @param l Length of the String
     * @return The String to be written into the connection table
     */
    private String formatMDLInt(int i, int l) {
        String s = "", fs = "";
        NumberFormat nf = getNumberInstance(ENGLISH);
        nf.setParseIntegerOnly(true);
        nf.setMinimumIntegerDigits(1);
        nf.setMaximumIntegerDigits(l);
        nf.setGroupingUsed(false);
        s = nf.format(i);
        l -= s.length();
        for (int f = 0; f < l; f++) {
            fs += " ";
        }
        fs += s;
        return fs;
    }
}
