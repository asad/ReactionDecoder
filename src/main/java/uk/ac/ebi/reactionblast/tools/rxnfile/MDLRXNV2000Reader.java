/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003-2007  The Chemistry Development Kit (CDK) project
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
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import static java.lang.Integer.valueOf;
import static java.lang.System.getProperty;
import java.util.StringTokenizer;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.TITLE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import static org.openscience.cdk.io.IChemObjectReader.Mode.RELAXED;
import org.openscience.cdk.io.formats.IResourceFormat;
import static org.openscience.cdk.io.formats.MDLRXNFormat.getInstance;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Reads a molecule from an MDL RXN file {
 *
 * @cdk.cite DAL92}. This MDL RXN reader uses the MDLV2000 reader to read each
 * mol file
 * @cdk.module io
 *
 *
 * @author Egon Willighagen
 * @author Thomas Kuhn
 * @cdk.created 2003-07-24
 *
 * @cdk.keyword file format, MDL RXN
 * @cdk.bug 1849923
 */
public class MDLRXNV2000Reader extends DefaultChemObjectReader {

    private static ILoggingTool LOGGER = createLoggingTool(MDLRXNV2000Reader.class);
    BufferedReader input = null;

    /**
     * Constructs a new MDLReader that can read AtomContainer from a given
     * Reader.
     *
     * @param in The Reader to read from
     */
    public MDLRXNV2000Reader(Reader in) {
        this(in, RELAXED);
    }

    /**
     *
     * @param in
     * @param mode
     */
    public MDLRXNV2000Reader(Reader in, Mode mode) {
        if (in instanceof BufferedReader) {
            input = (BufferedReader) in;
        } else {
            input = new BufferedReader(in);
        }
        super.mode = mode;
    }

    /**
     *
     * @param input
     */
    public MDLRXNV2000Reader(InputStream input) {
        this(input, RELAXED);
    }

    /**
     *
     * @param input
     * @param mode
     */
    public MDLRXNV2000Reader(InputStream input, Mode mode) {
        this(new InputStreamReader(input), mode);
    }

    /**
     *
     */
    public MDLRXNV2000Reader() {
        this(new StringReader(""));
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
     * @param input
     * @throws CDKException
     */
    @Override
    public void setReader(Reader input) throws CDKException {
        if (input instanceof BufferedReader) {
            this.input = (BufferedReader) input;
        } else {
            this.input = new BufferedReader(input);
        }
    }

    /**
     *
     * @param input
     * @throws CDKException
     */
    @Override
    public void setReader(InputStream input) throws CDKException {
        setReader(new InputStreamReader(input));
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
            if (IChemModel.class.equals(intf)) {
                return true;
            }
            if (IChemFile.class.equals(intf)) {
                return true;
            }
            if (IReaction.class.equals(intf)) {
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
     * Takes an object which subclasses IChemObject, e.g.AtomContainer, and will
     * read this (from file, database, internet etc). If the specific
     * implementation does not support a specific IChemObject it will throw an
     * Exception.
     *
     * @param <T>
     * @param object The object that subclasses IChemObject
     * @return The IChemObject read
     * @exception CDKException
     */
    @Override
    public <T extends IChemObject> T read(T object) throws CDKException {
        if (object instanceof IReaction) {
            return (T) readReaction(object.getBuilder());
        } else if (object instanceof IReactionSet) {
            IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
            reactionSet.addReaction(readReaction(object.getBuilder()));
            return (T) reactionSet;
        } else if (object instanceof IChemModel) {
            IChemModel model = object.getBuilder().newInstance(IChemModel.class);
            IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
            reactionSet.addReaction(readReaction(object.getBuilder()));
            model.setReactionSet(reactionSet);
            return (T) model;
        } else if (object instanceof IChemFile) {
            IChemFile chemFile = object.getBuilder().newInstance(IChemFile.class);
            IChemSequence sequence = object.getBuilder().newInstance(IChemSequence.class);
            sequence.addChemModel(read(object.getBuilder().newInstance(IChemModel.class)));
            chemFile.addChemSequence(sequence);
            return (T) chemFile;
        } else {
            throw new CDKException("Only supported are Reaction and ChemModel, and not "
                    + object.getClass().getName() + ".");
        }
    }

    /**
     *
     * @param object
     * @return
     */
    public boolean accepts(IChemObject object) {
        if (object instanceof IReaction) {
            return true;
        } else if (object instanceof IChemModel) {
            return true;
        } else if (object instanceof IChemFile) {
            return true;
        } else if (object instanceof IReactionSet) {
            return true;
        }
        return false;
    }

    /**
     * Read a Reaction from a file in MDL RXN format
     *
     * @return The Reaction that was read from the MDL file.
     */
    private IReaction readReaction(IChemObjectBuilder builder) throws CDKException {
        IReaction reaction = builder.newInstance(IReaction.class);
        try {
            input.readLine(); // first line should be $RXN
            input.readLine(); // second line
            input.readLine(); // third line
            input.readLine(); // fourth line
        } catch (IOException exception) {
            LOGGER.debug(exception);
            throw new CDKException("Error while reading header of RXN file", exception);
        }

        int reactantCount = 0;
        int productCount = 0;
        try {
            String countsLine = input.readLine();
            /* this line contains the number of reactants
             and products */
            StringTokenizer tokenizer = new StringTokenizer(countsLine);
            reactantCount = valueOf(tokenizer.nextToken());
            LOGGER.info("Expecting " + reactantCount + " reactants in file");
            productCount = valueOf(tokenizer.nextToken());
            LOGGER.info("Expecting " + productCount + " products in file");
        } catch (IOException | NumberFormatException exception) {
            LOGGER.debug(exception);
            throw new CDKException("Error while counts line of RXN file", exception);
        }

        // now read the reactants
        try {
            for (int i = 1; i <= reactantCount; i++) {
                StringBuilder molFile = new StringBuilder();
                input.readLine(); // announceMDLFileLine
                String molFileLine = "";
                do {
                    molFileLine = input.readLine();
                    molFile.append(molFileLine);
                    molFile.append(getProperty("line.separator"));
                } while (!molFileLine.equals("M  END"));

                // read MDL molfile content
                // Changed this to mdlv2000 reader
                MDLV2000Reader reader = new MDLV2000Reader(
                        new StringReader(molFile.toString()),
                        super.mode);
                IAtomContainer reactant = reader.read(
                        builder.newInstance(IAtomContainer.class));
                if (reactant == null) {
                    continue;
                }
                // add reactant mol ID
                String readMolID = (String) reactant.getProperty(TITLE);
                if (readMolID != null) {
                    reactant.setID(readMolID.trim());
                }
                // add reactant
                reaction.addReactant(reactant);
            }
        } catch (CDKException exception) {
            // rethrow exception from MDLReader
            throw exception;
        } catch (IOException | IllegalArgumentException exception) {
            LOGGER.debug(exception);
            throw new CDKException("Error while reading reactant", exception);
        }

        // now read the products
        try {
            for (int i = 1; i <= productCount; i++) {
                StringBuilder molFile = new StringBuilder();
                input.readLine(); // String announceMDLFileLine = 
                String molFileLine = "";
                do {
                    molFileLine = input.readLine();
                    molFile.append(molFileLine);
                    molFile.append(getProperty("line.separator"));
                } while (!molFileLine.equals("M  END"));

                // read MDL molfile content
                MDLV2000Reader reader = new MDLV2000Reader(
                        new StringReader(molFile.toString()));
                IAtomContainer product = reader.read(
                        builder.newInstance(IAtomContainer.class));

                if (product == null) {
                    continue;
                }

                // add product molID
                String readMolID = (String) product.getProperty(TITLE);
                if (readMolID != null) {
                    product.setID(readMolID.trim());
                }
                // add product
                reaction.addProduct(product);
            }
        } catch (CDKException exception) {
            // rethrow exception from MDLReader
            throw exception;
        } catch (IOException | IllegalArgumentException exception) {
            LOGGER.debug(exception);
            throw new CDKException("Error while reading products", exception);
        }

        // now try to map things, if wanted
        LOGGER.info("Reading atom-atom mapping from file");
        // distribute all atoms over two GraphAtomContainer's
        IAtomContainer reactingSide = builder.newInstance(IAtomContainer.class);
        java.util.Iterator<IAtomContainer> molecules = reaction.getReactants().atomContainers().iterator();
        while (molecules.hasNext()) {
            reactingSide.add(molecules.next());
        }
        IAtomContainer producedSide = builder.newInstance(IAtomContainer.class);
        molecules = reaction.getProducts().atomContainers().iterator();
        while (molecules.hasNext()) {
            producedSide.add(molecules.next());
        }

        // map the atoms
        int mappingCount = 0;
        for (int i = 0; i < reactingSide.getAtomCount(); i++) {
            for (int j = 0; j < producedSide.getAtomCount(); j++) {
                IAtom eductAtom = reactingSide.getAtom(i);
                IAtom productAtom = producedSide.getAtom(j);
                if (eductAtom.getProperty(ATOM_ATOM_MAPPING) != null
                        && eductAtom.getProperty(ATOM_ATOM_MAPPING).equals(productAtom.getProperty(ATOM_ATOM_MAPPING))) {
                    reaction.addMapping(
                            builder.newInstance(IMapping.class, eductAtom, productAtom));
                    mappingCount++;
                    break;
                }
            }
        }
        LOGGER.info("Mapped atom pairs: " + mappingCount);

        return reaction;
    }

    @Override
    public void close() throws IOException {
        input.close();
    }
}
