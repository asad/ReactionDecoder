/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad at ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.tools;

import static java.awt.Color.WHITE;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Image;
import static java.awt.RenderingHints.KEY_ANTIALIASING;
import static java.awt.RenderingHints.VALUE_ANTIALIAS_ON;
import java.awt.image.BufferedImage;
import static java.awt.image.BufferedImage.TYPE_INT_BGR;
import java.awt.image.RenderedImage;
import java.io.File;
import static java.io.File.separator;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.String.format;
import static java.lang.String.valueOf;
import static java.lang.System.exit;
import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.List;
import static java.util.logging.Level.INFO;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;

import static javax.imageio.ImageIO.write;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainer;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import static org.openscience.cdk.tools.CDKHydrogenAdder.getInstance;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockReactionCanoniser;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalReactionLabeller;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Reader {

     private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Reader.class);
     static final String NEW_LINE = getProperty("line.separator");

    /**
     *
     * @param AtomContainer
     * @return
     */
    public IAtomContainer layout(IAtomContainer AtomContainer) {
        try {
            StructureDiagramGenerator sdg = new StructureDiagramGenerator();
            sdg.setMolecule(AtomContainer, true);
            sdg.generateCoordinates();
            return sdg.getMolecule();
        } catch (CDKException e) {
            return AtomContainer;
        }
    }

    /**
     *
     * @param width
     * @param height
     * @return
     */
    public BufferedImage makeBlankImage(int width, int height) {
        BufferedImage image = new BufferedImage(width, height, TYPE_INT_BGR);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        return image;
    }

    /**
     *
     * @param name
     * @return
     * @throws IOException
     */
    public IReaction readReaction(String name) throws IOException {
        return readReaction(name, true);
    }

    /**
     *
     * @param name
     * @param reMap
     * @return
     * @throws IOException
     */
    public IReaction readReaction(String name, boolean reMap) throws IOException {
        return readReaction(new File(name).getName().split(".rxn")[0], new File(name).getParent(), reMap, false);
    }

    /**
     *
     * @param name
     * @param dir
     * @param reMap
     * @return
     * @throws IOException
     */
    public IReaction readReaction(String name, String dir, boolean reMap) throws IOException {
        return readReaction(name, dir, reMap, false);
    }

    /**
     *
     * @param fileNames
     * @return
     */
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
                LOGGER.error(WARNING, format("RXN file not found! %s", filepath.getName()));
                exit(1);
            }
            try {
                LOGGER.error(INFO, "Annotating Reaction {0}", filepath.getName());
                IReaction rxnReactions;
                try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath));) {
                    try {
                        rxnReactions = reader.read(new Reaction());
                        reader.close();
                        rxnReactions.setID(filepath.getName().split(".rxn")[0]);
                        reactions.add(rxnReactions);
                    } catch (IOException | CDKException ex) {
                        LOGGER.debug("ERROR in Reading Reaction file " + filepath + NEW_LINE + ex);
                    }
                }
            } catch (IOException ex) {
                LOGGER.debug("Failed to Read and Annotate RXN File ");
                LOGGER.error(SEVERE, null, ex);
            }
        }
        return reactions;
    }

    /**
     *
     * @param name
     * @param dir
     * @param reMap
     * @param removeHydrogens
     * @return
     * @throws IOException
     */
    public IReaction readReaction(String name, String dir, boolean reMap, boolean removeHydrogens) throws IOException {
        String filepath = dir + separator + name + ".rxn";

        List<IReaction> parseRXN = parseRXN(filepath);
        IReaction rxnReactions = parseRXN.iterator().hasNext() ? parseRXN.iterator().next() : null;
        if (reMap) {
            rxnReactions = map(rxnReactions);
            renumberMappingIDs(rxnReactions);
        }
        return rxnReactions;
    }

//    public IReaction readReaction(String name, String dir, boolean reMap, boolean removeHydrogens) throws IOException {
//        String filepath = dir + File.separator + name + ".rxn";
//        IReaction rxnReactions = null;
//        try (chemaxon.formats.MolImporter me = new chemaxon.formats.MolImporter(filepath)) {
//            chemaxon.struc.Molecule mol = me.read();
//            RxnMolecule reaction = chemaxon.struc.RxnMolecule.getReaction(mol);
//            reaction.valenceCheck();
//            me.close();
//            try {
//                reaction.calcHybridization();
//                Hydrogenize.convertImplicitHToExplicit(reaction);
//                rxnReactions = CDKChemaxonIOConveter.getCDKReaction(reaction);
//                try {
//                    rxnReactions.setID(name);
//                    if (reMap) {
//                        rxnReactions = map(rxnReactions);
//                        renumberMappingIDs(rxnReactions);
//                    }
//                } catch (Exception e) {
//                    System.LOGGER.debug("ERROR in Mapping " + MolExporter.exportToObject(reaction, "SMARTS"));
//                }
//            } catch (Exception ex) {
//                System.LOGGER.debug("ERROR in Reading Reaction file " + filepath);
//            }
//        }
//
//        return rxnReactions;
//    }
    /**
     *
     * @param name
     * @param dir
     * @return
     * @throws IOException
     */
    public IReaction readReaction(String name, String dir) throws IOException {
        return readReaction(name, dir, false, false);
    }

    /**
     *
     * @param name
     * @param dir
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IAtomContainer readMDLMolecule(String name, String dir) throws FileNotFoundException, CDKException {
        String filepath = dir + name + ".mol";
        MDLV2000Reader reader = new MDLV2000Reader(new FileReader(filepath));
        IAtomContainer AtomContainer = reader.read(new AtomContainer());
        AtomContainer.setID(name);
        return AtomContainer;
    }

    /**
     *
     * @param reaction
     * @return
     */
    public IReaction canonise(IReaction reaction) {
        ICanonicalReactionLabeller labeller = new BlockReactionCanoniser();
        return labeller.getCanonicalReaction(reaction);
    }

    /**
     *
     * @param reaction
     * @return
     */
    public IReaction map(IReaction reaction) {
        try {
            ReactionMechanismTool rmt = new ReactionMechanismTool(
                    reaction, true, false, false, new StandardizeReaction());
            return rmt.getSelectedSolution().getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        } catch (Exception e) {
            LOGGER.error(SEVERE, null, e);
            return reaction;
        }
    }

    /**
     *
     * @param image
     * @param outputDirPath
     * @param name
     * @throws IOException
     */
    public void writeToFile(Image image, String outputDirPath, String name) throws IOException {
        File outputDir = new File(outputDirPath);
        if (!outputDir.exists()) {
            outputDir.mkdir();
        }
        File outfile = new File(outputDir, name + ".png");
        if (!outfile.exists()) {
            outfile.createNewFile();
        }
        write((RenderedImage) image, "PNG", outfile);
    }

    /**
     *
     * @param reaction
     */
    public void addImplicitHydrogens(IReaction reaction) {
        addImplicitHydrogens(reaction.getReactants());
        addImplicitHydrogens(reaction.getProducts());
    }

    /**
     *
     * @param molSet
     */
    public void addImplicitHydrogens(IAtomContainerSet molSet) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            addImplicitHydrogens(atomContainer);
        }
    }

    /**
     *
     * @param atomContainer
     */
    public void addImplicitHydrogens(IAtomContainer atomContainer) {
        try {
            getInstance(getInstance()).addImplicitHydrogens(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            LOGGER.error(SEVERE, null, e);
        }
    }

    /**
     *
     * @param atomContainer
     */
    public void typeAtoms(IAtomContainer atomContainer) {
        try {
            percieveAtomTypesAndConfigureAtoms(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            LOGGER.error(SEVERE, null, e);
        }
    }

    /**
     *
     * @param reaction
     */
    protected void renumberMappingIDs(IReaction reaction) {
        int i = 1;
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            a0.setID(valueOf(i));
            a1.setID(valueOf(i));
            mapping.setID(valueOf(i));
            i++;
        }
    }

    /**
     *
     * @param molDrawer
     * @param AtomContainer
     * @param dirPath
     * @param name
     * @throws IOException
     */
    public void draw(DirectMoleculeDrawer molDrawer, IAtomContainer AtomContainer,
            String dirPath, String name) throws IOException {
        int width = 700;
        int height = 700;
        SingleMoleculeLayout layout = new SingleMoleculeLayout(molDrawer.getParams());
        layout.layout(AtomContainer, new Vector2d(width / 2, height / 2));
        ZoomToFitLayout layout2 = new ZoomToFitLayout(molDrawer);
        Image image = makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
//        if (molDrawer.getParams().useAntialias) {
//            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
//                               RenderingHints.VALUE_ANTIALIAS_ON);
//        }
//        molDrawer.drawAtomContainer(AtomContainer, g);
        layout2.layout(AtomContainer, new Dimension(width, height), g);
        File dir = new File(dirPath);
        if (!dir.exists()) {
            dir.mkdir();
        }
        write((RenderedImage) image, "PNG", new File(dir, name + ".png"));
    }
}
