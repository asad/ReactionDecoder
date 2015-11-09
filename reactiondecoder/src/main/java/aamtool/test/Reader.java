/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad@ebi.ac.uk>.
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
package aamtool.test;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.imageio.ImageIO;
import javax.vecmath.Vector2d;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockReactionCanoniser;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalReactionLabeller;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */

public class Reader {

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

    public BufferedImage makeBlankImage(int width, int height) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_BGR);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);
        return image;
    }

    public IReaction readReaction(String name) throws IOException {
        return readReaction(name, true);
    }

    public IReaction readReaction(String name, boolean reMap) throws IOException {
        return readReaction(new File(name).getName().split(".rxn")[0], new File(name).getParent(), reMap, false);
    }

    public IReaction readReaction(String name, String dir, boolean reMap) throws IOException {
        return readReaction(name, dir, reMap, false);
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
                Logger.getLogger(Reader.class.getName()).log(Level.WARNING, String.format("RXN file not found! %s", filepath.getName()));
                System.exit(1);
            }
            try {
                Logger.getLogger(Reader.class.getName()).log(Level.INFO, "Annotating Reaction {0}", filepath.getName());
                IReaction rxnReactions;
                try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath));) {
                    try {
                        rxnReactions = reader.read(new Reaction());
                        reader.close();
                        rxnReactions.setID(filepath.getName().split(".rxn")[0]);
                        reactions.add(rxnReactions);
                    } catch (IOException | CDKException ex) {
                        System.err.println("ERROR in Reading Reaction file " + filepath + "\n" + ex);
                    }
                }
            } catch (IOException ex) {
                System.err.println("Failed to Read and Annotate RXN File ");
                Logger.getLogger(Reader.class.getName()).log(Level.SEVERE, null, ex);
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
        String filepath = dir + File.separator + name + ".rxn";

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
//                    System.err.println("ERROR in Mapping " + MolExporter.exportToObject(reaction, "SMARTS"));
//                }
//            } catch (Exception ex) {
//                System.err.println("ERROR in Reading Reaction file " + filepath);
//            }
//        }
//
//        return rxnReactions;
//    }
    public IReaction readReaction(String name, String dir) throws IOException {
        return readReaction(name, dir, false, false);
    }

    public IAtomContainer readMDLMolecule(String name, String dir) throws FileNotFoundException, CDKException {
        String filepath = dir + name + ".mol";
        MDLV2000Reader reader = new MDLV2000Reader(new FileReader(filepath));
        IAtomContainer AtomContainer = (IAtomContainer) reader.read(new AtomContainer());
        AtomContainer.setID(name);
        return AtomContainer;
    }

    public IReaction canonise(IReaction reaction) {
        ICanonicalReactionLabeller labeller = new BlockReactionCanoniser();
        return labeller.getCanonicalReaction(reaction);
    }

    public IReaction map(IReaction reaction) {
        try {
            ReactionMechanismTool rmt = new ReactionMechanismTool(
                    reaction, true, false, false, new StandardizeReaction());
            return rmt.getSelectedSolution().getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        } catch (Exception e) {
            Logger.getLogger(Reader.class.getName()).log(Level.SEVERE, null, e);
            return reaction;
        }
    }

    public void writeToFile(Image image, String outputDirPath, String name) throws IOException {
        File outputDir = new File(outputDirPath);
        if (!outputDir.exists()) {
            outputDir.mkdir();
        }
        File outfile = new File(outputDir, name + ".png");
        if (!outfile.exists()) {
            outfile.createNewFile();
        }
        ImageIO.write((RenderedImage) image, "PNG", outfile);
    }

    public void addImplicitHydrogens(IReaction reaction) {
        addImplicitHydrogens(reaction.getReactants());
        addImplicitHydrogens(reaction.getProducts());
    }

    public void addImplicitHydrogens(IAtomContainerSet molSet) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            addImplicitHydrogens(atomContainer);
        }
    }

    public void addImplicitHydrogens(IAtomContainer atomContainer) {
        try {
            CDKHydrogenAdder.getInstance(
                    DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            Logger.getLogger(Reader.class.getName()).log(Level.SEVERE, null, e);
        }
    }

    public void typeAtoms(IAtomContainer atomContainer) {
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            Logger.getLogger(Reader.class.getName()).log(Level.SEVERE, null, e);
        }
    }

    protected void renumberMappingIDs(IReaction reaction) {
        int i = 1;
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            a0.setID(String.valueOf(i));
            a1.setID(String.valueOf(i));
            mapping.setID(String.valueOf(i));
            i++;
        }
    }

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
        ImageIO.write((RenderedImage) image, "PNG", new File(dir, name + ".png"));
    }
}
