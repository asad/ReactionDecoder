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
package mapping;

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
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockReactionCanoniser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.ImageGenerator;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalReactionLabeller;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BaseTest {

    final static String KEGG_RXN_DIR = "rxn/";
    final static String Rhea_RXN_DIR = "rxn/rhea/";
    final static String Brenda_RXN_DIR = "rxn/brenda/";
    final static String BUG_RXN_DIR = "rxn/bug/";
    final static String OTHER_RXN = "rxn/other/";
    final static String METRXN_RXN = "rxn/metrxn/";
    final static String INFORCHEM_RXN = "rxn/infochem/";
    final static String MACIE_RXN = "rxn/macie/";

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

    /**
     *
     * @param name
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IReaction readReaction(String name) throws FileNotFoundException, CDKException, Exception {
        return readReaction(name, true);
    }

    /**
     *
     * @param name
     * @param reMap
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IReaction readReaction(String name, boolean reMap) throws FileNotFoundException, CDKException, Exception {
        return readReaction(name, "Output/", reMap, false);
    }

    /**
     *
     * @param name
     * @param dir
     * @param reMap
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IReaction readReaction(String name, String dir, boolean reMap) throws FileNotFoundException, CDKException, Exception {
        return readReaction(name, dir, reMap, false);
    }

    /**
     *
     * @param name
     * @param dir
     * @param reMap
     * @param removeHydrogens
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IReaction readReaction(String name, String dir, boolean reMap, boolean removeHydrogens) throws Exception {
        String filepath = dir + name + ".rxn";

        IReaction reaction = null;
        try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(filepath))) {
            reaction = (IReaction) reader.read(new Reaction());
            reaction.setID(name);
        } catch (Exception ex) {
            Logger.getLogger(BaseTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (reMap) {
            reaction = map(reaction);
        }

        if (removeHydrogens) {
            // XXX WARNING : this may not work correctly!
            IReaction hydrogenFreeReaction = new Reaction();
            IAtomContainerSet reactants = new AtomContainerSet();
            for (IAtomContainer atomContainer : reaction.getReactants().atomContainers()) {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = ExtAtomContainerManipulator.convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID((String) atomContainer.getProperty(CDKConstants.TITLE));
                reactants.addAtomContainer(acMinusH);
            }
            hydrogenFreeReaction.setReactants(reactants);
            IAtomContainerSet products = new AtomContainerSet();
            for (IAtomContainer atomContainer : reaction.getProducts().atomContainers()) {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = ExtAtomContainerManipulator.convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID((String) atomContainer.getProperty(CDKConstants.TITLE));
                products.addAtomContainer(acMinusH);
            }
            hydrogenFreeReaction.setProducts(products);
            for (IMapping mapping : reaction.mappings()) {
                if (((IAtom) mapping.getChemObject(0)).getSymbol().equals("H")
                        || ((IAtom) mapping.getChemObject(1)).getSymbol().equals("H")) {
                    continue;
                }
                hydrogenFreeReaction.addMapping(mapping);
            }
            reaction = hydrogenFreeReaction;
        }
        renumberMappingIDs(reaction);
        return reaction;
    }

    /**
     *
     * @param name
     * @param dir
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public IReaction readReaction(String name, String dir) throws FileNotFoundException, CDKException, Exception {
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
        IAtomContainer AtomContainer = (IAtomContainer) reader.read(new AtomContainer());
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
                    //                    reaction, true, true, true, new CDKReactionStandardizer());
                    reaction, true, false, false, new StandardizeReaction());
            return rmt.getSelectedSolution().getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        } catch (Exception e) {
            e.printStackTrace();
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
        ImageIO.write((RenderedImage) image, "PNG", outfile);
    }

    private void setLonePairs(IReaction reaction) {
        LonePairElectronChecker checker = new LonePairElectronChecker();
        setLonePairs(reaction.getReactants(), checker);
        setLonePairs(reaction.getProducts(), checker);
    }

    private void setLonePairs(IAtomContainerSet AtomContainerSet, LonePairElectronChecker checker) {
        for (IAtomContainer atomContainer : AtomContainerSet.atomContainers()) {
            try {
                checker.saturate(atomContainer);
            } catch (CDKException c) {
                c.printStackTrace();
            }
        }
    }

    private void detectAromaticity(IReaction reaction) {
        detectAromaticity(reaction.getReactants());
        detectAromaticity(reaction.getProducts());
    }

    private void detectAromaticity(IAtomContainerSet molSet) {
        for (IAtomContainer ac : molSet.atomContainers()) {
            try {
                ExtAtomContainerManipulator.aromatizeCDK(ac);
            } catch (CDKException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
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
            e.printStackTrace();
        }
    }

    public void typeAtoms(IAtomContainer atomContainer) {
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private void setMappingIDs(IReaction reaction) {
        int i = 0;
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            Object mappingID = a0.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
//            Integer mappingID =
//                (Integer)a0.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
            if (mappingID != null) {
                a0.setID(String.valueOf(mappingID));
                a1.setID(String.valueOf(mappingID));
            } else {
                a0.setID(String.valueOf(i));
                a1.setID(String.valueOf(i));
            }
            i++;
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

    /**
     *
     * @param reactionID
     * @param directory
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public ReactionMechanismTool testReactions(String reactionID, String directory) throws FileNotFoundException, Exception {
        IReaction cdkReaction = readReaction(reactionID, directory, false);
        SmilesGenerator withAtomClasses = SmilesGenerator.unique().aromatic().withAtomClasses();
        System.out.println("Input reactions " + withAtomClasses.createReactionSMILES(cdkReaction));
        /*
         RMT for the reaction mapping
         */
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, true, true, false, new StandardizeReaction());
        MappingSolution s = rmt.getSelectedSolution();

        System.out.println("Mapped reactions " + withAtomClasses.createReactionSMILES(s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens()));

        System.out.println("Reaction ID: " + reactionID + ", Selected Algorithm: " + s.getAlgorithmID());
//        System.out.println("Cleaved/Formed " + s.getBondChangeCalculator().getFormedCleavedWFingerprint().toString());
//        System.out.println("Order Changed " + s.getBondChangeCalculator().getOrderChangesWFingerprint().toString());
//        System.out.println("Stereo Changed " + s.getBondChangeCalculator().getStereoChangesWFingerprint().toString());
//        System.out.println("RC Changed " + s.getBondChangeCalculator().getReactionCenterWFingerprint().toString());
//        System.out.println("BE " + s.getBondEnergyChange() + ", Fragment " + s.getTotalFragmentChanges());
        IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        new ImageGenerator().drawTopToBottomReactionLayout("Output", reactionWithCompressUnChangedHydrogens, (reactionID + s.getAlgorithmID()));

        int i = 1;
        for (MappingSolution m : rmt.getAllSolutions()) {
            System.out.println("--------------------------------------");
            BondChangeCalculator bcc = m.getBondChangeCalculator();
            System.out.println(m.getAlgorithmID() + ", fp " + bcc.getFormedCleavedWFingerprint().toString());
//            System.out.println(m.getAlgorithmID() + ", fp " + bcc.getOrderChangesWFingerprint().toString());

            System.out.println("BE " + m.getBondEnergySum() + ", Fragment " + m.getTotalFragmentChanges());
            new ImageGenerator().drawLeftToRightReactionLayout("Output", bcc.getReactionWithCompressUnChangedHydrogens(), ("Map_" + reactionID + m.getAlgorithmID()));
//            ImageGenerator.drawTopToBottomReactionLayout(bcc.getReactionWithCompressUnChangedHydrogens(), ("Map_" + reactionID + m.getAlgorithmID()));
            i++;
            System.out.println();
            System.out.println("--------------------------------------");
        }
        return rmt;

    }

    /**
     *
     * @param reactionID
     * @param directory
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public BondChangeCalculator testRCReactions(String reactionID, String directory) throws FileNotFoundException, Exception {
        String NEW_LINE = System.getProperty("line.separator");
        IReaction cdkReaction = readReaction(reactionID, directory, false);
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, true, true, false);
        MappingSolution s = rmt.getSelectedSolution();
        new ImageGenerator().drawLeftToRightReactionLayout("Output", s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens(), (reactionID + s.getAlgorithmID()));

        StringBuilder sb = new StringBuilder();
        sb.append("++++++++++++++++++++++++++++++++++++++++++");
        sb.append(NEW_LINE);
        sb.append("Bond Changes:");
        sb.append(NEW_LINE);
        sb.append("//");
        sb.append(NEW_LINE);
        sb.append(s.getBondChangeCalculator().getFormedCleavedWFingerprint().toString());
        sb.append(NEW_LINE);
        sb.append(s.getBondChangeCalculator().getOrderChangesWFingerprint().toString());
        sb.append(NEW_LINE);
        sb.append(s.getBondChangeCalculator().getStereoChangesWFingerprint().toString());
        sb.append(NEW_LINE);
        sb.append(s.getBondChangeCalculator().getReactionCenterWFingerprint().toString());
        sb.append(NEW_LINE);
        sb.append("//");
        sb.append(NEW_LINE);
        System.out.println(sb.toString());

        StringBuilder rcSteps = new StringBuilder();
        rcSteps.append("Formed Cleaved");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterFormedCleavedFingerprint = s.getBondChangeCalculator().getReactionCenterFormedCleavedFingerprint();
        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterFormedCleavedFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            rcSteps.append(m.getValue());
        }
        rcSteps.append("Order Change");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterOrderChangeFingerprint = s.getBondChangeCalculator().getReactionCenterOrderChangeFingerprint();
        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterOrderChangeFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            rcSteps.append(m.getValue());
        }
        rcSteps.append("Stereo Change");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterStereoChangeFingerprint = s.getBondChangeCalculator().getReactionCenterStereoChangeFingerprint();
        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterStereoChangeFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            rcSteps.append(m.getValue());
        }
        rcSteps.append(NEW_LINE);
        System.out.println(rcSteps.toString());
        return s.getBondChangeCalculator();
    }

    /**
     *
     * @param reactionID
     * @param directory
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public BondChangeCalculator map(String reactionID, String directory) throws FileNotFoundException, Exception {
        String NEW_LINE = System.getProperty("line.separator");
        IReaction cdkReaction = readReaction(reactionID, directory, false);
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, true, true, false);
        MappingSolution s = rmt.getSelectedSolution();
        new ImageGenerator().drawLeftToRightReactionLayout("Output", s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens(), (reactionID + s.getAlgorithmID()));
        return s.getBondChangeCalculator();
    }
}
