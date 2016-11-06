package uk.ac.ebi.reactionblast;

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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import javax.vecmath.Vector2d;
import org.openscience.cdk.AtomContainerSet;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.TITLE;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import static org.openscience.cdk.tools.CDKHydrogenAdder.getInstance;
import org.openscience.cdk.tools.LonePairElectronChecker;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import static uk.ac.ebi.reactionblast.mapping.helper.MappingHandler.cleanMapping;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeCDK;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.convertExplicitToImplicitHydrogens;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import java.net.URL;
import uk.ac.ebi.reactionblast.mechanism.helper.ReactionMappingUtility;
import static java.lang.String.valueOf;
import java.util.logging.Level;
import static java.util.logging.Logger.getLogger;
import static javax.imageio.ImageIO.write;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class TestUtility extends ReactionMappingUtility {

    public static final String KEGG_RXN_DIR = "rxn" + separator + "kegg" + separator;
    public static final String RHEA_RXN_DIR = "rxn" + separator + "rhea" + separator;
    public static final String BRENDA_RXN_DIR = "rxn" + separator + "brenda" + separator;
    public static final String BUG_RXN_DIR = "rxn" + separator + "bug" + separator;
    public static final String OTHER_RXN = "rxn" + separator + "other" + separator;
    public static final String METRXN_RXN = "rxn" + separator + "metrxn" + separator;
    public static final String INFORCHEM_RXN = "rxn" + separator + "infochem" + separator;
    public static final String MACIE_RXN = "rxn" + separator + "macie" + separator;
    public static final String USPTO_RXN = "rxn" + separator + "uspto" + separator;
    private static final Logger LOG = getLogger(TestUtility.class.getName());
    protected final URL resource;

    TestUtility() {
        resource = this.getClass().getClassLoader().getResource("");
    }

    /**
     *
     * @param reaction
     */
    protected static void setLonePairs(IReaction reaction) {
        LonePairElectronChecker checker = new LonePairElectronChecker();
        setLonePairs(reaction.getReactants(), checker);
        setLonePairs(reaction.getProducts(), checker);
    }

    /**
     *
     * @param AtomContainerSet
     * @param checker
     */
    protected static void setLonePairs(IAtomContainerSet AtomContainerSet, LonePairElectronChecker checker) {
        for (IAtomContainer atomContainer : AtomContainerSet.atomContainers()) {
            try {
                checker.saturate(atomContainer);
            } catch (CDKException c) {
                c.printStackTrace();
            }
        }
    }

    /**
     *
     * @param reaction
     */
    protected static void detectAromaticity(IReaction reaction) {
        detectAromaticity(reaction.getReactants());
        detectAromaticity(reaction.getProducts());
    }

    /**
     *
     * @param molSet
     */
    protected static void detectAromaticity(IAtomContainerSet molSet) {
        for (IAtomContainer ac : molSet.atomContainers()) {
            try {
                aromatizeCDK(ac);
            } catch (CDKException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }

    /**
     *
     * @param reaction
     */
    protected static void setMappingIDs(IReaction reaction) {
        int i = 0;
        for (IMapping mapping : reaction.mappings()) {
            IAtom a0 = (IAtom) mapping.getChemObject(0);
            IAtom a1 = (IAtom) mapping.getChemObject(1);
            Object mappingID = a0.getProperty(ATOM_ATOM_MAPPING);
            //            Integer mappingID =
            //                (Integer)a0.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
            if (mappingID != null) {
                a0.setID(valueOf(mappingID));
                a1.setID(valueOf(mappingID));
            } else {
                a0.setID(valueOf(i));
                a1.setID(valueOf(i));
            }
            i++;
        }
    }

    /**
     *
     * @param reaction
     */
    protected static void renumberMappingIDs(IReaction reaction) {
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
     * @param reaction
     */
    protected static void addImplicitHydrogens(IReaction reaction) {
        addImplicitHydrogens(reaction.getReactants());
        addImplicitHydrogens(reaction.getProducts());
    }

    /**
     *
     * @param molSet
     */
    protected static void addImplicitHydrogens(IAtomContainerSet molSet) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            addImplicitHydrogens(atomContainer);
        }
    }

    /**
     *
     * @param atomContainer
     */
    protected static void addImplicitHydrogens(IAtomContainer atomContainer) {
        try {
            getInstance(getInstance()).addImplicitHydrogens(atomContainer);
        } catch (CDKException ex) {
            Logger.getLogger(TestUtility.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @param atomContainer
     */
    protected static void typeAtoms(IAtomContainer atomContainer) {
        try {
            percieveAtomTypesAndConfigureAtoms(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /**
     *
     * @param AtomContainer
     * @return
     */
    protected static IAtomContainer layout(IAtomContainer AtomContainer) {
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
    protected static BufferedImage makeBlankImage(int width, int height) {
        BufferedImage image = new BufferedImage(width, height, TYPE_INT_BGR);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
        g.setColor(WHITE);
        g.fillRect(0, 0, width, height);
        return image;
    }

    /**
     *
     * @param molDrawer
     * @param AtomContainer
     * @param dirPath
     * @param name
     * @throws IOException
     */
    protected static void draw(DirectMoleculeDrawer molDrawer, IAtomContainer AtomContainer,
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
    protected static IReaction readReactionFile(String name, String dir, boolean reMap, boolean removeHydrogens) throws Exception {
        String filepath = dir + name + ".rxn";
        FileInputStream fileInputStream = new FileInputStream(new File(filepath));
        return readReactionFile(name, fileInputStream, reMap, removeHydrogens);
    }

    /**
     *
     * @param name
     * @param in
     * @param reMap
     * @param removeHydrogens
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    protected static IReaction readReactionFile(String name, InputStream in, boolean reMap, boolean removeHydrogens) throws Exception {
        IReaction reaction = null;
        try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(in)) {
            reaction = reader.read(new Reaction());
            reaction.setID(name);
        } catch (Exception ex) {
            getLogger(TestUtility.class.getName()).log(SEVERE, null, ex);
        }

        if (removeHydrogens && reaction != null) {
            // XXX WARNING : this may not work correctly!
            IReaction hydrogenFreeReaction = new Reaction();
            hydrogenFreeReaction.setID(name);
            IAtomContainerSet hydrogenFreeReactants = new AtomContainerSet();

            for (IAtomContainer atomContainer : reaction.getReactants().atomContainers()) {
                setNullHCountToZero(atomContainer);
                percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID((String) atomContainer.getProperty(TITLE));
                hydrogenFreeReactants.addAtomContainer(acMinusH);
            }
            hydrogenFreeReaction.setReactants(hydrogenFreeReactants);
            IAtomContainerSet hydrogenFreeProducts = new AtomContainerSet();
            for (IAtomContainer atomContainer : reaction.getProducts().atomContainers()) {
                setNullHCountToZero(atomContainer);
                percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID((String) atomContainer.getProperty(TITLE));
                hydrogenFreeProducts.addAtomContainer(acMinusH);
            }

            hydrogenFreeReaction.setProducts(hydrogenFreeProducts);
            for (IMapping mapping : reaction.mappings()) {
                if (((IElement) mapping.getChemObject(0)).getSymbol().equals("H")
                        || ((IElement) mapping.getChemObject(1)).getSymbol().equals("H")) {
                    continue;
                }
                hydrogenFreeReaction.addMapping(mapping);
            }
            reaction = hydrogenFreeReaction;
        }

        if (reMap) {
            cleanMapping(reaction);
        } else {
            renumberMappingIDs(reaction);
        }

        return reaction;
    }

    /**
     * Set all null hydrogen counts to 0. Generally hydrogen counts are present
     * and if not we add them. However the molecule being tested can't include
     * hydrogen counts as then fingerprints don't line up (substructure
     * filtering). The previous behaviour of the SMARTS matching was to treat
     * null hydrogens as 0 - the new behaviour is to complain about it.
     *
     * @param mol molecule to zero out hydrogen counts
     */
    static void setNullHCountToZero(IAtomContainer mol) {
        for (IAtom a : mol.atoms()) {
            if (a.getImplicitHydrogenCount() == null) {
                a.setImplicitHydrogenCount(0);
            }
        }
    }

    /**
     *
     * @param reactionSet
     * @param remap override existing mappings
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    protected List<IReaction> mapReactions(IReactionSet reactionSet, boolean remap) throws FileNotFoundException, Exception {
        List<IReaction> mappedReactionList = new ArrayList<>();
        for (IReaction cdkReaction : reactionSet.reactions()) {

            IReaction mappedReaction = mapReaction(cdkReaction, remap);
            /*
            Add mapped reaction to the list
             */ mappedReactionList.add(mappedReaction);
        }
        return mappedReactionList;
    }

    /**
     *
     * @param cdkReaction reaction for be mapped
     * @param remap override existing mappings
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    protected IReaction mapReaction(IReaction cdkReaction, boolean remap) throws FileNotFoundException, Exception {

        String reactionName = cdkReaction.getID();
        IReaction cleanReaction = cleanReaction(cdkReaction, reactionName);
        /*
        * RMT for the reaction mapping
         */
        ReactionMechanismTool rmt = new ReactionMechanismTool(cleanReaction, remap, new StandardizeReaction());

        /*
        Reaction with hydrogens mapped but unchanged hydrogens suppressed
         */
        //IReaction reactionWithCompressUnChangedHydrogens = rmt.getSelectedSolution().getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
        /*
        Reaction with hydrogens mapped
         */
        IReaction mappedReaction = rmt.getSelectedSolution().getReaction();

        /*
        optional step: Renumber the atoms as per mapping
         */
        renumberMappingIDs(mappedReaction);

        return mappedReaction;
    }

    /**
     *
     * @param reaction
     * @param reactionName
     * @return
     * @throws FileNotFoundException
     */
    protected IReaction cleanReaction(IReaction reaction, String reactionName) throws Exception {
        //write code to fix reactions (Atom type , hydrogens etc.)
        //TO DO
        reaction.setID(reactionName);
        return reaction;
    }
}
