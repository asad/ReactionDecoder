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
import java.io.IOException;
import javax.imageio.ImageIO;
import javax.vecmath.Vector2d;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.graphics.direct.DirectMoleculeDrawer;
import uk.ac.ebi.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import uk.ac.ebi.reactionblast.graphics.direct.layout.ZoomToFitLayout;
import uk.ac.ebi.reactionblast.mapping.blocks.BlockReactionCanoniser;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalReactionLabeller;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class TestUtility {

    static final String KEGG_RXN_DIR = "rxn/kegg/";
    static final String Rhea_RXN_DIR = "rxn/rhea/";
    static final String Brenda_RXN_DIR = "rxn/brenda/";
    static final String BUG_RXN_DIR = "rxn/bug/";
    static final String OTHER_RXN = "rxn/other/";
    static final String METRXN_RXN = "rxn/metrxn/";
    static final String INFORCHEM_RXN = "rxn/infochem/";
    static final String MACIE_RXN = "rxn/macie/";

    protected void setLonePairs(IReaction reaction) {
        LonePairElectronChecker checker = new LonePairElectronChecker();
        setLonePairs(reaction.getReactants(), checker);
        setLonePairs(reaction.getProducts(), checker);
    }

    protected void setLonePairs(IAtomContainerSet AtomContainerSet, LonePairElectronChecker checker) {
        for (IAtomContainer atomContainer : AtomContainerSet.atomContainers()) {
            try {
                checker.saturate(atomContainer);
            } catch (CDKException c) {
                c.printStackTrace();
            }
        }
    }

    protected void detectAromaticity(IReaction reaction) {
        detectAromaticity(reaction.getReactants());
        detectAromaticity(reaction.getProducts());
    }

    protected void detectAromaticity(IAtomContainerSet molSet) {
        for (IAtomContainer ac : molSet.atomContainers()) {
            try {
                ExtAtomContainerManipulator.aromatizeCDK(ac);
            } catch (CDKException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }

    protected void setMappingIDs(IReaction reaction) {
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

    protected void addImplicitHydrogens(IReaction reaction) {
        addImplicitHydrogens(reaction.getReactants());
        addImplicitHydrogens(reaction.getProducts());
    }

    protected void addImplicitHydrogens(IAtomContainerSet molSet) {
        for (IAtomContainer atomContainer : molSet.atomContainers()) {
            addImplicitHydrogens(atomContainer);
        }
    }

    protected void addImplicitHydrogens(IAtomContainer atomContainer) {
        try {
            CDKHydrogenAdder.getInstance(
                    DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    protected void typeAtoms(IAtomContainer atomContainer) {
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
        } catch (CDKException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /**
     *
     * @param reaction
     * @return
     */
    protected IReaction canonise(IReaction reaction) {
        ICanonicalReactionLabeller labeller = new BlockReactionCanoniser();
        return labeller.getCanonicalReaction(reaction);
    }

    protected IAtomContainer layout(IAtomContainer AtomContainer) {
        try {
            StructureDiagramGenerator sdg = new StructureDiagramGenerator();
            sdg.setMolecule(AtomContainer, true);
            sdg.generateCoordinates();
            return sdg.getMolecule();
        } catch (CDKException e) {
            return AtomContainer;
        }
    }

    protected BufferedImage makeBlankImage(int width, int height) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_BGR);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);
        return image;
    }

    protected void draw(DirectMoleculeDrawer molDrawer, IAtomContainer AtomContainer,
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
