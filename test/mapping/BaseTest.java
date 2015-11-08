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

import java.awt.Image;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;
import uk.ac.ebi.reactionblast.tools.ImageGenerator;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BaseTest extends TestUtility {


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
