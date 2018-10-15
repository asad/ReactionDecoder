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

import java.awt.Image;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.System.out;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Reader;
import static java.lang.System.getProperty;
import static java.util.logging.Level.SEVERE;

import static javax.imageio.ImageIO.write;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.tools.ImageGenerator.LeftToRightReactionCenterImage;
import static uk.ac.ebi.reactionblast.tools.ImageGenerator.TopToBottomReactionLayoutImage;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MappingUtility extends TestUtility {

    private final static ILoggingTool LOGGER
            = createLoggingTool(MappingUtility.class);

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
        return readReactionFile(name, "Output/", reMap, false);
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
        return readReactionFile(name, dir, reMap, false);
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
        return readReactionFile(name, dir, false, false);
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
    public IReaction map(IReaction reaction) {
        try {
            ReactionMechanismTool rmt = new ReactionMechanismTool(
                    //                    reaction, true, true, true, new CDKReactionStandardizer());
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
     * @param reactionID
     * @param directory
     * @return
     * @throws FileNotFoundException
     * @throws Exception
     */
    public ReactionMechanismTool testReactions(String reactionID, String directory) throws FileNotFoundException, Exception {
        System.out.println("Mapping Reaction " + reactionID);
        IReaction cdkReaction = readReaction(reactionID, directory, false);
        ExtReactionManipulatorTool.addExplicitH(cdkReaction);
//        SmilesGenerator withAtomClasses = unique().aromatic().withAtomClasses();
//        out.println("Input reactions " + withAtomClasses.createReactionSMILES(cdkReaction));
        ReactionMechanismTool annotation = getAnnotation(cdkReaction);
        MappingSolution s = annotation.getSelectedSolution();
//        out.println("Mapped reactions " + withAtomClasses.createReactionSMILES(s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens()));
        return annotation;

    }

    /**
     *
     * @param cdkReaction
     * @return
     * @throws AssertionError
     * @throws Exception
     */
    public ReactionMechanismTool getAnnotation(IReaction cdkReaction) throws AssertionError, Exception {
        /*
         RMT for the reaction mapping
         */
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, true, true, false, new StandardizeReaction());
        MappingSolution s = rmt.getSelectedSolution();

//        out.println("Reaction ID: " + s.getReaction().getID() + ", Selected Algorithm: " + s.getAlgorithmID());
//        System.out.println("Cleaved/Formed " + s.getBondChangeCalculator().getFormedCleavedWFingerprint().toString());
//        System.out.println("Order Changed " + s.getBondChangeCalculator().getOrderChangesWFingerprint().toString());
//        System.out.println("Stereo Changed " + s.getBondChangeCalculator().getStereoChangesWFingerprint().toString());
//        System.out.println("RC Changed " + s.getBondChangeCalculator().getReactionCenterWFingerprint().toString());
//        System.out.println("BE " + s.getBondEnergySum() + ", Fragment " + s.getTotalFragmentChanges());
        try {
            IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();

            /*
             * Code for Image generation
             */
            LeftToRightReactionCenterImage(reactionWithCompressUnChangedHydrogens, (s.getReaction().getID() + s.getAlgorithmID() + "RC"), "Output");
            TopToBottomReactionLayoutImage(reactionWithCompressUnChangedHydrogens, (s.getReaction().getID() + s.getAlgorithmID()), "Output");
        } catch (Exception e) {
            LOGGER.error(SEVERE, null, e);
        }
//        int i = 1;
//        for (MappingSolution m : rmt.getAllSolutions()) {
//            out.println("--------------------------------------");
//            BondChangeCalculator bcc = m.getBondChangeCalculator();
//            out.println(m.getAlgorithmID() + ", fp " + bcc.getFormedCleavedWFingerprint().toString());
//            out.println(m.getAlgorithmID() + ", fp " + bcc.getOrderChangesWFingerprint().toString());
//
//            out.println("BE " + m.getBondEnergySum() + ", Fragment " + m.getTotalFragmentChanges());
//            new ImageGenerator().drawLeftToRightReactionLayout("Output", bcc.getReactionWithCompressUnChangedHydrogens(), ("Map_" + s.getReaction().getID() + m.getAlgorithmID()));
//            i++;
//            out.println();
//            out.println("--------------------------------------");
//        }
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
        String NEW_LINE = getProperty("line.separator");
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
        out.println(sb.toString());

        StringBuilder rcSteps = new StringBuilder();
        rcSteps.append("Formed Cleaved");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterFormedCleavedFingerprint = s.getBondChangeCalculator().getReactionCenterFormedCleavedFingerprint();
        reactionCenterFormedCleavedFingerprint.entrySet().stream().filter((m) -> !(m.getKey() == -1)).forEach((m) -> {
            rcSteps.append(m.getValue());
        });
        rcSteps.append("Order Change");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterOrderChangeFingerprint = s.getBondChangeCalculator().getReactionCenterOrderChangeFingerprint();
        reactionCenterOrderChangeFingerprint.entrySet().stream().filter((m) -> !(m.getKey() == -1)).forEach((m) -> {
            rcSteps.append(m.getValue());
        });
        rcSteps.append("Stereo Change");
        rcSteps.append(NEW_LINE);
        Map<Integer, IPatternFingerprinter> reactionCenterStereoChangeFingerprint = s.getBondChangeCalculator().getReactionCenterStereoChangeFingerprint();
        reactionCenterStereoChangeFingerprint.entrySet().stream().filter((m) -> !(m.getKey() == -1)).forEach((m) -> {
            rcSteps.append(m.getValue());
        });
        rcSteps.append(NEW_LINE);
        out.println(rcSteps.toString());
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
        IReaction cdkReaction = readReaction(reactionID, directory, false);
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, true, true, false);
        MappingSolution s = rmt.getSelectedSolution();
        new ImageGenerator().drawLeftToRightReactionLayout("Output", s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens(), (reactionID + s.getAlgorithmID()));
        return s.getBondChangeCalculator();
    }
}
