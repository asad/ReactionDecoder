/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.aamtool;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getAtomCount;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.ReactionFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.mechanism.helper.MoleculeMoleculePair;
import static uk.ac.ebi.reactionblast.tools.ReactionSimilarityTool.getSimilarity;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Annotator extends Helper {

    static final String NEW_LINE = getProperty("line.separator");
    static final String TAB = "\t";
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Annotator.class);

    /**
     *
     */
    protected boolean REPORT_ALL_MAPPINGS;

    /**
     *
     */
    protected boolean GENERATE_IMAGE;

    /**
     *
     */
    protected boolean GENERATE_AAMIMAGE;

    /**
     *
     */
    protected boolean REPORT_MMP;

    /**
     *
     */
    protected boolean REPORT_PATTERNS;

    /**
     *
     */
    protected boolean REMAP;

    /**
     *
     */
    protected String PREFIX;

    Annotator() {
        this.REPORT_ALL_MAPPINGS = false;
        this.GENERATE_IMAGE = false;
        this.GENERATE_AAMIMAGE = false;
        this.REPORT_MMP = false;
        this.REPORT_PATTERNS = false;
        this.REMAP = true;
        this.PREFIX = "";
    }

    /**
     *
     * @param cdkReaction
     * @param reMap remap the reaction
     * @return
     * @throws Exception
     */
    protected ReactionMechanismTool getReactionMechanismTool(IReaction cdkReaction, boolean reMap) throws Exception {
        ReactionMechanismTool rmt;
        /*
         Check if the reaction is already mapped
         */
        if (getAtomCount(cdkReaction.getReactants()) == cdkReaction.getMappingCount()) {
            cdkReaction.setFlag(MAPPED, true);
        } else {
            cdkReaction.setFlag(MAPPED, false);
        }
        rmt = new ReactionMechanismTool(cdkReaction, reMap, true, false, new StandardizeReaction());
//        IPatternFingerprinter formedCleavedWFingerprint = rmt
//                .getSelectedSolution()
//                .getBondChangeCalculator()
//                .getFormedCleavedWFingerprint();
//        System.out.println("formedCleavedWFingerprint " + formedCleavedWFingerprint);
        return rmt;
    }

    /**
     *
     * @param reactionID
     * @param mech
     * @return
     * @throws IOException
     * @throws CDKException
     * @throws Exception
     */
    protected boolean writeFiles(String reactionID, ReactionMechanismTool mech) throws IOException, CDKException, Exception {

        MappingSolution s = mech.getSelectedSolution();
        if (s == null) {
            return false;
        }
        File writeRXNMappedFile = writeRXNMappedFile(new File(".").getCanonicalPath(), s.getBondChangeCalculator().getReaction(), reactionID);
        out.println("Mapped RXN File " + writeRXNMappedFile.getAbsolutePath());

        if (GENERATE_IMAGE) {
            try {
                File generateImage = generateImage(new File(".").getCanonicalPath(), s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens(), reactionID);
                out.println("Annotated RXN Image " + generateImage.getAbsolutePath());
            } catch (Exception e) {
                LOGGER.error(SEVERE, "Unable to generate AAM image", e);
            }
        } else if (!GENERATE_IMAGE && GENERATE_AAMIMAGE) {
            try {
                File generateImage = generateAAMImage(new File(".").getCanonicalPath(), s.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens(), reactionID);
                out.println("Annotated RXN Image " + generateImage.getAbsolutePath());
            } catch (Exception e) {
                LOGGER.error(SEVERE, "Unable to generate AAM image", e);
            }
        }
        return true;
    }

    /**
     *
     * @param results
     * @param jobID
     * @throws IOException
     */
    private void writeSimilarityMatrix(List<SimilarityResult> results, String jobID) throws IOException {
        String rootPath = new File(".").getCanonicalPath();
        File bcMatrix = new File(rootPath, jobID + "_Bond_Change" + ".mat");
        File rcMatrix = new File(rootPath, jobID + "_Reaction_Centre" + ".mat");
        File stMatrix = new File(rootPath, jobID + "_Structure_Similarity" + ".mat");

        FileWriter writerBC = new FileWriter(bcMatrix);
        BufferedWriter bufferedWriterBC = new BufferedWriter(writerBC);

        FileWriter writerRC = new FileWriter(rcMatrix);
        BufferedWriter bufferedWriterRC = new BufferedWriter(writerRC);

        FileWriter writerST = new FileWriter(stMatrix);
        BufferedWriter bufferedWriterST = new BufferedWriter(writerST);

        bufferedWriterBC.newLine();
        try {
            for (SimilarityResult s : results) {
                if (s.getSimilarityReactions().containsKey("BC")) {
                    bufferedWriterBC.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + s.getSimilarityReactions().get("BC"));
                    bufferedWriterBC.newLine();
                } else {
                    bufferedWriterBC.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + "NA");
                    bufferedWriterBC.newLine();
                }
            }
        } finally {
            bufferedWriterBC.close();
        }

        try {
            for (SimilarityResult s : results) {
                if (s.getSimilarityReactions().containsKey("RC")) {
                    bufferedWriterRC.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + s.getSimilarityReactions().get("RC"));
                    bufferedWriterRC.newLine();
                } else {
                    bufferedWriterRC.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + "NA");
                    bufferedWriterRC.newLine();
                }
            }
        } finally {
            bufferedWriterRC.close();
        }

        try {
            for (SimilarityResult s : results) {
                if (s.getSimilarityReactions().containsKey("ST")) {
                    bufferedWriterST.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + s.getSimilarityReactions().get("ST"));
                    bufferedWriterST.newLine();
                } else {
                    bufferedWriterST.write("\"" + s.getQuery() + "\"" + TAB + "\"" + s.getTarget() + "\"" + TAB + "NA");
                    bufferedWriterST.newLine();
                }
            }
        } finally {
            bufferedWriterST.close();
        }

    }

    private void printRPAIRPatternAsText(MappingSolution s, StringBuilder sb) throws CloneNotSupportedException {
        Map<String, Collection<String>> moleculeMoleculeTransformationPairs = s.getBondChangeCalculator().getMoleculeMoleculeTransformationPairs();

        StringBuilder sbcomp = new StringBuilder();
        int index = 1;
        for (String m : moleculeMoleculeTransformationPairs.keySet()) {
            StringBuilder mmp = new StringBuilder(m);
            mmp.append("\t");
            mmp.append(moleculeMoleculeTransformationPairs.get(m));
            sbcomp.append(index).append(": ").append(mmp);
            sbcomp.append(NEW_LINE);
            index++;
        }

        Collection<MoleculeMoleculePair> reactionTransform = s.getBondChangeCalculator().getReactionCentreTransformationPairs();

        StringBuilder pair1 = new StringBuilder();
        index = 1;
        for (MoleculeMoleculePair m : reactionTransform) {
            pair1.append(index).append(": ").append(m.getSmirks1());
            pair1.append(NEW_LINE);
            index++;
        }

        StringBuilder pair2 = new StringBuilder();
        index = 1;
        for (MoleculeMoleculePair m : reactionTransform) {
            pair2.append(index).append(": ").append(m.getSmirks2());
            pair2.append(NEW_LINE);
            index++;
        }

        StringBuilder pair3 = new StringBuilder();
        index = 1;
        for (MoleculeMoleculePair m : reactionTransform) {
            pair3.append(index).append(": ").append(m.getSmirks3());
            pair3.append(NEW_LINE);
            index++;
        }

        StringBuilder sbFC = new StringBuilder();
        StringBuilder sbOC = new StringBuilder();
        StringBuilder sbST = new StringBuilder();

        Map<Integer, IPatternFingerprinter> reactionCenterFormedCleavedFingerprint = s.getBondChangeCalculator().getReactionCenterFormedCleavedFingerprint();

        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterFormedCleavedFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            IPatternFingerprinter value = m.getValue().clone();
            value.setFingerprintID("Reaction Center at Level: " + (m.getKey()));
            sbFC.append(value).append(NEW_LINE);
        }

        Map<Integer, IPatternFingerprinter> reactionCenterOrderChangeFingerprint = s.getBondChangeCalculator().getReactionCenterOrderChangeFingerprint();

        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterOrderChangeFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            IPatternFingerprinter value = m.getValue().clone();
            value.setFingerprintID("Reaction Center at Level: " + (m.getKey()));
            sbOC.append(value).append(NEW_LINE);
        }

        Map<Integer, IPatternFingerprinter> reactionCenterStereoChangeFingerprint = s.getBondChangeCalculator().getReactionCenterStereoChangeFingerprint();
        for (Map.Entry<Integer, IPatternFingerprinter> m : reactionCenterStereoChangeFingerprint.entrySet()) {
            if (m.getKey() == -1) {
                continue;
            }
            IPatternFingerprinter value = m.getValue().clone();
            value.setFingerprintID("Reaction Center at Level: " + (m.getKey()));
            sbST.append(value).append(NEW_LINE);
        }

        sb.append(NEW_LINE);
        sb.append("//");
        sb.append(NEW_LINE);
        sb.append("Reaction Centre Formed/Cleaved");
        sb.append(NEW_LINE);
        sb.append(sbFC.toString());
        sb.append(NEW_LINE);
        sb.append("Reaction Centre Order Changed");
        sb.append(NEW_LINE);
        sb.append(sbOC.toString());
        sb.append(NEW_LINE);
        sb.append("Reaction Centre Stereo Changes");
        sb.append(NEW_LINE);
        sb.append(sbST.toString());
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("//").append(NEW_LINE);
        sb.append("TRANSFORMATIONS");
        sb.append(NEW_LINE);
        sb.append("MMP Level 1");
        sb.append(NEW_LINE);
        sb.append(pair1.toString());
        sb.append(NEW_LINE);
        sb.append("MMP Level 2");
        sb.append(NEW_LINE);
        sb.append(pair2.toString());
        sb.append(NEW_LINE);
        sb.append("MMP Level 3");
        sb.append(NEW_LINE);
        sb.append(pair3.toString());
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("//").append(NEW_LINE);
        sb.append("REACTION MMP (RPAIR)");
        sb.append(NEW_LINE);
        sb.append(sbcomp.toString());
        sb.append(NEW_LINE);
    }

    private void printRPAIRPatternAsXML(MappingSolution s, org.w3c.dom.Document doc, org.w3c.dom.Element rootElement) {

        Map<Integer, IPatternFingerprinter> reactionCenterFormedCleavedFingerprint = s.getBondChangeCalculator().getReactionCenterFormedCleavedFingerprint();
        Map<Integer, IPatternFingerprinter> reactionCenterOrderChangeFingerprint = s.getBondChangeCalculator().getReactionCenterOrderChangeFingerprint();
        Map<Integer, IPatternFingerprinter> reactionCenterStereoChangeFingerprint = s.getBondChangeCalculator().getReactionCenterStereoChangeFingerprint();

        /*
         Get level Information
         */
        Set<Integer> levels = reactionCenterFormedCleavedFingerprint.keySet();

        levels.stream().filter((i) -> !(i == -1)).forEachOrdered((i) -> {
            //Start of Fingerprint elements
            org.w3c.dom.Element rc = doc.createElement("ReactionCenters");
            //Start of Fingerprint elements
            rootElement.appendChild(rc);
            //Start of BC as child node of Fingerprint elements
            org.w3c.dom.Attr attr = doc.createAttribute("LEVEL");
            attr.setValue(i + "");
            rc.setAttributeNode(attr);
            if (reactionCenterFormedCleavedFingerprint.containsKey(i)) {
                // FC elements
                org.w3c.dom.Element fp_FORMED_CLEAVED = doc.createElement("FC");
                fp_FORMED_CLEAVED.appendChild(doc.createTextNode(reactionCenterFormedCleavedFingerprint.get(i).getFeatures().toString()));
                rc.appendChild(fp_FORMED_CLEAVED);
            }
            if (reactionCenterOrderChangeFingerprint.containsKey(i)) {
                // OC elements
                org.w3c.dom.Element fp_ORDER_CHANGED = doc.createElement("OC");
                fp_ORDER_CHANGED.appendChild(doc.createTextNode(reactionCenterOrderChangeFingerprint.get(i).getFeatures().toString()));
                rc.appendChild(fp_ORDER_CHANGED);
            }
            if (reactionCenterStereoChangeFingerprint.containsKey(i)) {
                // ST elements
                org.w3c.dom.Element fp_STEREO_CHANGED = doc.createElement("ST");
                fp_STEREO_CHANGED.appendChild(doc.createTextNode(reactionCenterStereoChangeFingerprint.get(i).getFeatures().toString()));
                rc.appendChild(fp_STEREO_CHANGED);
            }
        });

        Collection<MoleculeMoleculePair> reactionTransform = s.getBondChangeCalculator().getReactionCentreTransformationPairs();

        //Start of Fingerprint elements
        org.w3c.dom.Element transform = doc.createElement("TRANSFORMATION");
        rootElement.appendChild(transform);

        //Start of RPAIR as child node of Fingerprint elements
        org.w3c.dom.Attr attr = doc.createAttribute("LEVEL");
        attr.setValue(1 + "");
        transform.setAttributeNode(attr);

        int index = 1;

        for (MoleculeMoleculePair m : reactionTransform) {
            // RAIR elements
            org.w3c.dom.Element rpairMATCH = doc.createElement("MMP" + index);
            rpairMATCH.appendChild(doc.createTextNode(m.getSmirks1()));
            transform.appendChild(rpairMATCH);
            index++;
        }

        //Start of Fingerprint elements
        transform = doc.createElement("TRANSFORMATION");
        rootElement.appendChild(transform);

        //Start of RPAIR as child node of Fingerprint elements
        attr = doc.createAttribute("LEVEL");
        attr.setValue(2 + "");
        transform.setAttributeNode(attr);
        index = 1;

        for (MoleculeMoleculePair m : reactionTransform) {
            // RAIR elements
            org.w3c.dom.Element rpairMATCH = doc.createElement("MMP" + index);
            rpairMATCH.appendChild(doc.createTextNode(m.getSmirks2()));
            transform.appendChild(rpairMATCH);
            index++;
        }

        //Start of Fingerprint elements
        transform = doc.createElement("TRANSFORMATION");
        rootElement.appendChild(transform);

        //Start of RPAIR as child node of Fingerprint elements
        attr = doc.createAttribute("LEVEL");
        attr.setValue(3 + "");
        transform.setAttributeNode(attr);
        index = 1;

        for (MoleculeMoleculePair m : reactionTransform) {
            // RAIR elements
            org.w3c.dom.Element rpairMATCH = doc.createElement("MMP" + index);
            rpairMATCH.appendChild(doc.createTextNode(m.getSmirks3()));
            transform.appendChild(rpairMATCH);
            index++;
        }

        Map<String, Collection<String>> moleculeMoleculeTransformationPairs = s.getBondChangeCalculator().getMoleculeMoleculeTransformationPairs();

        index = 1;
        for (String m : moleculeMoleculeTransformationPairs.keySet()) {

            //Start of Fingerprint elements
            org.w3c.dom.Element rpair = doc.createElement("RPAIR");
            rootElement.appendChild(rpair);

            //Start of RPAIR as child node of Fingerprint elements
            attr = doc.createAttribute("COUNT");
            attr.setValue(index + "");
            rpair.setAttributeNode(attr);

            // RAIR elements
            org.w3c.dom.Element rpairMATCH = doc.createElement("MMP");
            Collection<String> mmp = moleculeMoleculeTransformationPairs.get(m);
            StringBuilder sb = new StringBuilder(m);
            sb.append("\t");
            sb.append(mmp);
            rpairMATCH.appendChild(doc.createTextNode(sb.toString()));
            rpair.appendChild(rpairMATCH);
            index++;
        }
    }

    /**
     *
     * @param rmt
     * @param reactionID
     * @param sb
     * @throws java.lang.CloneNotSupportedException
     */
    protected void annotateReactionAsText(ReactionMechanismTool rmt, String reactionID, StringBuilder sb) throws CloneNotSupportedException {
        DecimalFormatSymbols instance = DecimalFormatSymbols.getInstance();
        instance.setExponentSeparator("E");//x10^
        DecimalFormat df = new DecimalFormat("##E00", instance);
        NumberFormat myFormatter = NumberFormat.getInstance();
        myFormatter.setMinimumFractionDigits(2);
        myFormatter.setMaximumFractionDigits(2);
        try {
            MappingSolution s = rmt.getSelectedSolution();
            if (s == null) {
                out.println("No valid solution found");
                return;
            }

            if (REPORT_PATTERNS) {
                sb.append(NEW_LINE);
                sb.append("//");
                sb.append(NEW_LINE);
                //Start of Fingerprint elements
                sb.append("FINGERPRINTS BC");
                sb.append(NEW_LINE);
                if (!s.getBondChangeCalculator().getFormedCleavedWFingerprint()
                        .getFeatures().isEmpty()) {
                    sb.append(NEW_LINE);
                    sb.append("FORMED_CLEAVED");
                    sb.append(NEW_LINE);
                    sb.append(s.getBondChangeCalculator().getFormedCleavedWFingerprint()
                            .getFeatures().toString()).append(NEW_LINE);
                }
                if (!s.getBondChangeCalculator().getOrderChangesWFingerprint()
                        .getFeatures().isEmpty()) {
                    sb.append(NEW_LINE);
                    sb.append("ORDER_CHANGED");
                    sb.append(NEW_LINE);
                    sb.append(s.getBondChangeCalculator().getOrderChangesWFingerprint()
                            .getFeatures().toString()).append(NEW_LINE);
                }
                //
                if (!s.getBondChangeCalculator().getStereoChangesWFingerprint()
                        .getFeatures().isEmpty()) {
                    // fp_STEREO_CHANGED elements
                    sb.append(NEW_LINE);
                    sb.append("STEREO_CHANGED");
                    sb.append(NEW_LINE);
                    sb.append(s.getBondChangeCalculator().getStereoChangesWFingerprint()
                            .getFeatures().toString()).append(NEW_LINE);
                }
                //Start of Fingerprint elements
                sb.append(NEW_LINE);
                sb.append("//");
                sb.append(NEW_LINE);
                sb.append("FINGERPRINTS RC");
                sb.append(NEW_LINE);
                sb.append(s.getBondChangeCalculator().getReactionCenterWFingerprint()
                        .getFeatures().toString()).append(NEW_LINE);
                /*
                 Call RPAIR type Transformations
                 */

                if (REPORT_MMP) {
                    printRPAIRPatternAsText(s, sb);
                }
            }

            /*
             * Selected AAM solution
             */
            SmilesGenerator smileGenerator = new SmilesGenerator(
                    SmiFlavor.Unique
                    | SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo
            );
            //Start of Fingerprint elements
            sb.append(NEW_LINE);
            sb.append("//");
            sb.append(NEW_LINE);
            sb.append("SELECTED AAM MAPPING");
            sb.append(NEW_LINE);
            //Start of Fingerprint elements
            try {
                IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator()
                        .getReactionWithCompressUnChangedHydrogens();
                sb.append(smileGenerator.create(reactionWithCompressUnChangedHydrogens));
            } catch (CDKException e) {
                LOGGER.info("Error in creating reaction SMILES ");
                LOGGER.error(SEVERE, null, e);
            }
            sb.append(NEW_LINE);
            sb.append(NEW_LINE);

            /*
             * Old atom index mapped to new mapping Index
             */
            if (REMAP) {
                //Start of Mapping Information
                sb.append(NEW_LINE);
                sb.append("//");
                sb.append(NEW_LINE);
                sb.append("REACTANT INPUT ATOM INDEX<->AAM ID");
                sb.append(NEW_LINE);
                //Start of Fingerprint elements
                sb.append(s.getReactor().getInputRankLabelledAtomsReactant());
                sb.append(NEW_LINE);
                sb.append("PRODUCT INPUT ATOM INDEX<->AAM ID");
                sb.append(NEW_LINE);
                //Start of Mapping Information
                sb.append(s.getReactor().getInputRankLabelledAtomsProduct());
                sb.append(NEW_LINE);
                sb.append(NEW_LINE);
            }
            if (REPORT_ALL_MAPPINGS) {
                int index = 1;
                for (MappingSolution m : rmt.getAllSolutions()) {
                    //Start of Fingerprint elements
                    sb.append("//");
                    sb.append(NEW_LINE);
                    sb.append(index).append(") AAM MAPPING ALGORITHM: ").append(m.getAlgorithmID().description());
                    sb.append(NEW_LINE);
                    //Start of Fingerprint elements
                    sb.append(smileGenerator.create(m.getBondChangeCalculator()
                            .getReactionWithCompressUnChangedHydrogens()));
                    sb.append(NEW_LINE);
                    //Start of Fingerprint elements
                    sb.append("SCORE: ").append((m.getTotalBondChanges() + m.getTotalFragmentChanges()));
                    sb.append(", CHAOS: ").append(m.getTotalFragmentChanges()).append(" <=> ").append(m.getSmallestFragmentCount());
                    sb.append(", SIGMA: ").append(m.getTotalCarbonBondChanges());
                    sb.append(", ENERGY: ").append(df.format(m.getBondEnergySum()));
                    sb.append(", DELTA: ").append(df.format(m.getEnergyDelta()));
                    sb.append(NEW_LINE);
                    index++;
                }
            }
        } catch (CDKException ex) {
            LOGGER.debug("Invalid RXN File " + reactionID);
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     *
     * @param rmt
     * @param reactionID
     * @param doc
     * @param rootElement
     */
    protected void annotateReactionAsXML(ReactionMechanismTool rmt, String reactionID, Document doc, Element rootElement) {
        DecimalFormatSymbols instance = DecimalFormatSymbols.getInstance();
        instance.setExponentSeparator("E");//x10^
        DecimalFormat df = new DecimalFormat("##E00", instance);

        NumberFormat myFormatter = NumberFormat.getInstance();
        myFormatter.setMinimumFractionDigits(2);
        myFormatter.setMaximumFractionDigits(2);
        Element annot = doc.createElement("ANNOTATION");
        rootElement.appendChild(annot);
        try {
            MappingSolution s = rmt.getSelectedSolution();
            if (s == null) {
                out.println("No valid solution found");
                return;
            }

            if (REPORT_PATTERNS) {
                //Start of Fingerprint elements
                Element fp = doc.createElement("FINGERPRINTS");
                annot.appendChild(fp);
                //Start of BC as child node of Fingerprint elements
                Attr attr = doc.createAttribute("BC");
                attr.setValue("1");
                fp.setAttributeNode(attr);
                if (!s.getBondChangeCalculator().getFormedCleavedWFingerprint().getFeatures().isEmpty()) {
                    // fp_Reaction_Centre elements
                    Element fp_Formed_Cleaved = doc.createElement("FORMED_CLEAVED");
                    fp_Formed_Cleaved.appendChild(doc.createTextNode(s.getBondChangeCalculator()
                            .getFormedCleavedWFingerprint().getFeatures().toString()));
                    fp.appendChild(fp_Formed_Cleaved);
                }
                if (!s.getBondChangeCalculator().getOrderChangesWFingerprint().getFeatures().isEmpty()) {
                    // fp_STEREO_CHANGED elements
                    Element fp_ORDER_CHANGED = doc.createElement("ORDER_CHANGED");
                    fp_ORDER_CHANGED.appendChild(doc.createTextNode(s.getBondChangeCalculator()
                            .getOrderChangesWFingerprint().getFeatures().toString()));
                    fp.appendChild(fp_ORDER_CHANGED);
                }
                //
                if (!s.getBondChangeCalculator().getStereoChangesWFingerprint().getFeatures().isEmpty()) {
                    // fp_STEREO_CHANGED elements
                    Element fp_STEREO_CHANGED = doc.createElement("STEREO_CHANGED");
                    fp_STEREO_CHANGED.appendChild(doc.createTextNode(s.getBondChangeCalculator()
                            .getStereoChangesWFingerprint().getFeatures().toString()));
                    fp.appendChild(fp_STEREO_CHANGED);
                }
                //Start of Fingerprint elements
                fp = doc.createElement("FINGERPRINTS");
                annot.appendChild(fp);
                //Start of RC as child node of Fingerprint elements
                attr = doc.createAttribute("RC");
                attr.setValue("2");
                fp.setAttributeNode(attr);
                /*
                 fp_Reaction_Centre elements
                 */
                Element fp_Reaction_Centre = doc.createElement("CENTRE");
                fp_Reaction_Centre.appendChild(doc.createTextNode(s.getBondChangeCalculator()
                        .getReactionCenterWFingerprint().getFeatures().toString()));
                fp.appendChild(fp_Reaction_Centre);

                /*
                 Call RPAIR type Transformations
                 */
                if (REPORT_MMP) {
                    printRPAIRPatternAsXML(s, doc, annot);
                }
            }

            /*
             * Selected AAM solution
             */
            SmilesGenerator smileGenerator = new SmilesGenerator(
                    SmiFlavor.Unique
                    | SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo
            );
            //Start of Fingerprint elements
            Element aam = doc.createElement("MAPPING");
            annot.appendChild(aam);

            //Start of BEST SOL as child node of AAM elements
            Attr attr = doc.createAttribute("STATUS");
            attr.setValue("SELECTED");
            aam.setAttributeNode(attr);
            // AAM elements
            Element selected_AAM = doc.createElement("AAM");

            //Start of Fingerprint elements
            try {
                IReaction reactionWithCompressUnChangedHydrogens = s.getBondChangeCalculator()
                        .getReactionWithCompressUnChangedHydrogens();
                selected_AAM.appendChild(doc.createTextNode(smileGenerator.create(reactionWithCompressUnChangedHydrogens)));
            } catch (CDKException e) {
                LOGGER.info("Error in creating reaction SMILES ");
                LOGGER.error(SEVERE, null, e);
            }
            aam.appendChild(selected_AAM);

            //OLD RANK
            /*
             * Old atom index mapped to new mapping Index
             */
            if (REMAP) {
                String reactant_atom_rank = s.getReactor().getInputRankLabelledAtomsProduct().toString();
                String product_atom_rank = s.getReactor().getInputRankLabelledAtomsProduct().toString();

                Element selected_AAM_RANK_R = doc.createElement("RANK_REACTANT");
                selected_AAM_RANK_R.appendChild(doc.createTextNode(reactant_atom_rank));
                aam.appendChild(selected_AAM_RANK_R);

                Element selected_AAM_RANK_P = doc.createElement("RANK_PRODUCT");
                selected_AAM_RANK_P.appendChild(doc.createTextNode(product_atom_rank));
                aam.appendChild(selected_AAM_RANK_P);
            }

            if (REPORT_ALL_MAPPINGS) {
                for (MappingSolution m : rmt.getAllSolutions()) {
                    //Start of Fingerprint elements
                    aam = doc.createElement("MAPPING");
                    annot.appendChild(aam);
                    //Start of BEST SOL as child node of AAM elements
                    attr = doc.createAttribute("ALGORTIHM");
                    attr.setValue(m.getAlgorithmID().description());
                    aam.setAttributeNode(attr);
                    // AAM elements
                    Element solutionAAM = doc.createElement("AAM");
                    solutionAAM.appendChild(doc.createTextNode(smileGenerator.create(m.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens())));
                    aam.appendChild(solutionAAM);
                    // AAM elements
                    Element score = doc.createElement("SCORE");
                    score.appendChild(doc.createTextNode((m.getTotalBondChanges() + m.getTotalFragmentChanges()) + ""));
                    aam.appendChild(score);
                    // AAM elements
                    score = doc.createElement("CHAOS");
                    score.appendChild(doc.createTextNode(m.getTotalFragmentChanges() + " <=> " + m.getSmallestFragmentCount()));
                    aam.appendChild(score);
                    // AAM elements
                    score = doc.createElement("SIGMA");
                    score.appendChild(doc.createTextNode(m.getTotalCarbonBondChanges() + ""));
                    aam.appendChild(score);
                    // AAM elements
                    score = doc.createElement("ENERGY");
                    score.appendChild(doc.createTextNode(df.format(m.getBondEnergySum()) + ""));
                    aam.appendChild(score);
                    // AAM elements
                    score = doc.createElement("DELTA");
                    score.appendChild(doc.createTextNode(df.format(m.getEnergyDelta()) + ""));
                    aam.appendChild(score);
                }
            }
        } catch (CDKException ex) {
            LOGGER.debug("Invalid RXN File " + reactionID);
            LOGGER.error(SEVERE, null, ex);
        }
    }

    /**
     *
     * @param annotateRXNQ
     * @param reactionQID
     * @param annotateRXNT
     * @param reactionTID
     * @param doc
     * @param rootElement
     * @throws Exception
     */
    protected void compareRXNXML(ReactionMechanismTool annotateRXNQ, String reactionQID, ReactionMechanismTool annotateRXNT, String reactionTID, Document doc, Element rootElement) throws Exception {
        NumberFormat myFormatter = NumberFormat.getInstance();
        myFormatter.setMinimumFractionDigits(2);
        myFormatter.setMaximumFractionDigits(2);
        Element element = doc.createElement("COMPARISON");
        rootElement.appendChild(element);
        //Start of Fingerprint elements
        Element query = doc.createElement("QUERY");
        element.appendChild(query);
        annotateReactionAsXML(annotateRXNQ, reactionQID, doc, query);
        BondChangeCalculator bondChangeCalculatorQ = annotateRXNQ.getSelectedSolution().getBondChangeCalculator();
        //Start of Fingerprint elements
        Element target = doc.createElement("TARGET");
        element.appendChild(target);
        annotateReactionAsXML(annotateRXNT, reactionTID, doc, target);
        BondChangeCalculator bondChangeCalculatorT = annotateRXNT.getSelectedSolution().getBondChangeCalculator();
        IPatternFingerprinter fpQ = new PatternFingerprinter();
        fpQ.add(bondChangeCalculatorQ.getFormedCleavedWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getOrderChangesWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getStereoChangesWFingerprint());
        IPatternFingerprinter fpT = new PatternFingerprinter();
        fpT.add(bondChangeCalculatorT.getFormedCleavedWFingerprint());
        fpT.add(bondChangeCalculatorT.getOrderChangesWFingerprint());
        fpT.add(bondChangeCalculatorT.getStereoChangesWFingerprint());
        double similarityBondChanges = getSimilarity(fpQ, fpT);
        //Start of Fingerprint elements
        Element sim = doc.createElement("SIMILARITY");
        element.appendChild(sim);
        //Start of RC as child node of Fingerprint elements
        Attr attr = doc.createAttribute("BC");
        attr.setValue("1");
        sim.setAttributeNode(attr);
        // AAM elements
        Element score = doc.createElement("SCORE");
        score.appendChild(doc.createTextNode(myFormatter.format(similarityBondChanges)));
        sim.appendChild(score);
        double similarityReactionCentres = getSimilarity(bondChangeCalculatorQ.getReactionCenterWFingerprint(), bondChangeCalculatorT.getReactionCenterWFingerprint());
        //Start of Fingerprint elements
        sim = doc.createElement("SIMILARITY");
        element.appendChild(sim);
        //Start of RC as child node of Fingerprint elements
        attr = doc.createAttribute("RC");
        attr.setValue("2");
        sim.setAttributeNode(attr);
        // AAM elements
        score = doc.createElement("SCORE");
        score.appendChild(doc.createTextNode(myFormatter.format(similarityReactionCentres)));
        sim.appendChild(score);
        ReactionFingerprinter rfQ = new ReactionFingerprinter(bondChangeCalculatorQ.getReaction());
        ReactionFingerprinter rfT = new ReactionFingerprinter(bondChangeCalculatorT.getReaction());
        double similarityReactionStructure = getSimilarity(rfQ.getReactionStruturalFingerprint(), rfT.getReactionStruturalFingerprint());
        //Start of Fingerprint elements
        sim = doc.createElement("SIMILARITY");
        element.appendChild(sim);
        //Start of RC as child node of Fingerprint elements
        attr = doc.createAttribute("ST");
        attr.setValue("3");
        sim.setAttributeNode(attr);
        // AAM elements
        score = doc.createElement("SCORE");
        score.appendChild(doc.createTextNode(myFormatter.format(similarityReactionStructure)));
        sim.appendChild(score);
    }

    /**
     *
     * @param annotateRXNQ
     * @param reactionQID
     * @param annotateRXNT
     * @param reactionTID
     * @return
     * @throws Exception
     */
    protected Map<String, String> similarityReactions(ReactionMechanismTool annotateRXNQ, String reactionQID, ReactionMechanismTool annotateRXNT, String reactionTID) throws Exception {

        Map<String, String> scores = new HashMap<>();

        NumberFormat myFormatter = NumberFormat.getInstance();
        myFormatter.setMinimumFractionDigits(2);
        myFormatter.setMaximumFractionDigits(2);
        BondChangeCalculator bondChangeCalculatorQ = annotateRXNQ.getSelectedSolution().getBondChangeCalculator();
        BondChangeCalculator bondChangeCalculatorT = annotateRXNT.getSelectedSolution().getBondChangeCalculator();
        IPatternFingerprinter fpQ = new PatternFingerprinter();
        fpQ.add(bondChangeCalculatorQ.getFormedCleavedWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getOrderChangesWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getStereoChangesWFingerprint());
        IPatternFingerprinter fpT = new PatternFingerprinter();
        fpT.add(bondChangeCalculatorT.getFormedCleavedWFingerprint());
        fpT.add(bondChangeCalculatorT.getOrderChangesWFingerprint());
        fpT.add(bondChangeCalculatorT.getStereoChangesWFingerprint());
        double similarityBondChanges = getSimilarity(fpQ, fpT);
        scores.put("BC", myFormatter.format(similarityBondChanges));
        double similarityReactionCentres = getSimilarity(bondChangeCalculatorQ.getReactionCenterWFingerprint(), bondChangeCalculatorT.getReactionCenterWFingerprint());
        scores.put("RC", myFormatter.format(similarityReactionCentres));
        ReactionFingerprinter rfQ = new ReactionFingerprinter(bondChangeCalculatorQ.getReaction());
        ReactionFingerprinter rfT = new ReactionFingerprinter(bondChangeCalculatorT.getReaction());
        double similarityReactionStructure = getSimilarity(rfQ.getReactionStruturalFingerprint(), rfT.getReactionStruturalFingerprint());
        scores.put("ST", myFormatter.format(similarityReactionStructure));

        return scores;
    }

    /**
     *
     * @param annotateRXNQ
     * @param reactionQID
     * @param annotateRXNT
     * @param reactionTID
     * @param sb StreactionWithLayouting buildereactionWithLayout
     * @throws Exception
     */
    protected void compareRXNText(ReactionMechanismTool annotateRXNQ, String reactionQID, ReactionMechanismTool annotateRXNT, String reactionTID, StringBuilder sb) throws Exception {
        NumberFormat myFormatter = NumberFormat.getInstance();
        myFormatter.setMinimumFractionDigits(2);
        myFormatter.setMaximumFractionDigits(2);
        sb.append(NEW_LINE).append("//");
        sb.append(NEW_LINE);
        sb.append("Annotating Query Reaction ").append(reactionQID).append(NEW_LINE);
        annotateReactionAsText(annotateRXNQ, reactionQID, sb);
        BondChangeCalculator bondChangeCalculatorQ = annotateRXNQ.getSelectedSolution().getBondChangeCalculator();
        sb.append(NEW_LINE).append("//");
        sb.append(NEW_LINE);
        sb.append("Annotating Target Reaction ").append(reactionTID).append(NEW_LINE);
        annotateReactionAsText(annotateRXNT, reactionTID, sb);
        BondChangeCalculator bondChangeCalculatorT = annotateRXNT.getSelectedSolution().getBondChangeCalculator();
        IPatternFingerprinter fpQ = new PatternFingerprinter();
        fpQ.add(bondChangeCalculatorQ.getFormedCleavedWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getOrderChangesWFingerprint());
        fpQ.add(bondChangeCalculatorQ.getStereoChangesWFingerprint());
        IPatternFingerprinter fpT = new PatternFingerprinter();
        fpT.add(bondChangeCalculatorT.getFormedCleavedWFingerprint());
        fpT.add(bondChangeCalculatorT.getOrderChangesWFingerprint());
        fpT.add(bondChangeCalculatorT.getStereoChangesWFingerprint());
        sb.append(NEW_LINE);
        sb.append("//");
        sb.append(NEW_LINE).append("REACTION SIMILARITY METRICS (Min:0, Max:1.0)");
        sb.append(NEW_LINE);
        double similarityBondChanges = getSimilarity(fpQ, fpT);
        sb.append("Bond Change Similarity (BC): ").append(myFormatter.format(similarityBondChanges));
        sb.append(NEW_LINE);
        double similarityReactionCentres = getSimilarity(bondChangeCalculatorQ.getReactionCenterWFingerprint(), bondChangeCalculatorT.getReactionCenterWFingerprint());
        sb.append("Reaction Centre Similarity (RC): ").append(myFormatter.format(similarityReactionCentres));
        sb.append(NEW_LINE);
        ReactionFingerprinter rfQ = new ReactionFingerprinter(bondChangeCalculatorQ.getReaction());
        ReactionFingerprinter rfT = new ReactionFingerprinter(bondChangeCalculatorT.getReaction());
        double similarityReactionStructure = getSimilarity(rfQ.getReactionStruturalFingerprint(), rfT.getReactionStruturalFingerprint());
        sb.append("Reaction Structure Similarity (ST): ").append(myFormatter.format(similarityReactionStructure));
        sb.append(NEW_LINE);
    }
}
