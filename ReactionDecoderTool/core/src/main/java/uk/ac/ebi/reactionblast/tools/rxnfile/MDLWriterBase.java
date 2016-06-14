/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.BufferedWriter;
import java.io.IOException;
import static java.lang.String.valueOf;
import java.text.NumberFormat;
import static java.text.NumberFormat.getNumberInstance;
import java.util.List;
import static java.util.Locale.ENGLISH;
import static org.openscience.cdk.CDKConstants.COMMENT;
import static org.openscience.cdk.config.Isotopes.getInstance;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.DefaultChemObjectWriter;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class MDLWriterBase extends DefaultChemObjectWriter {

    private static final String M = "M  ";

    /**
     * Formats an integer to fit into the connection table and changes it to a
     * String.
     *
     * @param i The int to be formated
     * @param l Length of the String
     * @return The String to be written into the connection table
     */
    protected static String formatMDLInt(int i, int l) {
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

    /**
     * Formats a float to fit into the connection table and changes it to a
     * String.
     *
     * @param fl The float to be formated
     * @return The String to be written into the connection table
     */
    protected static String formatMDLFloat(float fl) {
        String s = "", fs = "";
        int l;
        NumberFormat nf = getNumberInstance(ENGLISH);
        nf.setMinimumIntegerDigits(1);
        nf.setMaximumIntegerDigits(4);
        nf.setMinimumFractionDigits(4);
        nf.setMaximumFractionDigits(4);
        nf.setGroupingUsed(false);
        s = nf.format(fl);
        l = 10 - s.length();
        for (int f = 0; f < l; f++) {
            fs += " ";
        }
        fs += s;
        return fs;
    }

    /**
     * Formats a String to fit into the connection table.
     *
     * @param s The String to be formated
     * @param le The length of the String
     * @return The String to be written in the connection table
     */
    protected static String formatMDLString(String s, int le) {
        s = s.trim();
        if (s.length() > le) {
            return s.substring(0, le);
        }
        int l;
        l = le - s.length();
        for (int f = 0; f < l; f++) {
            s += " ";
        }
        return s;
    }
    /**
     * SuperAtomContainer
     */
    protected SuperAtomContainer superAtomContainer;

    /**
     * Returns a <code>String</code> to which is appended enough space (" ")
     * characters to make the total length equal to <code>length</code>. If the
     * length of <code>string</code> is greater than <code>length</code>, then
     * <code>string</code> is returned with no modifications.
     *
     * @param string a string to operate on
     * @param length the desired final length of the string
     * @return a string with spaces appended to it, or just <code>string</code>
     */
    public String padRight(String string, int length) {
        String result = string;
        int spaceCount = length - string.length();

        for (int i = 0; i < spaceCount; i++) {
            result += ' ';
        }

        return result;
    }

    /**
     * Returns a <code>String</code> to which is prepended enough space (" ")
     * characters to make the total length equal to <code>length</code>. If the
     * length of <code>string</code> is greater than <code>length</code>, then
     * <code>string</code> is returned with no modifications.
     *
     * @param string a string to operate on
     * @param length the desired final length of the string
     * @return a string with spaces prepended to it, or just <code>string</code>
     */
    public String padLeft(String string, int length) {
        String result = string;
        int spaceCount = length - string.length();

        for (int i = 0; i < spaceCount; i++) {
            result = ' ' + result;
        }

        return result;
    }

    /**
     *
     * @param writer
     * @param container
     * @throws IOException
     */
    protected void writeAtomValue(BufferedWriter writer, IAtomContainer container) throws IOException {
        // Write Atom Value
        for (int i = 0; i < container.getAtomCount(); i++) {
            IAtom atom = container.getAtom(i);
            if (atom.getProperty(COMMENT) != null && atom.getProperty(COMMENT) instanceof String && !((String) atom.getProperty(COMMENT)).trim().isEmpty()) {
                writer.write("V  ");
                writer.write(formatMDLInt(i + 1, 3));
                writer.write(" ");
                writer.write((String) atom.getProperty(COMMENT));
                writer.newLine();
            }
        }
    }

    /**
     *
     * @param writer
     * @param container
     * @throws IOException
     */
    protected void writeFormalAtomicCharges(BufferedWriter writer, IAtomContainer container) throws IOException {
        // write formal atomic charges
        for (int i = 0; i < container.getAtomCount(); i++) {
            IAtom atom = container.getAtom(i);
            Integer charge = atom.getFormalCharge();
            if (charge != null && charge != 0) {
                writer.write("M  CHG  1 ");
                writer.write(formatMDLInt(i + 1, 3));
                writer.write(" ");
                writer.write(formatMDLInt(charge, 3));
                writer.newLine();
            }
        }
    }

    /**
     *
     * @param writer
     * @param container
     * @throws IOException
     */
    protected void writeFormalIsotope(BufferedWriter writer, IAtomContainer container) throws IOException {
        // write formal isotope information
        for (int i = 0; i < container.getAtomCount(); i++) {
            IAtom atom = container.getAtom(i);
            if (!(atom instanceof IPseudoAtom)) {
                Integer atomicMass = atom.getMassNumber();
                if (atomicMass != null) {
                    int majorMass = getInstance().getMajorIsotope(atom.getSymbol()).getMassNumber();
                    if (atomicMass != majorMass) {
                        writer.write("M  ISO  1 ");
                        writer.write(formatMDLInt(i + 1, 3));
                        writer.write(" ");
                        writer.write(formatMDLInt(atomicMass, 3));
                        writer.newLine();
                    }
                }
            }
        }
    }

    /**
     *
     * @param writer
     * @param rgroupList
     * @throws IOException
     */
    protected void writeRGroups(BufferedWriter writer, List<Integer> rgroupList) throws IOException {
        //write RGP line (max occurrence is 16 data points per line)
        if (rgroupList != null) {
            StringBuffer rgpLine = new StringBuffer();
            int cnt = 0;
            for (int i = 1; i <= rgroupList.size(); i++) {
                rgpLine.append(formatMDLInt(rgroupList.get(i - 1), 4));
                i++;
                rgpLine.append(formatMDLInt(rgroupList.get(i - 1), 4));
                cnt++;
                if (i == rgroupList.size() || i == 16) {
                    rgpLine.insert(0, "M  RGP" + formatMDLInt(cnt, 3));
                    writer.write(rgpLine.toString());
                    writer.newLine();
                    rgpLine = new StringBuffer();
                    cnt = 0;
                }
            }
        }
    }

    private void writeSgroupCount(IAtomContainer molecule, BufferedWriter writer) throws IOException {
        int substructureCount = superAtomContainer.countSuperatoms();
        writer.write(M + "STY" + padLeft(valueOf(substructureCount), 3));
        for (int i = 0; i < substructureCount; i++) {
            writer.write(padLeft(valueOf(i + 1), 4));
            writer.write(padLeft("SUP", 4));
        }
        writer.newLine();

        writer.write(M + "SLB" + padLeft(valueOf(substructureCount), 3));
        for (int i = 0; i < substructureCount; i++) {
            writer.write(padLeft(valueOf(i + 1), 4));
            writer.write(padLeft(valueOf(i + 1), 4));
        }
        writer.newLine();

    }

    /**
     *
     * @param molecule
     * @param writer
     * @throws IOException
     */
    protected void writeSgroupProperty(IAtomContainer molecule, BufferedWriter writer) throws IOException {
        if (superAtomContainer.countSuperatoms() > 0) {
            writeSgroupCount(molecule, writer);
            for (int i = 0; i < superAtomContainer.countSuperatoms(); i++) {
                writeSingleSgroup(superAtomContainer.getSubstructures(i), writer);
            }
        }
    }

    private void writeSingleSgroup(Substructure substructure, BufferedWriter writer) throws IOException {
        writer.write(M + "SAL" + padLeft(valueOf(substructure.getIndex() + 1), 4)
                + padLeft(valueOf(substructure.getSuperAtomCount()), 3));
        for (int i = 0; i < substructure.getSuperAtomCount(); i++) {
            writer.write(padLeft(valueOf(substructure.getSuperAtom(i).getIndex() + 1), 4));
        }
        writer.newLine();
        writer.write(M + "SBL" + padLeft(valueOf(substructure.getIndex() + 1), 4) + padLeft(valueOf(substructure.getSuperBondCount()), 3));
        for (int i = 0; i < substructure.getSuperBondCount(); i++) {
            writer.write(padLeft(valueOf(substructure.getSuperBond(i).getBondIndex() + 1), 4));
        }
        writer.newLine();
        for (int i = 0; i < substructure.getSuperAtomCount(); i++) {
            writer.write(M + "SMT" + padLeft(valueOf(substructure.getIndex() + 1), 4)
                    + " " + substructure.getSuperAtom(i).getPseudoAtom().getLabel());
        }
        writer.newLine();
        /*
         * Presently not supported
         //        writer.write(M + "SCL" + padLeft(String.valueOf(substructure.getIndex() + 1), 4) + "  ");
         //        writer.newLine();
         //        for (int i = 0; i < substructure.getSuperBondCount(); i++) {
         //            SuperBonds crossingBond = substructure.getSuperBond(i);
         //            writer.write(M + "SBV" + padLeft(String.valueOf(substructure.getIndex() + 1), 4));
         //            writer.write(padLeft(String.valueOf(crossingBond.getBondIndex() + 1), 4));
         //            writer.write(padLeft(formatMDLFloat((float) crossingBond.getX()), 10));
         //            writer.write(padLeft(formatMDLFloat((float) crossingBond.getY()), 10));
         //            writer.newLine();
         //        }
         * 
         */
    }
}
