/* Copyright (C) 1997-2007  Christoph Steinbeck <steinbeck@users.sourceforge.net>
 *                    2011  Egon Willighagen <egonw@users.sourceforge.net>
 *                    2011  Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.BufferedReader;
import java.io.IOException;
import static java.lang.Float.parseFloat;
import static java.lang.Integer.parseInt;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.replaceAtomByAtom;

/**
 * @cdk.module io
 * 
 * @cdk.keyword file format, MDL molfile
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class MDLReaderBase extends DefaultChemObjectReader {

    private static final ILoggingTool logger
            = createLoggingTool(MDLReaderBase.class);

    /**
     *
     * @param molecule
     * @param prevAtom
     * @param pseudoAtom
     */
    public static void replaceAtom(IAtomContainer molecule, IAtom prevAtom, IPseudoAtom pseudoAtom) {
        if (prevAtom.getPoint2d() != null) {
            pseudoAtom.setPoint2d(prevAtom.getPoint2d());
        }
        if (prevAtom.getPoint3d() != null) {
            pseudoAtom.setPoint3d(prevAtom.getPoint3d());
        }
        replaceAtomByAtom(molecule, prevAtom, pseudoAtom);
    }

    /**
     *
     */
    protected SuperAtomContainer superAtomContainer;

    /*
     * M STY (Sgroup Type)
     * M SLB (Sgroup Label)
     * M SAL (Sgroup Atom list)
     * M SBL (Sgroup Bond list)
     * M SMT (Sgroup SubScript)
     * M SBV
     */

    /**
     *
     * @param molecule
     * @param line
     * @param linecount
     * @throws CDKException
     */

    protected void createAtomProperty(IAtomContainer molecule, String line, int linecount) throws CDKException {
        if (line.startsWith("M  CHG")) {
            // FIXME: if this is encountered for the first time, all
            // atom charges should be set to zero first!
            int infoCount = parseInt(extractString(line, 6, 9));
            StringTokenizer st = new StringTokenizer(line.substring(9));
            for (int i = 1; i <= infoCount; i++) {
                String token = st.nextToken();
                int atomNumber = parseInt(token.trim());
                token = st.nextToken();
                int charge = parseInt(token.trim());
                molecule.getAtom(atomNumber - 1).setFormalCharge(charge);
            }
        } else if (line.startsWith("M  ISO")) {
            try {
                String countString = extractString(line, 6, 10);
                int infoCount = parseInt(countString);
                StringTokenizer st = new StringTokenizer(line.substring(10));
                for (int i = 1; i <= infoCount; i++) {
                    int atomNumber = parseInt(st.nextToken().trim());
                    int absMass = parseInt(st.nextToken().trim());
                    if (absMass != 0) {
                        IAtom isotope = molecule.getAtom(atomNumber - 1);
                        isotope.setMassNumber(absMass);
                    }
                }
            } catch (NumberFormatException exception) {
                String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": " + line + " in property block.";
                logger.error(error);
                handleError("NumberFormatException in isotope information.", linecount, 7, 11, exception);
            }
        } else if (line.startsWith("M  RAD")) {
            try {
                String countString = extractString(line, 6, 9);
                int infoCount = parseInt(countString);
                StringTokenizer st = new StringTokenizer(line.substring(9));
                for (int i = 1; i <= infoCount; i++) {
                    int atomNumber = parseInt(st.nextToken().trim());
                    int spinMultiplicity = parseInt(st.nextToken().trim());
                    if (spinMultiplicity > 1) {
                        IAtom radical = molecule.getAtom(atomNumber - 1);
                        for (int j = 2; j <= spinMultiplicity; j++) {
                            // 2 means doublet -> one unpaired electron
                            // 3 means triplet -> two unpaired electron
                            molecule.addSingleElectron(molecule.getBuilder().newInstance(ISingleElectron.class, radical));
                        }
                    }
                }
            } catch (NumberFormatException exception) {
                String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": " + line + " in property block.";
                logger.error(error);
                handleError("NumberFormatException in radical information", linecount, 7, 10, exception);
            }
        }
    }

    /**
     *
     * @param input
     * @param outputContainer
     * @param molecule
     * @param line
     * @param linecount
     * @throws IOException
     * @throws CDKException
     */
    protected void createGroupOldVersion(BufferedReader input, IAtomContainer outputContainer, IAtomContainer molecule, String line, int linecount) throws IOException, CDKException {
        try {
            String atomNumberString = extractString(line, 3, 6);
            int atomNumber = parseInt(atomNumberString);
            //String whatIsThisString = line.substring(6,9).trim();
            String atomName = input.readLine();
            // convert Atom into a PseudoAtom
            IAtom prevAtom = outputContainer.getAtom(atomNumber - 1);
            IPseudoAtom pseudoAtom = molecule.getBuilder().newInstance(IPseudoAtom.class, prevAtom);
            pseudoAtom.setSymbol(prevAtom.getSymbol());
            pseudoAtom.setLabel(atomName);
            if (prevAtom.getPoint2d() != null) {
                pseudoAtom.setPoint2d(prevAtom.getPoint2d());
            }
            if (prevAtom.getPoint3d() != null) {
                pseudoAtom.setPoint3d(prevAtom.getPoint3d());
            }
            replaceAtomByAtom(molecule, prevAtom, pseudoAtom);
        } catch (NumberFormatException exception) {
            String error = "Error (" + exception.toString() + ") while parsing line " + linecount + ": " + line + " in property block.";
            logger.error(error);
            handleError("NumberFormatException in group information", linecount, 4, 7, exception);
        }
    }

    /*
     * M STY (Sgroup Type)
     * M SLB (Sgroup Label)
     * M SAL (Sgroup Atom list)
     * M SBL (Sgroup Bond list)
     * M SMT (Sgroup SubScript)
     * M SBV
     */

    /**
     *
     * @param molecule
     * @param outputContainer
     * @param line
     * @param linecount
     * @throws CDKException
     */

    protected void createSgroupProperty(IAtomContainer molecule, IAtomContainer outputContainer, String line, int linecount) throws CDKException {
        try {
            String property = extractString(line, 0, 6);
            switch (property) {
                case "M  STY": {
                    int entryCount = extractInt(line, 6, 9);
                    superAtomContainer = new SuperAtomContainer(entryCount);
                    for (int i = 0; i < entryCount; i++) {
                        int offset = i * 8;
                        int sgroupIndex = extractInt(line, offset + 10, offset + 13);
                        String sgroupType = extractString(line, offset + 14, offset + 17);
                        if (!"SUP".equals(sgroupType)) {
                            throw new RuntimeException("Error parsing Substructure, only Superatom is supported.");
                        }
                        Substructure sub = new Substructure((sgroupIndex - 1), sgroupType);
                        superAtomContainer.add(sub);
                    }
                    break;
                }
                case "M  SLB": {
                    //        Ignore Identifier
                    int entryCount = extractInt(line, 6, 9);
                    for (int i = 0; i < entryCount; i++) {
                        int offset = i * 8;
                        int sgroupIndex = extractInt(line, offset + 10, offset + 13);
                        int sgroupIdentifier = extractInt(line, offset + 14, offset + 17);
                        Substructure substructure = superAtomContainer.getSubstructures(sgroupIndex - 1);
                        substructure.setSGroupIdentifier(sgroupIdentifier);
                    }
                    break;
                }
                case "M  SAL": {
                    int sgroupIndex = extractInt(line, 6, 10);
                    Substructure substructure = superAtomContainer.getSubstructures(sgroupIndex - 1);
                    int entryCount = extractInt(line, 10, 13);
                    for (int i = 0; i < entryCount; i++) {
                        int offset = i * 4;
                        int atomIndex = extractInt(line, offset + 13, offset + 17);
                        SuperAtoms superAtom = new SuperAtoms();
                        superAtom.setIndex(atomIndex - 1);
                        superAtom.setPrevAtom(outputContainer.getAtom(atomIndex - 1));
                        substructure.add(superAtom);
                    }
                    break;
                }
                case "M  SBL": {
                    int sgroupIndex = extractInt(line, 6, 10);
                    Substructure substructure = superAtomContainer.getSubstructures(sgroupIndex - 1);
                    int entryCount = extractInt(line, 10, 13);
                    for (int i = 0; i < entryCount; i++) {
                        int offset = i * 4;
                        int bondIndex = extractInt(line, offset + 13, offset + 17);
                        substructure.addCrossingBond((bondIndex - 1), outputContainer.getBond(bondIndex - 1));
                    }
                    break;
                }
                case "M  SMT": {
                    int sgroupIndex = extractInt(line, 6, 10);
                    Substructure substructure = superAtomContainer.getSubstructures(sgroupIndex - 1);
                    String label = extractString(line, 10, line.length());
                    for (SuperAtoms sa : substructure.atoms()) {
                        IPseudoAtom pseudoAtom = molecule.getBuilder().newInstance(IPseudoAtom.class);
                        pseudoAtom.setLabel(label);
                        pseudoAtom.setSymbol(sa.getPrevAtom().getSymbol());
                        sa.setPseudoAtom(pseudoAtom);
                        replaceAtom(outputContainer, sa.getPrevAtom(), sa.getPseudoAtom());
                    }
                    break;
                }
                case "M  SBV": {
                    int sgroupIndex = extractInt(line, 6, 10);
                    Substructure substructure = superAtomContainer.getSubstructures(sgroupIndex - 1);
                    int bondIndex = extractInt(line, 10, 14);
                    double x = extractFloat(line, 14, 24);
                    double y = extractFloat(line, 24, 34);
                    substructure.setCrossingVector(outputContainer.getBond(bondIndex - 1), x, y);
                    break;
                }
            }
        } catch (NumberFormatException exception) {
            String error = "Error (" + exception.getMessage() + ") while parsing line "
                    + linecount + ": " + line + " in property block.";
            logger.error(error);
            handleError(
                    "NumberFormatException in radical information",
                    linecount, 7, 10,
                    exception);
        }
    }

    /**
     *
     * @param molecule
     * @param outputContainer
     * @param line
     * @param RGroupCounter
     * @param aliasAtomNumber
     */
    protected void createPseudoAtomProperty(IAtomContainer molecule, IAtomContainer outputContainer, String line, int RGroupCounter, int aliasAtomNumber) {
        // Reads the pseudo atom property from the mol file
        String[] aliasArray = line.split("\\\\");
        // name of the alias atom like R1 odr R2 etc.
        String alias = "";
        for (String aliasArray1 : aliasArray) {
            alias += aliasArray1;
        }
        IAtom aliasAtom = outputContainer.getAtom(aliasAtomNumber);
        // skip if already a pseudoatom
        if (aliasAtom instanceof IPseudoAtom) {
            ((IPseudoAtom) aliasAtom).setLabel(alias);
            return;
        }
        IAtom newPseudoAtom = molecule.getBuilder().newInstance(IPseudoAtom.class, alias);
        if (aliasAtom.getPoint2d() != null) {
            newPseudoAtom.setPoint2d(aliasAtom.getPoint2d());
        }
        if (aliasAtom.getPoint3d() != null) {
            newPseudoAtom.setPoint3d(aliasAtom.getPoint3d());
        }
        outputContainer.addAtom(newPseudoAtom);
        List<IBond> bondsOfAliasAtom = outputContainer.getConnectedBondsList(aliasAtom);
        for (int i = 0; i < bondsOfAliasAtom.size(); i++) {
            IBond bondOfAliasAtom = bondsOfAliasAtom.get(i);
            IAtom connectedToAliasAtom = bondOfAliasAtom.getOther(aliasAtom);
            IBond newBond = bondOfAliasAtom.getBuilder().newInstance(IBond.class);
            newBond.setAtoms(new IAtom[]{connectedToAliasAtom, newPseudoAtom});
            newBond.setOrder(bondOfAliasAtom.getOrder());
            outputContainer.addBond(newBond);
            outputContainer.removeBond(aliasAtom, connectedToAliasAtom);
        }
        outputContainer.removeAtom(aliasAtom);
        RGroupCounter++;
    }

    /**
     *
     * @param line
     * @param rAtoms
     * @param Rnumber
     */
    protected void createRGroupAtomProperties(String line, Map<Integer, IPseudoAtom> rAtoms, int Rnumber) {
        StringTokenizer st = new StringTokenizer(line);
        //Ignore first 3 tokens (overhead).
        st.nextToken();
        st.nextToken();
        st.nextToken();
        //Process the R group numbers as defined in RGP line.
        while (st.hasMoreTokens()) {
            Integer position = new Integer(st.nextToken());
            Rnumber = new Integer(st.nextToken());
            IPseudoAtom pseudoAtom = rAtoms.get(position);
            if (pseudoAtom != null) {
                pseudoAtom.setLabel("R" + Rnumber);
            }
        }
    }

    /**
     *
     */
    public abstract void customizeJob();

    /**
     * Returns the float value of the specified substring.
     *
     * @param string the string
     * @param start the start index
     * @param stop the end index
     *
     * @return the float value
     */
    public float extractFloat(String string, int start, int stop) {
        string = string.substring(start, stop);
        string = string.trim();
        float value = 0;
        if ("".equals(string)) {
            return 0;
        }
        try {
            value = parseFloat(string);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Error trying to parse: " + string + " as a float.");
        }
        return value;
    }

    /**
     * Returns the integer value of the specified substring.
     *
     * @param string the string
     * @param start the start index
     * @param stop the end index
     *
     * @return the int value
     */
    public int extractInt(String string, int start, int stop) {
        string = string.substring(start, stop);
        string = string.trim();
        int value = 0;
        if ("".equals(string)) {
            return 0;
        }
        try {
            value = parseInt(string);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Error trying to parse: " + string + " as an integer.");
        }
        return value;
    }

    /**
     * Returns the specified substring, trimming of whitespace.
     *
     * @param string the string
     * @param start the start index
     * @param stop the end index
     *
     * @return the substring
     */
    public String extractString(String string, int start, int stop) {
        string = string.substring(start, stop);
        string = string.trim();
        return string;
    }

}
