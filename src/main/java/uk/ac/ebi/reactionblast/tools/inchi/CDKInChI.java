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
 /*
 * CDKInChI.java
 *
 * Created on 12 January 2007, 14:55
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast.tools.inchi;

import java.util.List;
import net.sf.jniinchi.INCHI_OPTION;
import net.sf.jniinchi.INCHI_RET;
import static net.sf.jniinchi.INCHI_RET.OKAY;
import static net.sf.jniinchi.INCHI_RET.WARNING;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import static org.openscience.cdk.inchi.InChIGeneratorFactory.getInstance;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @author Syed Asad Rahman, EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
public class CDKInChI {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CDKInChI.class);
    /*
    On suggestion from D. Schomburg as 'At' a radioactive halogen that never appears in nature
     */

    /**
     *
     */
    public static final String R_Group_replacement_String = "At";

    /**
     *
     */
    public static final String[] METALS = {"At", "Th", "Pa", "U", "Np", "Pu", "Am",
        "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};

    private InChIGenerator _genInchi;
    private final InChIToStructure _intostruct;
    private final InChIGeneratorFactory factory;

    /**
     *
     */
    protected IAtomContainer molecule;

    /**
     * Creates a new instance of CDKInChI
     *
     * @throws CDKException
     */
    public CDKInChI() throws CDKException {
        factory = getInstance();
        _genInchi = null;
        _intostruct = null;
        molecule = null;
    }

    private synchronized IAtomContainer convertRGroupsToMetals(IAtomContainer mol) {
        try {
            IAtomContainer convertedMol = mol.clone();
            for (IAtom atom : convertedMol.atoms()) {
                if (isR(atom)) {
                    atom.setSymbol(R_Group_replacement_String);
                }
            }
            return convertedMol;
        } catch (CloneNotSupportedException c) {
            return mol;
        }
    }

    private synchronized boolean isR(IAtom atom) {
        String symbol = atom.getSymbol();
        return symbol.charAt(0) == 'R'
                && !(symbol.equals("Re")
                || symbol.equals("Rh")
                || symbol.equals("Rn")
                || symbol.equals("Ru"));
    }

    /**
     *
     * @param mol
     * @throws org.openscience.cdk.exception.CDKException
     * @return
     */
    public synchronized String getInChI(IAtomContainer mol) throws CDKException {
        String inchi = "";
        try {
            _genInchi = factory.getInChIGenerator(convertRGroupsToMetals(mol), "");
            inchi();
            inchi = _genInchi.getInchi();
        } catch (CDKException e) {
            LOGGER.debug("Error in generating InChI code " + e);
        }
        return inchi;
    }

    /**
     *
     * @param mol
     * @param options
     * @throws org.openscience.cdk.exception.CDKException
     * @return
     */
    public synchronized String getInChI(IAtomContainer mol, String options) throws CDKException {
        String inchi = "";
        // Generate factory - throws CDKException if native code does not load
        try {
            _genInchi = factory.getInChIGenerator(convertRGroupsToMetals(mol), options);
            inchi();
            inchi = _genInchi.getInchi();
        } catch (CDKException e) {
            LOGGER.debug("Error in generating InChI code " + e);
        }
        return inchi;
    }

    /**
     *
     * @param mol
     * @param options
     * @throws org.openscience.cdk.exception.CDKException
     * @return
     */
    public synchronized String getInChI(IAtomContainer mol, List<INCHI_OPTION> options) throws CDKException {
        String inchi = "";

        // Generate factory - throws CDKException if native code does not load
        try {
            _genInchi = factory.getInChIGenerator(convertRGroupsToMetals(mol), options);
            inchi();
            inchi = _genInchi.getInchi();
        } catch (CDKException e) {
            LOGGER.debug("Error in generating InChI code " + e);
        }
        return inchi;
    }

    private synchronized void structure() throws CDKException {

        INCHI_RET ret = _intostruct.getReturnStatus();

        if (ret == WARNING) {
            // Structure generated, but with warning message
            LOGGER.debug("InChI warning: " + _intostruct.getMessage());
        } else if (ret != OKAY) {
            // Structure generation failed
            throw new CDKException("Structure generation failed failed: " + ret.toString() + " [" + _intostruct.getMessage() + "]");
        }
    }

    private synchronized void inchi() throws CDKException {

        INCHI_RET ret = _genInchi.getReturnStatus();
        if (ret == WARNING) {
            // CDKInChI generated, but with warning message
//            System.LOGGER.debug("InChI warning: " + _genInchi.getMessage());
        } else if (ret != OKAY) {
            // CDKInChI generation failed
            throw new CDKException("InChI failed: " + ret.toString() + " [" + _genInchi.getMessage() + "]");
        }

    }

    /**
     *
     * @return returns Auxilary information about the CDKInChI
     */
    public synchronized String getAuxinfo() {
        return _genInchi.getAuxInfo();
    }
}
