package com.bioinceptionlabs.reactionblast.tools;

/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static com.bioinceptionlabs.reactionblast.mapping.MappingHandler.cleanMapping;
import static org.openscience.smsd.ExtAtomContainerManipulator.convertExplicitToImplicitHydrogens;
import com.bioinceptionlabs.reactionblast.tools.MDLRXNV2000Reader;
import static java.lang.String.valueOf;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class TestUtility {

    public static final String KEGG_RXN_DIR = "rxn/kegg/";
    public static final String RHEA_RXN_DIR = "rxn/rhea/";
    public static final String BRENDA_RXN_DIR = "rxn/brenda/";
    public static final String BUG_RXN_DIR = "rxn/bug/";
    public static final String OTHER_RXN = "rxn/other/";
    public static final String METRXN_RXN = "rxn/metrxn/";
    public static final String INFORCHEM_RXN = "rxn/infochem/";
    public static final String MACIE_RXN = "rxn/macie/";
    private final static ILoggingTool LOGGER
            = createLoggingTool(TestUtility.class);

    /**
     *
     * @param reaction
     */
    protected void renumberMappingIDs(IReaction reaction) {
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

    private InputStream getFileWithUtil(String fileName) throws IOException {
        ClassLoader classLoader = getClass().getClassLoader();
        return classLoader.getResourceAsStream(fileName);
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
    protected IReaction readReactionFile(String name, String dir, boolean reMap, boolean removeHydrogens) throws Exception {
        String filepath = dir + name + ".rxn";

        IReaction reaction = null;
        try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(getFileWithUtil(filepath))) {
            reaction = reader.read(new Reaction());
            reaction.setID(name);
            LOGGER.debug("Read Reaction ");
            for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                LOGGER.debug("r " + ac.getTitle() + ":" + reaction.getReactantCoefficient(ac));
            }
            for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                LOGGER.debug("p " + ac.getTitle() + ":" + reaction.getProductCoefficient(ac));
            }
        } catch (Exception ex) {
            LOGGER.error(SEVERE, "Unable to parse the RXN file", ex.getMessage());
        }

        if (removeHydrogens && reaction != null) {
            // XXX WARNING : this may not work correctly!
            IReaction hydrogenFreeReaction = new Reaction();
            IAtomContainerSet hydrogenFreeReactants = new AtomContainerSet();
            for (IAtomContainer atomContainer : reaction.getReactants().atomContainers()) {
                setNullHCountToZero(atomContainer);
                percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID(atomContainer.getTitle());
                hydrogenFreeReactants.addAtomContainer(acMinusH);
            }
            hydrogenFreeReaction.setReactants(hydrogenFreeReactants);
            IAtomContainerSet hydrogenFreeProducts = new AtomContainerSet();
            for (IAtomContainer atomContainer : reaction.getProducts().atomContainers()) {
                setNullHCountToZero(atomContainer);
                percieveAtomTypesAndConfigureAtoms(atomContainer);
                IAtomContainer acMinusH = convertExplicitToImplicitHydrogens(atomContainer);
                acMinusH.setID(atomContainer.getTitle());
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
}
