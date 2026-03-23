/*
 * Copyright (C) 2003-2020 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
 * GraphMatching.java
 *
 * Created on February 6, 2006, 8:20 PM
 *
 *
 * @author Syed Asad Rahman, BioInception
 * @contact asad.rahman@bioinceptionlabs.com
 *
 */
package uk.ac.ebi.reactionblast.mapping.graph;

import java.io.IOException;
import java.io.Serializable;
import java.util.BitSet;
import static java.util.Collections.unmodifiableMap;
import java.util.Map;
import java.util.logging.Level;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.mapping.interfaces.AbstractGraphMatching;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.cloneWithIDs;

/**
 *
 * @author Syed Asad Rahman, BioInception
 * @contact asad.rahman@bioinceptionlabs.com
 */
public class GraphMatching extends AbstractGraphMatching implements Serializable {

    private final static ILoggingTool LOGGER = createLoggingTool(GraphMatching.class);
    private static final long serialVersionUID = 0xf06b2d5f9L;
    private final IAtomContainer educt;
    private final IAtomContainer product;
    private IAtomContainer matchedPart = null;
    private Map<IAtom, IAtom> bestAtomMappingList;
    private int fragmentCount = 0;
    private final static boolean DEBUG = false;

    /**
     * Creates a new instance of GraphMatching
     *
     * @param reaction_ID
     * @param eductOrg
     * @param productOrg
     * @param suffix
     * @param removeHydrogen
     * @throws Exception
     */
    public GraphMatching(String reaction_ID, IAtomContainer eductOrg, IAtomContainer productOrg, String suffix, boolean removeHydrogen) throws Exception {

        try {

            educt = eductOrg;
            product = productOrg;
            educt.setID(eductOrg.getID());
            product.setID(productOrg.getID());

            if (educt.getAtomCount() > 0 && product.getAtomCount() > 0) {
                setMatchedPart(cloneWithIDs(educt));
            }
        } catch (CloneNotSupportedException e) {
            throw new CDKException("Error: In GraphMatching Class" + e);
        }

    }

    /**
     *
     * @param holder
     * @param removeHydrogen
     * @param substrateIndex
     * @param productIndex
     * @param eductFP
     * @param prodFP
     * @return
     */
    @Override
    public synchronized boolean mcsMatch(Holder holder,
            boolean removeHydrogen,
            Integer substrateIndex,
            Integer productIndex,
            BitSet eductFP,
            BitSet prodFP) {

        if (educt.getAtomCount() <= 0 && product.getAtomCount() <= 0) {
            return false;
        }

        try {
            try {
                setMCSUpdationFlags(holder, substrateIndex, productIndex);
            } catch (Exception ex) {
                LOGGER.error(Level.SEVERE, null, ex);
            }
            BestMatch initMCSAtom = holder.getBestMatchContainer();
            if (initMCSAtom.containsKey(substrateIndex, productIndex)) {
                this.bestAtomMappingList = initMCSAtom.getAtomMatch(substrateIndex, productIndex).getMappingsByAtoms();
                this.fragmentCount = initMCSAtom.getTotalFragmentCount(substrateIndex, productIndex);
                if (this.bestAtomMappingList != null && !this.bestAtomMappingList.isEmpty()) {
                    return true;
                }
            }
        } catch (IOException ex) {
            LOGGER.debug("Files: " + educt.getID() + ", " + product.getID());
            LOGGER.debug(SEVERE, null, ex);
        }
        return false;
    }

    /**
     *
     * @param reaction
     * @return
     */
    @Override
    public synchronized int removeMatchedAtomsAndUpdateAAM(IReaction reaction) {
        int delta = 0;

        LOGGER.debug("Before removing Mol Size E: " + educt.getAtomCount()
                + " , Before removing Mol Size P: " + product.getAtomCount());
        int beforeESize = educt.getAtomCount();

        if (bestAtomMappingList != null) {
            for (Map.Entry<IAtom, IAtom> map : bestAtomMappingList.entrySet()) {
                String eID = map.getKey().getID();
                IAtom eAtom = getAtomByID(educt, eID);
                String pID = map.getValue().getID();
                LOGGER.debug("eID " + eID + ",pID " + pID);
                IAtom pAtom = getAtomByID(product, pID);

                if (eAtom != null && pAtom != null) {
                    IMapping im = SilentChemObjectBuilder.getInstance().newInstance(IMapping.class, eAtom, pAtom);
                    reaction.addMapping(im);
                }
                if (eAtom != null) {
                    educt.removeAtom(eAtom);
                }
                if (pAtom != null) {
                    product.removeAtom(pAtom);
                }
                delta = fragmentCount;
            }
        }

        for (IAtom atom : educt.atoms()) {
            IAtom matchedAtom = getAtomByID(matchedPart, atom.getID());
            if (matchedAtom != null) {
                matchedPart.removeAtom(matchedAtom);
            }
        }

        LOGGER.debug("After removing Mol Size E: " + educt.getAtomCount()
                + " , After removing Mol Size P: " + product.getAtomCount());

        if (beforeESize == educt.getAtomCount()) {
            try {
                if (DEBUG) {
                    LOGGER.debug(educt.getID() + ": SMILES " + SmilesGenerator.generic().create(educt));
                    LOGGER.debug(product.getID() + ": SMILES " + SmilesGenerator.generic().create(product));
                }
                throw new CDKException("Failed to remove matched parts between " + educt.getID() + ": "
                        + educt.getAtomCount() + " , " + product.getID() + " : " + product.getAtomCount()
                        + ", Mapping count: " + bestAtomMappingList.size() + "...atom ids did not matched!");
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, "Failed to remove matched parts between " + educt.getID() + ": "
                        + educt.getAtomCount() + " , " + product.getID() + " : " + product.getAtomCount()
                        + ", Mapping count: " + bestAtomMappingList.size() + "...atom ids did not matched!", ex);

                Runtime.getRuntime().exit(1);
            }
        }
        return delta;
    }

    private synchronized IAtom getAtomByID(IAtomContainer ac, String ID) {
        if (ID == null) {
            return null;
        }
        for (IAtom atom : ac.atoms()) {
            if (ID.equalsIgnoreCase(atom.getID())) {
                return atom;
            }
        }
        return null;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized IAtomContainer getRemainingEduct() {
        return educt;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized IAtomContainer getRemainingProduct() {
        return product;
    }

    /**
     *
     * @return
     */
    protected synchronized Map<IAtom, IAtom> getFirstAtomMapping() {
        return unmodifiableMap(bestAtomMappingList);
    }

    /**
     * @return the matchedPart
     */
    @Override
    public synchronized IAtomContainer getMatchedPart() {
        return matchedPart;
    }

    /**
     * @param aMatchedPart the matchedPart to set
     */
    private synchronized void setMatchedPart(IAtomContainer aMatchedPart) {
        matchedPart = aMatchedPart;
    }
}
