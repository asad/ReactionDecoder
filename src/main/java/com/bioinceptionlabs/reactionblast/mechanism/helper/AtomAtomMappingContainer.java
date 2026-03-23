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
package com.bioinceptionlabs.reactionblast.mechanism.helper;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.reactionblast.mapping.Reactor;

/**
 * This class maintain the atom-atom mapping of a reaction between reactants and
 * products
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * @author Lorenzo Baldacci {lorenzo@ebi.ac.uk|lbaldacc@csr.unibo.it}
 */
public class AtomAtomMappingContainer extends Object implements Serializable {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(AtomAtomMappingContainer.class);

    private static final String NEW_LINE = System.lineSeparator();
    private static final long serialVersionUID = 17879096958755L;

    private List<IAtom> reactantAtomArray = new ArrayList<>();
    private List<IAtom> productAtomArray = new ArrayList<>();
//    private Reactor myReaction = null;

    /**
     * Class constructor. Creates the mapping of a given reaction.
     *
     * @param reactor Reactor for which the AtomAtomMappingContainer is required
     * @param withoutH Store Mapping without H
     * @throws Exception
     */
    public AtomAtomMappingContainer(Reactor reactor, boolean withoutH) throws Exception {
        this(reactor.getReactionWithAtomAtomMapping(), withoutH);
    }

    /**
     *
     * @param reactants
     * @param products
     * @param withoutH
     */
    public AtomAtomMappingContainer(IAtomContainerSet reactants, IAtomContainerSet products, boolean withoutH) {
        int atomNo = 0;
        int mappedAtomsR = 0;
        int mappedAtomsP = 0;
        IAtom[] atVect = null;

        //REACTANTS
        if (withoutH) {
            for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                IAtomContainer M = reactants.getAtomContainer(i);
                for (IAtom a : M.atoms()) {
                    if (!a.getSymbol().equalsIgnoreCase("H")) {
                        atomNo++;
                    }
                }
            }
        } else {
            for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                atomNo += reactants.getAtomContainer(i).getAtomCount();
            }
        }
        atVect = new IAtom[atomNo];
        for (int i = 0; i < atVect.length; i++) {
            atVect[i] = null;
        }
        for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
            for (int j = 0; j < reactants.getAtomContainer(i).getAtomCount(); j++) {
                IAtom at = reactants.getAtomContainer(i).getAtom(j);
                if (withoutH && at.getSymbol().equalsIgnoreCase("H")) {
                    continue;
                }
                if (at.getID() == null) {
                    continue;
                }
                int atomID;
                try {
                    atomID = Integer.valueOf(at.getID());
                } catch (NumberFormatException e) {
                    continue;
                }
                if (atomID <= 0) {
                    continue;
                }
                if (atomID - 1 >= atVect.length) {
                    continue;
                }
                atVect[atomID - 1] = at;
                mappedAtomsR++;
            }
        }
        for (int i = 0; i < mappedAtomsR; i++) {
            reactantAtomArray.add(atVect[i]);
        }
        //Checking for holes in the vector. 
        boolean findNull = false;
        boolean error = false;
        for (IAtom atVect1 : atVect) {
            if (findNull && (atVect1 != null)) {
                error = true;
            }
            if (atVect1 == null) {
                findNull = true;
            }
        }
        if (error) {
            StringBuilder sb = new StringBuilder();
            sb.append("ERROR in AtomAtomMapping-found hole in the mapping (reactants atomIDs)");
            for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
                sb.append(NEW_LINE).append("Mol:").append(reactants.getAtomContainer(i).getID());
                for (int j = 0; j < reactants.getAtomContainer(i).getAtomCount(); j++) {
                    IAtom at = reactants.getAtomContainer(i).getAtom(j);
                    sb.append(NEW_LINE).append(at.getSymbol()).append(at.getID());
                }
            }
            LOGGER.error(sb.toString());
        }
        //end of checking statements

        //PRODUCTS
        atomNo = 0;
        if (withoutH) {
            for (int i = 0; i < products.getAtomContainerCount(); i++) {
                IAtomContainer M = products.getAtomContainer(i);
                for (IAtom a : M.atoms()) {
                    if (!a.getSymbol().equalsIgnoreCase("H")) {
                        atomNo++;
                    }
                }
            }
        } else {
            for (int i = 0; i < products.getAtomContainerCount(); i++) {
                atomNo += products.getAtomContainer(i).getAtomCount();
            }
        }
        atVect = new IAtom[atomNo];
        for (int i = 0; i < atVect.length; i++) {
            atVect[i] = null;
        }
        for (int i = 0; i < products.getAtomContainerCount(); i++) {
            for (int j = 0; j < products.getAtomContainer(i).getAtomCount(); j++) {
                IAtom at = products.getAtomContainer(i).getAtom(j);
                if (withoutH && at.getSymbol().equalsIgnoreCase("H")) {
                    continue;
                }

                if (at.getID() == null) {
                    continue;
                }
                int atomID;
                try {
                    atomID = Integer.valueOf(at.getID());
                } catch (NumberFormatException e) {
                    continue;
                }
                if (atomID <= 0) {
                    continue;
                }
                if (atomID - 1 >= atVect.length) {
                    continue;
                }
                atVect[atomID - 1] = at;
                mappedAtomsP++;
            }
        }
        for (int i = 0; i < mappedAtomsP; i++) {
            productAtomArray.add(atVect[i]);
        }
        //Checking for holes in the vector. 
        findNull = false;
        error = false;
        for (IAtom atVect1 : atVect) {
            if (findNull && (atVect1 != null)) {
                error = true;
            }
            if (atVect1 == null) {
                findNull = true;
            }
        }
        if (mappedAtomsP != mappedAtomsR) {
            error = true;
        }
        if (error) {
            StringBuilder sb = new StringBuilder();
            sb.append("ERROR in AtomAtomMapping-found hole in the mapping (products atomIDs)");
            sb.append("mapped reactants atoms: ").append(mappedAtomsR).append(", mapped products atoms: ").append(mappedAtomsP);
            for (int i = 0; i < products.getAtomContainerCount(); i++) {
                sb.append(NEW_LINE).append("Mol:").append(products.getAtomContainer(i).getID());
                for (int j = 0; j < products.getAtomContainer(i).getAtomCount(); j++) {
                    IAtom at = products.getAtomContainer(i).getAtom(j);
                    sb.append(NEW_LINE).append(at.getSymbol()).append(at.getID());
                }
            }
            LOGGER.error(sb.toString());
        }
        //end of checking statements
    }

    /**
     *
     * @param reaction
     * @param withoutH
     */
    public AtomAtomMappingContainer(IReaction reaction, boolean withoutH) {
        for (IMapping m : reaction.mappings()) {
            IAtom rAtom = (IAtom) m.getChemObject(0);
            IAtom pAtom = (IAtom) m.getChemObject(1);
            if (withoutH && rAtom != null && pAtom != null
                    && (rAtom.getSymbol().equalsIgnoreCase("H")
                    || pAtom.getSymbol().equalsIgnoreCase("H"))) {
            } else {
                reactantAtomArray.add(rAtom);
                productAtomArray.add(pAtom);
            }
        }
    }

    /**
     * This method prints the matrix to the standard output
     *
     * @return
     */
    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append(reactantAtomArray.size()).append(NEW_LINE);
        for (int i = 0; i < reactantAtomArray.size(); i++) {
            result.append(i).append("\t");
        }
        result.append(NEW_LINE);
        for (int i = 0; i < reactantAtomArray.size(); i++) {
            result.append((reactantAtomArray.get(i)).getSymbol()).append((reactantAtomArray.get(i)).getID()).append("\t");
        }
        result.append(NEW_LINE);
        for (int i = 0; i < productAtomArray.size(); i++) {
            result.append((productAtomArray.get(i)).getSymbol()).append((productAtomArray.get(i)).getID()).append("\t");
        }
        result.append(NEW_LINE);
        return result.toString();
    }

    /**
     * The method returns the product atom mapped to the reactant atom passed as
     * parameter. Returns null if the reactantAtom is not mapped to any product
     * atom.
     *
     * @param reactantAtom The IAtom for which the product atom is required.
     * @return The product atom mapped to the given reactant atom.
     */
    public IAtom getMappedProductAtom(IAtom reactantAtom) {
        IAtom a = null;
        int reactantIdx = -1;
        for (int i = 0; i < reactantAtomArray.size(); i++) {
            if (reactantAtomArray.get(i).getID().equals(reactantAtom.getID())) {
                reactantIdx = i;
            }
        }
        if (reactantIdx != -1) {
            a = productAtomArray.get(reactantIdx);
        }
        return a;
    }

    /**
     * The method returns the idx-th reactant atom which has been mapped.
     *
     * @param idx The index of the reactant atom which is required.
     * @return The idx-th a reactant atom mapped.
     */
    public IAtom getReactantAtom(int idx) {
        IAtom ret = null;
        if ((idx < reactantAtomArray.size()) && (idx > -1)) {
            ret = reactantAtomArray.get(idx);
        }
        return ret;
    }

    /**
     *
     * @param idx
     * @return
     */
    public IAtom getProductAtom(int idx) {
        IAtom ret = null;
        if ((idx < productAtomArray.size()) && (idx > -1)) {
            ret = productAtomArray.get(idx);
        }
        return ret;
    }

    /**
     * Returns the number of mappings which the AtomAtomMappingContainer
     * contains.
     *
     * @return the number of mappings which the AtomAtomMappingContainer
     * contains.
     */
    public int getSize() {
        return reactantAtomArray.size();
    }

    /**
     *
     * @return
     */
    public int getSizeNoHydrogens() {
        int count = 0;
        count = reactantAtomArray.stream().filter((a)
                -> (!a.getSymbol().equals("H"))).map((_item) -> 1)
                .reduce(count, Integer::sum);
        return count;
    }

    /**
     * Returns true if the reactant atom is present
     *
     * @param atom
     * @return
     */
    public boolean isReactantAtomPresent(IAtom atom) {
        return reactantAtomArray.contains(atom) == true;
    }

    /**
     * Return true if the product atom is present
     *
     * @param atom
     * @return
     */
    public boolean isProductAtomPresent(IAtom atom) {
        return productAtomArray.contains(atom) == true;
    }
}
