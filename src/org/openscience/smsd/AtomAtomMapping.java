/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received query copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * 
 */
package org.openscience.smsd;

import com.google.common.collect.HashBiMap;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 * Holds atom-atom mappings information between source and target molecules
 *
 *  
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class AtomAtomMapping implements Serializable {

    private static final long serialVersionUID = 1223637237262778L;
    private static final Logger LOG = Logger.getLogger(AtomAtomMapping.class.getName());
    private final IAtomContainer query;
    private final IAtomContainer target;
    private final HashBiMap<IAtom, IAtom> mapping;
    private final Map<Integer, Integer> mappingIndex;

    /**
     *
     * @param query source molecule
     * @param target target molecule
     */
    public AtomAtomMapping(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        this.mapping = HashBiMap.create(new HashMap<IAtom, IAtom>());
        this.mappingIndex = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final AtomAtomMapping other = (AtomAtomMapping) obj;
        if (this.getQuery() != other.getQuery() && (this.getQuery() == null || !this.query.equals(other.query))) {
            return false;
        }
        if (this.getTarget() != other.getTarget() && (this.getTarget() == null || !this.target.equals(other.target))) {
            return false;
        }
        if (this.mapping != other.getMappingsByAtoms() && (this.mapping == null || !this.mapping.equals(other.mapping))) {
            return false;
        }
        return this.mappingIndex == other.getMappingsByIndex() || (this.mappingIndex != null && this.mappingIndex.equals(other.mappingIndex));
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 67 * hash + (this.getQuery() != null ? this.getQuery().hashCode() : 0);
        hash = 67 * hash + (this.getTarget() != null ? this.getTarget().hashCode() : 0);
        hash = 67 * hash + (this.mapping != null ? this.getMappingsByAtoms().hashCode() : 0);
        hash = 67 * hash + (this.mappingIndex != null ? this.getMappingsByIndex().hashCode() : 0);
        return hash;
    }

    /**
     *
     * @param atom1
     * @param atom2
     */
    public synchronized void put(IAtom atom1, IAtom atom2) {
        mapping.put(atom1, atom2);
        mappingIndex.put(getQuery().getAtomNumber(atom1), getTarget().getAtomNumber(atom2));
    }

    /**
     * Returns String.
     *
     * @return string
     */
    @Override
    public synchronized String toString() {
        String s = "[";
        for (IAtom key : mapping.keySet()) {
            int keyIndex = getQuery().getAtomNumber(key);
            int valueIndex = getTarget().getAtomNumber(mapping.get(key));
            s += keyIndex + ":" + valueIndex + "|";
        }
        return s + "]";
    }

    /**
     * Returns true if 'query' is not isomorphic of 'target'.
     *
     * @return true if 'query' is not isomorphic of 'target'
     */
    public synchronized boolean isEmpty() {
        return mapping.isEmpty();
    }

    /**
     *
     * Clear mappings
     */
    public synchronized void clear() {
        mapping.clear();
        mapping.clear();
    }

    /**
     *
     * Returns mapping size.
     *
     * @return mapping size
     */
    public synchronized int getCount() {
        return mapping.isEmpty() ? 0 : mapping.size();
    }

    /**
     * Returns atom-atom mappings
     *
     * @return atom-atom mappings
     */
    public synchronized Map<IAtom, IAtom> getMappingsByAtoms() {
        return Collections.unmodifiableMap(new HashMap<>(mapping));
    }

    /**
     * Returns atom-atom index mappings
     *
     * @return atom-atom index mappings
     */
    public synchronized Map<Integer, Integer> getMappingsByIndex() {
        return Collections.unmodifiableSortedMap(new TreeMap<>(mappingIndex));
    }

    /**
     * Returns atom index of the given atom in the query molecule
     *
     * @param atom
     * @return
     */
    public synchronized int getQueryIndex(IAtom atom) {
        return getQuery().getAtomNumber(atom);
    }

    /**
     * Returns atom index of the given atom in the target molecule
     *
     * @param atom
     * @return
     */
    public synchronized int getTargetIndex(IAtom atom) {
        return getTarget().getAtomNumber(atom);
    }

    /**
     * Returns query molecule
     *
     * @return the query
     */
    public synchronized IAtomContainer getQuery() {
        return query;
    }

    /**
     * Returns target molecule
     *
     * @return the target
     */
    public synchronized IAtomContainer getTarget() {
        return target;
    }

    /**
     * Returns common mapped fragment in the query molecule.
     *
     * @return common mapped fragment in the query molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainer getMapCommonFragmentOnQuery() throws CloneNotSupportedException {
        IAtomContainer ac = getQuery().clone();
        List<IAtom> uniqueAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getQuery().atoms()) {
            if (!mapping.containsKey(atom)) {
                uniqueAtoms.add(ac.getAtom(getQueryIndex(atom)));
            }
        }

        for (IAtom atom : uniqueAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }

        return ac;
    }

    /**
     * Returns common mapped fragment in the target molecule.
     *
     * @return common mapped fragment in the target molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainer getMapCommonFragmentOnTarget() throws CloneNotSupportedException {

        IAtomContainer ac = getTarget().clone();
        List<IAtom> uniqueAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getTarget().atoms()) {
            if (!mapping.containsValue(atom)) {
                uniqueAtoms.add(ac.getAtom(getTargetIndex(atom)));
            }
        }

        for (IAtom atom : uniqueAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }
        return ac;
    }

    /**
     * Returns common mapped fragment in the query molecule.
     *
     * @return common mapped fragment in the query molecule
     * @throws CloneNotSupportedException
     */
    public synchronized IAtomContainer getCommonFragment() throws CloneNotSupportedException {
        IAtomContainer ac = getQuery().clone();
        List<IAtom> uniqueAtoms = Collections.synchronizedList(new ArrayList<IAtom>());
        for (IAtom atom : getQuery().atoms()) {
            if (!mapping.containsKey(atom)) {
                uniqueAtoms.add(ac.getAtom(getQueryIndex(atom)));
            }
        }

        /*
         Remove bond(s) from the query molecule if they are not present in the target.
         As we are mapping/projecting atoms, it might happen that a bond may or maynot 
         exist between atoms.
         */
        for (IBond bond : getQuery().bonds()) {
            IAtom atom1ForBondInTarget = mapping.get(bond.getAtom(0));
            IAtom atom2ForBondInTarget = mapping.get(bond.getAtom(1));
            IBond bondInTarget = getTarget().getBond(atom1ForBondInTarget, atom2ForBondInTarget);
            if (bondInTarget == null) {
                IAtom atom1InCommonContainer = ac.getAtom(getQueryIndex(bond.getAtom(0)));
                IAtom atom2InCommonContainer = ac.getAtom(getQueryIndex(bond.getAtom(1)));
                ac.removeBond(ac.getBond(atom1InCommonContainer, atom2InCommonContainer));
            }
        }

        for (IAtom atom : uniqueAtoms) {
            ac.removeAtomAndConnectedElectronContainers(atom);
        }

        /*
         Get canonicalised by fixing hydrogens 
         o/p i.e. atom type + CDKHydrogenManipulator
         */
        try {
            CDKHydrogenAdder.getInstance(ac.getBuilder()).addImplicitHydrogens(ac);
        } catch (CDKException ex) {
            Logger.getLogger(AtomAtomMapping.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        } catch (CDKException ex) {
            Logger.getLogger(AtomAtomMapping.class.getName()).log(Level.SEVERE, null, ex);
        }

        return ac;
    }

    /**
     * Returns Maximum Common Fragment between Query and Target as SMILES
     *
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public synchronized String getCommonFragmentAsSMILES() throws CloneNotSupportedException, CDKException {
        SmilesGenerator aromatic = SmilesGenerator.unique().withAtomClasses();
        return aromatic.create(getCommonFragment());
    }
}
