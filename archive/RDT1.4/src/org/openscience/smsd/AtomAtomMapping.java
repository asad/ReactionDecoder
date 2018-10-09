/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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

import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 * Holds atom-atom mappings information between source and target molecules
 *
 *
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class AtomAtomMapping implements Serializable {

    class MyQueryIAtomComp implements Comparator<IAtom> {

        @Override
        public int compare(IAtom o1, IAtom o2) {
            int atomNumber1 = getQuery().indexOf(o1);
            int atomNumber2 = getQuery().indexOf(o2);
            return atomNumber1 - atomNumber2;
        }
    }
    private static final long serialVersionUID = 1223637237262778L;
    private final IAtomContainer query;
    private final IAtomContainer target;
    private final Map<IAtom, IAtom> mapping;
    private final Map<Integer, Integer> mappingIndex;

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
     * @param query source molecule
     * @param target target molecule
     */
    public AtomAtomMapping(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        mapping = new TreeMap<>(new MyQueryIAtomComp());
        this.mappingIndex = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
    }

    /**
     *
     * @param atom1
     * @param atom2
     */
    public synchronized void put(IAtom atom1, IAtom atom2) {
        mapping.put(atom1, atom2);
        mappingIndex.put(getQuery().indexOf(atom1), getTarget().indexOf(atom2));
    }

    /**
     * Returns String with MMP and AAM.
     *
     * @return string
     */
    @Override
    public synchronized String toString() {
        StringBuilder s = new StringBuilder();
        try {
            IReaction reaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
            reaction.addReactant(getQuery().clone(), 1.0);
            reaction.addProduct(getTarget().clone(), 1.0);

            int counter = 1;
            for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                for (IAtom a : ac.atoms()) {
                    IAtom refAtom = getQuery().getAtom(ac.indexOf(a));
                    if (mapping.containsKey(refAtom)) {
                        a.setProperty(ATOM_ATOM_MAPPING, counter);
                        a.setFlag(MAPPED, true);
                        IAtom mappedAtom = mapping.get(refAtom);
                        int mappedAtomIndex = getTarget().indexOf(mappedAtom);
                        IAtom b = reaction.getProducts().getAtomContainer(0).getAtom(mappedAtomIndex);
                        b.setProperty(ATOM_ATOM_MAPPING, counter);
                        b.setFlag(MAPPED, true);
                        IMapping aammapping
                                = reaction.getBuilder().newInstance(IMapping.class, a, b);
                        reaction.addMapping(aammapping);
                    }
                    counter++;
                }
            }

            String createReactionSMILES = "NA";
            try {
                SmilesGenerator withAtomClasses = SmilesGenerator.generic().aromatic().withAtomClasses();
                createReactionSMILES = withAtomClasses.createReactionSMILES(reaction);
            } catch (CDKException ex) {
                Logger.getLogger(AtomAtomMapping.class.getName()).log(Level.SEVERE, null, ex);
            }

            s.append("MMP: ").append(createReactionSMILES).append(", AAM:[");
            for (IAtom firstAtom : mapping.keySet()) {
                int keyIndex = getQuery().indexOf(firstAtom) + 1;
                int valueIndex = getTarget().indexOf(mapping.get(firstAtom)) + 1;
                s.append(keyIndex).append(":").append(valueIndex).append("|");
            }

            s.append("]");

            try {
                s.append(", MCS: ").append(getCommonFragmentAsSMILES());
            } catch (CDKException ex) {
                Logger.getLogger(AtomAtomMapping.class.getName()).log(Level.SEVERE, null, ex);
            }

        } catch (CloneNotSupportedException ex) {
            Logger.getLogger(AtomAtomMapping.class.getName()).log(Level.SEVERE, null, ex);
        }
        return s.toString();
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
        mappingIndex.clear();
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
        return Collections.unmodifiableMap(new LinkedHashMap<>(mapping));
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
        return getQuery().indexOf(atom);
    }

    /**
     * Returns atom index of the given atom in the target molecule
     *
     * @param atom
     * @return
     */
    public synchronized int getTargetIndex(IAtom atom) {
        return getTarget().indexOf(atom);
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

        uniqueAtoms.stream().forEach((atom) -> {
            ac.removeAtom(atom);
        });

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

        uniqueAtoms.stream().forEach((atom) -> {
            ac.removeAtom(atom);
        });
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

        uniqueAtoms.stream().forEach((atom) -> {
            ac.removeAtom(atom);
        });

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

    /*
     * Java method to sort Map in Java by value e.g. HashMap or Hashtable
     * throw NullPointerException if Map contains null values
     * It also sort values even if they are duplicates
     */
    public Map<IAtom, IAtom> sortByValues(Map<IAtom, IAtom> map) {
        List<Map.Entry<IAtom, IAtom>> entries = new LinkedList<>(map.entrySet());
        Collections.sort(entries, (Entry<IAtom, IAtom> o1, Entry<IAtom, IAtom> o2) -> {
            int atomNumber1 = getQuery().indexOf(o1.getKey());
            int atomNumber2 = getQuery().indexOf(o2.getKey());
            return atomNumber1 - atomNumber2;
        });

        //LinkedHashMap will keep the keys in the order they are inserted
        //which is currently sorted on natural ordering
        Map<IAtom, IAtom> sortedMap = new LinkedHashMap<>();

        entries.stream().forEach((entry) -> {
            sortedMap.put(entry.getKey(), entry.getValue());
        });

        return sortedMap;
    }

}
