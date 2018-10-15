/* Copyright (C) 2011  Gilleain Torrance <gilleain.torrance@gmail.com>
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package uk.ac.ebi.reactionblast.tools.labelling;

import static java.lang.Integer.valueOf;
import static java.lang.System.out;
import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.AtomContainerSet;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class AbstractReactionLabeller {

    /**
     * A nasty hack necessary to get around a bug in the CDK
     */
    private final boolean fixAtomMappingCastType = false;

    private void fixAtomMapping(IAtomContainer canonicalForm) {
        for (IAtom a : canonicalForm.atoms()) {
            String v = (String) a.getProperty(ATOM_ATOM_MAPPING);
            if (v != null) {
                a.setProperty(ATOM_ATOM_MAPPING, valueOf(v));
            }
        }
    }

    /**
     *
     * @param original
     * @param clone
     * @param permutationMap
     * @param atomAtom
     */
    protected void atomAtomMap(IAtomContainerSet original,
            IAtomContainerSet clone, Map<IAtomContainer, int[]> permutationMap,
            Map<IAtom, IAtom> atomAtom) {
        for (int i = 0; i < original.getAtomContainerCount(); ++i) {
            IAtomContainer mol = original.getAtomContainer(i);
            IAtomContainer mol2 = clone.getAtomContainer(i);
            int[] permutation = permutationMap.get(mol2);
            for (int j = 0; j < mol.getAtomCount(); ++j) {
                atomAtom.put(mol.getAtom(j), mol2.getAtom(permutation[j]));
            }
        }
    }

    /**
     *
     * @param reaction
     * @param clone
     * @param permutationMap
     * @return
     */
    protected Map<IAtom, IAtom> atomAtomMap(
            IReaction reaction, IReaction clone,
            Map<IAtomContainer, int[]> permutationMap) {
        // create a Map of corresponding atoms for molecules 
        // (key: original Atom, value: clone Atom)
        Map<IAtom, IAtom> atomAtom = new Hashtable<>();
        IAtomContainerSet reactants = reaction.getReactants();
        IAtomContainerSet clonedReactants = clone.getReactants();
        atomAtomMap(reactants, clonedReactants, permutationMap, atomAtom);

        IAtomContainerSet products = reaction.getProducts();
        IAtomContainerSet clonedProducts = clone.getProducts();
        atomAtomMap(products, clonedProducts, permutationMap, atomAtom);

        return atomAtom;
    }

    /**
     *
     * @param reaction
     * @param atomAtomMap
     * @return
     */
    protected List<IMapping> cloneMappings(
            IReaction reaction, Map<IAtom, IAtom> atomAtomMap) {
        // clone the mappings
        int numberOfMappings = reaction.getMappingCount();
        List<IMapping> map = new ArrayList<>();
        for (int mappingIndex = 0; mappingIndex < numberOfMappings; mappingIndex++) {
            IMapping mapping = reaction.getMapping(mappingIndex);
            map.add(cloneMapping(mapping, atomAtomMap));
        }
        return map;
    }

    /**
     *
     * @param mapping
     * @param atomAtomMap
     * @return
     */
    protected IMapping cloneMapping(
            IMapping mapping, Map<IAtom, IAtom> atomAtomMap) {
        IChemObject keyChemObj0 = mapping.getChemObject(0);
        IChemObject keyChemObj1 = mapping.getChemObject(1);
        IChemObject co0 = atomAtomMap.get(keyChemObj0);
        IChemObject co1 = atomAtomMap.get(keyChemObj1);

        // THIS IS STUPID : BLAME THE IDIOT WHO FAILED TO PUT SET METHODS IN
        // IMAPPING (OR IREACTION, FOR THAT MATTER)
        if (co0 == null) {
            co0 = new ChemObject();
        }
        if (co1 == null) {
            co1 = new ChemObject();
        }
        return new Mapping(co0, co1);
    }

    /**
     *
     * @param reaction
     * @return
     */
    protected Map<IChemObject, Integer> makeIndexMap(IReaction reaction) {
        Map<IChemObject, Integer> indexMap
                = new HashMap<>();
        List<IAtomContainer> all
                = getAllAtomContainers(reaction);
        int globalIndex = 0;
        for (IAtomContainer ac : all) {
            for (IAtom atom : ac.atoms()) {
                indexMap.put(atom, globalIndex);
                globalIndex++;
            }
        }
        return indexMap;
    }

    /**
     * Clone and Sort the mappings based on the order of the first object in the
     * mapping (which is assumed to be the reactant).
     *
     * @param reaction
     * @param copyOfReaction
     * @param permutationMap
     */
    protected void cloneAndSortMappings(
            IReaction reaction, IReaction copyOfReaction,
            Map<IAtomContainer, int[]> permutationMap) {

        Map<IAtom, IAtom> atomAtomMap = atomAtomMap(
                reaction, copyOfReaction, permutationMap);
        List<IMapping> map = cloneMappings(reaction, atomAtomMap);
        sortMappings(copyOfReaction, map);
    }

    /**
     *
     * @param reaction
     * @param map
     */
    protected void sortMappings(IReaction reaction, List<IMapping> map) {
        // make a lookup for the indices of the atoms 
        final Map<IChemObject, Integer> indexMap = makeIndexMap(reaction);
        Comparator<IMapping> mappingSorter = (IMapping o1, IMapping o2) -> {
            IChemObject o10 = o1.getChemObject(0);
            IChemObject o20 = o2.getChemObject(0);
            if (o20 == null || o10 == null) {
                return 0;
            }
            Integer o10i = indexMap.get(o10);
            Integer o20i = indexMap.get(o20);
            if (o10i == null || o20i == null) {
                return 0;
            }
            return o10i.compareTo(o20i);
        };
        sort(map, mappingSorter);
        int mappingIndex = 0;
        for (IMapping mapping : map) {
            IChemObject o0 = mapping.getChemObject(0);
            if (o0 != null) {
                o0.setProperty(ATOM_ATOM_MAPPING, mappingIndex);
            }
            IChemObject o1 = mapping.getChemObject(1);
            if (o1 != null) {
                o1.setProperty(ATOM_ATOM_MAPPING, mappingIndex);
            }
            reaction.addMapping(mapping);
            mappingIndex++;
        }
    }

    /**
     *
     * @param moleculeSet
     * @param labeller
     * @param permutationMap
     * @return
     */
    protected IAtomContainerSet canoniseAtomContainerSet(
            IAtomContainerSet moleculeSet,
            ICanonicalMoleculeLabeller labeller,
            Map<IAtomContainer, int[]> permutationMap) {
        IAtomContainerSet canonicalMolecules = new AtomContainerSet();
        for (IAtomContainer atomContainer : moleculeSet.atomContainers()) {
            IAtomContainer canonicalForm
                    = labeller.getCanonicalMolecule(atomContainer);
            if (fixAtomMappingCastType) {
                fixAtomMapping(canonicalForm);
            }
            canonicalForm.setID(atomContainer.getID());
            permutationMap.put(
                    canonicalForm,
                    labeller.getCanonicalPermutation(atomContainer));
            canonicalMolecules.addAtomContainer(canonicalForm);
        }
        return canonicalMolecules;
    }

    /**
     *
     * @param reaction
     * @param labeller
     * @return
     */
    public IReaction labelReaction(
            IReaction reaction, ICanonicalMoleculeLabeller labeller) {
        out.println("labelling " + reaction.getID());
        IReaction canonReaction = new Reaction();

        Map<IAtomContainer, int[]> permutationMap = new HashMap<>();

        canonReaction.setProducts(
                canoniseAtomContainerSet(
                        reaction.getProducts(), labeller, permutationMap));
        canonReaction.setReactants(
                canoniseAtomContainerSet(
                        reaction.getReactants(), labeller, permutationMap));

        cloneAndSortMappings(reaction, canonReaction, permutationMap);
        canonReaction.setID(reaction.getID());
        return canonReaction;
    }
}
