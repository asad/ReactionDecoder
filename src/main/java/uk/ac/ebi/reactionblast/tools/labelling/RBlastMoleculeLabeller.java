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
package uk.ac.ebi.reactionblast.tools.labelling;

import java.util.ArrayList;
import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.signature.MoleculeSignature;
import signature.AbstractVertexSignature;
import static uk.ac.ebi.reactionblast.tools.labelling.AtomContainerAtomPermutor.permute;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class RBlastMoleculeLabeller implements ICanonicalMoleculeLabeller {

    /**
     *
     * @param container
     * @return
     */
    @Override
    public IAtomContainer getCanonicalMolecule(IAtomContainer container) {
        return permute(
                getCanonicalPermutation(container), container);
    }

    /**
     *
     * @param container
     * @return
     */
    @Override
    public int[] getCanonicalPermutation(IAtomContainer container) {
        int atomCount = container.getAtomCount();
        MoleculeSignature moleculeSignature = new MoleculeSignature(container);
        String canonicalSignatureString = null;
        AbstractVertexSignature canonicalSignature = null;
        Map<String, SortableOrbit> orbitMap = new HashMap<>();
        for (int atomIndex = 0; atomIndex < atomCount; atomIndex++) {
            AbstractVertexSignature signatureForI = moleculeSignature.signatureForVertex(atomIndex);
            String signatureStringForVertexI = signatureForI.toCanonicalString();
            if (canonicalSignatureString == null
                    || canonicalSignatureString.compareTo(signatureStringForVertexI) < 0) {
                canonicalSignatureString = signatureStringForVertexI;
                canonicalSignature = signatureForI;
            }
            if (orbitMap.containsKey(signatureStringForVertexI)) {
                orbitMap.get(signatureStringForVertexI).members.add(atomIndex);
            } else {
                String elementSymbol = container.getAtom(atomIndex).getSymbol();
                SortableOrbit orbit
                        = new SortableOrbit(elementSymbol, signatureStringForVertexI);
                orbit.members.add(atomIndex);
                orbitMap.put(signatureStringForVertexI, orbit);
            }
        }
        int[] canonicalLabels = canonicalSignature == null ? new int[0] : canonicalSignature.getCanonicalLabelling(atomCount);
        List<SortableOrbit> orbits = new ArrayList<>();
        orbitMap.values().stream().map((orbit) -> {
            orbit.sortMembers(canonicalLabels);
            return orbit;
        }).forEach((orbit) -> {
            orbits.add(orbit);
        });
        sort(orbits);
        int[] permutation = new int[atomCount];
        int newIndex = 0;
        for (SortableOrbit orbit : orbits) {
            for (int index : orbit.members) {
                permutation[index] = newIndex;
                newIndex++;
            }
        }
        return permutation;
    }

    private class OrbitLabelComparator implements Comparator<Integer> {

        public int[] canonicalLabels;

        OrbitLabelComparator(int[] canonicalLabels) {
            this.canonicalLabels = canonicalLabels;
        }

        @Override
        public int compare(Integer label0, Integer label1) {
            if (canonicalLabels[label0] < canonicalLabels[label1]) {
                return -1;
            } else if (canonicalLabels[label0] > canonicalLabels[label1]) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    private class SortableOrbit implements Comparable<SortableOrbit> {

        public String elementLabel;
        public String signatureLabel;
        public List<Integer> members;

        SortableOrbit(String elementLabel, String label) {
            this.signatureLabel = label;
            this.elementLabel = elementLabel;
            this.members = new ArrayList<>();
        }

        public void sortMembers(int[] canonicalLabels) {
            sort(members, new OrbitLabelComparator(canonicalLabels));
        }

        @Override
        public int compareTo(SortableOrbit other) {
            int elementCompare = elementLabel.compareTo(other.elementLabel);
            if (elementCompare != 0) {
                return elementCompare;
            } else {
                return -signatureLabel.compareTo(other.signatureLabel);
            }
        }
    }
}
