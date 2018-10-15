/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.ebi.centres.priority;

import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import uk.ac.ebi.centres.Comparison;
import uk.ac.ebi.centres.Descriptor;
import static uk.ac.ebi.centres.Descriptor.Type.NON_STEREOGENIC;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.LigandSorter;
import uk.ac.ebi.centres.Priority;
import uk.ac.ebi.centres.PriorityRule;

/**
 * A simple insertion sort for ligands. The number of ligands is not likely to
 * be very larger as such doing a merge sort would have little benefit.
 *
 * @author John May
 * @param <A>
 */
public class InsertionSorter<A> implements LigandSorter<A> {

    private final List<PriorityRule<A>> rules = new ArrayList<>(5);

    /**
     *
     * @param comparator
     */
    public InsertionSorter(PriorityRule<A> comparator) {
        this.rules.add(comparator);
    }

    /**
     *
     * @param comparators
     * @param restrict
     */
    public InsertionSorter(List<PriorityRule<A>> comparators, PriorityRule.Type restrict) {
        for (PriorityRule<A> rule : comparators) {
            if (rule.getRuleType().equals(restrict)) {
                this.rules.add(rule);
            }
        }
    }

    /**
     *
     * @param comparators
     */
    public InsertionSorter(List<PriorityRule<A>> comparators) {
        rules.addAll(comparators);
    }

    /**
     * Sorts in descending order and indicates whether all elements are unique
     * and the type of descriptor used.
     *
     * @inheritDoc
     */
    @Override
    public Priority prioritise(List<Ligand<A>> ligands) {

        Boolean unique = TRUE;
        Descriptor.Type type = NON_STEREOGENIC;
//        Set<Set<Integer>> duplicates = null;

        for (int i = 0; i < ligands.size(); i++) {
            for (int j = i; j > 0; j--) {

                Comparison comparison = compareLigands(ligands.get(j - 1),
                        ligands.get(j));

                type = comparison.getType().ordinal() > type.ordinal()
                        ? comparison.getType() : type;

                if (comparison.getOrder() < 0) {
                    swap(ligands, j, j - 1);
                } else {
                    if (comparison.getOrder() == 0) {
                        unique = FALSE;
                    }
                    break;
                }

            }
        }

        return new Priority(unique, type);

    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    public Comparison compareLigands(Ligand<A> first, Ligand<A> second) {
        for (PriorityRule<A> rule : rules) {
            Comparison comparison = rule.compareLigands(first, second);
            if (comparison.getOrder() != 0) {
                return comparison;
            }
        }
        return new Comparison() {
            @Override
            public Integer getOrder() {
                return 0;
            }

            @Override
            public Descriptor.Type getType() {
                return NON_STEREOGENIC;
            }
        };
    }

    /**
     *
     * @param list
     * @param i
     * @param j
     */
    public void swap(List list, int i, int j) {
        Object tmp = list.get(i);
        list.set(i, list.get(j));
        list.set(j, tmp);
    }

    /**
     *
     * @param sorted
     * @return
     */
    @Override
    public List<List<Ligand<A>>> getGroups(List<Ligand<A>> sorted) {

        // would be nice to have this integrated whilst sorting - may provide a small speed increase
        // but as most of our lists are small we take use ugly sort then group approach
        LinkedList<List<Ligand<A>>> groups = new LinkedList<>();

        sorted.stream().map((ligand) -> {
            if (groups.isEmpty()
                    || compareLigands(groups.getLast().iterator().next(),
                            ligand).getOrder() != 0) {
                groups.add(new ArrayList<>());
            }
            return ligand;
        }).forEachOrdered((ligand) -> {
            groups.getLast().add(ligand);
        });

        return groups;

    }
}
