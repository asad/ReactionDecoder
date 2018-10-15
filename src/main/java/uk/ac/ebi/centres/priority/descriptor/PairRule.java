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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.priority.descriptor;

import static com.google.common.collect.Collections2.permutations;
import static com.google.common.collect.Lists.newLinkedList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NavigableSet;
import java.util.Queue;
import java.util.TreeSet;

import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.Priority;
import static uk.ac.ebi.centres.PriorityRule.Type.GEOMETRICAL;
import static uk.ac.ebi.centres.descriptor.General.NONE;
import static uk.ac.ebi.centres.descriptor.General.UNKNOWN;
import static uk.ac.ebi.centres.descriptor.General.UNSPECIFIED;
import uk.ac.ebi.centres.priority.AbstractPriorityRule;
import uk.ac.ebi.centres.priority.access.DescriptorAccessor;

/**
 * A descriptor pair rule. This rule defines that like descriptor pairs have
 * priority over unlike descriptor pairs.
 *
 * @author John May
 * @param <A>
 */
public class PairRule<A>
        extends AbstractPriorityRule<A> {

    /**
     * Ugly piece of code to generate permutation of the given ligand groups.
     * There may be a much better way to do this. This method converts lists
     * with duplicates into all possible combinations. A, {B1, B2}, C would
     * permutate to A, B1, B2, C and A, B2, B1, C.
     * <p/>
     * This method was adapted from http://goo.gl/s6R7E
     *
     * @see <a href="http://goo.gl/s6R7E">http://www.daniweb.com/</a>
     */
    private static <T> List<List<T>> permutate(List<List<T>> uncombinedList) {
        List<List<T>> list = new ArrayList<>();

        // permeate the sublist
        for (List sublist : uncombinedList) {
            if (sublist.size() > 1) {
                Collection<List> tmp;
                tmp = permutations(sublist);
                sublist.clear();
                sublist.addAll(tmp);
            }
        }

        int index[] = new int[uncombinedList.size()];
        int combinations = combinations(uncombinedList) - 1;
        // Initialize index
        for (int i = 0; i < index.length; i++) {
            index[i] = 0;
        }
        // First combination is always valid
        List<T> combination = new ArrayList<>();
        for (int m = 0; m < index.length; m++) {
            Object o = uncombinedList.get(m).get(index[m]);
            if (o instanceof Collection) {
                combination.addAll((Collection) o);
            } else {
                combination.add((T) o);
            }
        }
        list.add(combination);

        for (int k = 0; k < combinations; k++) {
            combination = new ArrayList<>();
            boolean found = false;
            // We Use reverse order
            for (int l = index.length - 1; l >= 0 && found == false; l--) {
                int currentListSize = uncombinedList.get(l).size();
                if (index[l] < currentListSize - 1) {
                    index[l] += 1;
                    found = true;
                } else {
                    // Overflow
                    index[l] = 0;
                }
            }
            for (int m = 0; m < index.length; m++) {
                Object o = uncombinedList.get(m).get(index[m]);
                if (o instanceof Collection) {
                    combination.addAll((Collection) o);
                } else {
                    combination.add((T) o);
                }
            }
            list.add(combination);
        }
        return list;
    }

    private static <T> int combinations(List<List<T>> list) {
        int count = 1;
        for (List<T> current : list) {
            count *= current.size();
        }
        return count;
    }

    private final DescriptorAccessor<A> accessor;

    /**
     * Construct a pair rule with a given accessor. The accessor is used to
     * retrieve the required descriptor (e.g. primary, auxiliary, arc, etc.) on
     * each ligand.
     *
     * @param accessor access to a descriptor on a ligand
     */
    public PairRule(DescriptorAccessor<A> accessor) {
        super(GEOMETRICAL);
        this.accessor = accessor;
    }

    /**
     * Generates a set of descriptor lists that maintain the like/unlike pairing
     * whilst descriptors are added. The set is navigable and maintains priority
     * ordering when multiple lists are present. This method is a wrapper for
     * adding the seeding ligand to the queue.
     *
     * @param ligand the ligand on which to generate the descriptor lists for
     *
     * @return navigable set of descriptor lists
     */
    protected NavigableSet<DescriptorList> generate(Ligand<A> ligand) {
        // would be good to give an expected size
        Queue<Ligand<A>> queue = newLinkedList();
        queue.add(ligand);
        return generate(queue);
    }

    /**
     *
     * @param o1
     * @param o2
     * @return
     */
    @Override
    public int recursiveCompare(Ligand<A> o1, Ligand<A> o2) {
        // can't/don't need to do recursive on the pair rule
        return compare(o1, o2);
    }

    /**
     * Generates a set of descriptor lists that maintain the like/unlike pairing
     * whilst descriptors are added. The set is navigable and maintains priority
     * ordering when multiple lists are present.
     *
     * @param queue a queue of ligands for which to get descriptors and expand
     *
     * @return navigable set of descriptor lists
     */
    protected NavigableSet<DescriptorList> generate(Queue<Ligand<A>> queue) {

        NavigableSet<DescriptorList> lists = new TreeSet<>();

        // create a descriptor list with given exclusions
        DescriptorList descriptors = new DescriptorList(null, NONE, UNSPECIFIED, UNKNOWN);

        while (!queue.isEmpty()) {

            Ligand<A> ligand = queue.poll();
            descriptors.add(accessor.getDescriptor(ligand));

            List<Ligand<A>> ligands = ligand.getLigands();
            Priority priority = prioritise(ligands);
            if (priority.isUnique()) {

                // unique
                for (Ligand<A> child : ligands) {
                    queue.add(child);
                }

            } else {
                // non unique need to subdivide and combine
                for (List<Ligand<A>> combinated : permutate(getSorter().getGroups(ligands))) {

                    Queue<Ligand<A>> subqueue = new LinkedList<>(queue);
                    subqueue.addAll(combinated);

                    // add to current descriptor list
                    lists.addAll(descriptors.append(generate(subqueue)));

                }

                // queue was copied and delegated so we clear this instance
                queue.clear();

            }

        }

        if (lists.isEmpty()) {
            lists.add(descriptors);
        }

        return lists;

    }

    /**
     * Compares the two ligands based on their generated descriptor pairs.
     *
     * @param o1 first ligand
     * @param o2 second ligand
     *
     * @return the value of the comparison
     */
    @Override
    public int compare(Ligand<A> o1, Ligand<A> o2) {

        // produced pair lists are in order
        Iterator<DescriptorList> list1It = generate(o1).iterator();
        Iterator<DescriptorList> list2It = generate(o2).iterator();

        while (list1It.hasNext() && list2It.hasNext()) {
            int value = list1It.next().compareTo(list2It.next());
            if (value != 0) {
                return value;
            }
        }

        // there may be a different is list size but normally you'd have a
        // constitutional rule (which would find this) before this pairing rule
        // we don't go to the next level on this rule. We've already
        // exhaustively create pair lists (generate) for each ligand.
        return 0;
    }

}
