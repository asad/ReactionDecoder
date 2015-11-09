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
import java.util.Iterator;
import java.util.List;
import uk.ac.ebi.centres.Comparison;
import uk.ac.ebi.centres.Descriptor;
import static uk.ac.ebi.centres.Descriptor.Type.ASYMMETRIC;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.LigandComparison;
import uk.ac.ebi.centres.LigandSorter;
import uk.ac.ebi.centres.Priority;
import uk.ac.ebi.centres.PriorityRule;

/**
 * An abstract comparator that provides construction of the {@link Comparison} wrapper allowing subclasses to focus on
 * the actual comparison of ligands.
 *
 * @author John May
 * @param <A>
 */
public abstract class AbstractPriorityRule<A>
        implements PriorityRule<A> {

    private LigandSorter<A> sorter;
    private boolean halted = FALSE;
    /**
     * The type is store here and appended with the {@link
     * #compareLigands(uk.ac.ebi.centres.Ligand, uk.ac.ebi.centres.Ligand)}
     */
    private final Descriptor.Type reflection;
    private final Type ordering;

    /**
     * Default constructor creates an {@link Descriptor.Type#ASYMMETRIC} comparator.
     * @param ordering
     */
    public AbstractPriorityRule(Type ordering) {
        this(ASYMMETRIC, ordering);
    }

    /**
     * Constructor creates a comparator with the specified type.
     * @param ordering
     */
    public AbstractPriorityRule(Descriptor.Type reflection, Type ordering) {
        this.ordering = ordering;
        this.reflection = reflection;

    }

    @Override
    public void setHalt(boolean halt) {
        this.halted = halt;
    }

    /**
     *
     * @param o1
     * @param o2
     * @return
     */
    public int recursiveCompare(Ligand<A> o1, Ligand<A> o2) {

        int value = compare(o1, o2);
        return value != 0 || halted ? value
                : compare(o1.getLigands(), o2.getLigands());

    }

    /**
     * @inheritDoc
     */
    @Override
    public Comparison compareLigands(Ligand<A> o1, Ligand<A> o2) {
        return new LigandComparison(recursiveCompare(o1, o2), reflection);
    }

    /**
     * @inheritDoc
     */
    public void setSorter(LigandSorter<A> sorter) {
        this.sorter = sorter;
    }

    /**
     * Access the ligand sorter, if the sorter is null a default insertion sorter ({@link InsertionSorter}) is created
     * using 'this; rule as the comparator.
     *
     * @return a set ligand sorter or a newly created insertion sorter
     */
    public LigandSorter<A> getSorter() {
        if (sorter == null) {
            sorter = new InsertionSorter<>(this);
        }
        return sorter;
    }

    /**
     * Uses the injected ligand sorter to order the ligands.
     *
     * @param ligands the ligands that are to be sorted
     *
     * @return whether the ligands are unique
     */
    public Priority prioritise(List<Ligand<A>> ligands) {
        return getSorter().prioritise(ligands);
    }

    /**
     * Compares two lists of ligands. The ligands are first sorted and then iteratively compared by the sub-class
     * comparator. If no different is found whilst iterating through the list the larger of the two lists wins or a tie
     * is determined.
     *
     * @param first first list of ligands
     * @param second second list of ligands
     *
     * @return the value of the comparison
     */
    public int compare(List<Ligand<A>> first, List<Ligand<A>> second) {

        if (halted) {
            return 0;
        }

        // prioritise the ligands, unique isn't required
        prioritise(first);
        prioritise(second);

        // the iterators allow us iterate over the list
        Iterator<Ligand<A>> firstIt = first.iterator();
        Iterator<Ligand<A>> secondIt = second.iterator();

        // compare each element - at the first difference that ligand
        // has priority
        while (firstIt.hasNext() && secondIt.hasNext()) {
            Ligand<A> firstLigand = firstIt.next();
            Ligand<A> secondLigand = secondIt.next();
            int value = compare(firstLigand, secondLigand);
            if (value != 0) {
                return value;
            }
        }

        // no difference found yet, check for different size
        int sizediff = first.size() - second.size();

        if (sizediff != 0) {
            return sizediff;
        }

        // reiterate with recursive compare
        firstIt = first.iterator();
        secondIt = second.iterator();

        // compare each element - at the first difference that ligand
        // has priority
        while (firstIt.hasNext() && secondIt.hasNext()) {
            Ligand<A> firstLigand = firstIt.next();
            Ligand<A> secondLigand = secondIt.next();
            int value = recursiveCompare(firstLigand, secondLigand);
            if (value != 0) {
                return value;
            }
        }

        return 0;
    }

    /**
     *
     * @return
     */
    public boolean isHalted() {
        return halted;
    }

    /**
     * @inheritDoc
     */
    @Override
    public Descriptor.Type getType() {
        return reflection;
    }

    /**
     * Indicates whether the rule is conditional etc.
     *
     * @return
     */
    public Type getRuleType() {
        return ordering;
    }
}
