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
package uk.ac.ebi.centres.priority;

import static com.google.common.collect.Lists.newArrayListWithExpectedSize;
import java.util.List;

import uk.ac.ebi.centres.Comparison;
import uk.ac.ebi.centres.Descriptor;
import static uk.ac.ebi.centres.Descriptor.Type.NON_STEREOGENIC;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.LigandComparison;
import uk.ac.ebi.centres.LigandSorter;
import uk.ac.ebi.centres.PriorityRule;
import static uk.ac.ebi.centres.PriorityRule.Type.COMBINED;
import static uk.ac.ebi.centres.PriorityRule.Type.CONSTITUTIONAL;

/**
 * A priority rule made up of other rules. Each sub-rule is used exhaustively on
 * the digraph before the next one is applied.
 *
 * @author John May
 * @param <A>
 */
public final class CombinedRule<A> extends AbstractPriorityRule<A> {

    /**
     * Rule storage
     */
    private final List<PriorityRule<A>> rules;

    /**
     * Default constructor creates a combined rule with no sub-rules.
     */
    public CombinedRule() {
        super(COMBINED);
        rules = newArrayListWithExpectedSize(8);
    }

    /**
     * Creates a combined rule from several provided sub-rules.
     *
     * @param rules the rules to combined
     */
    public CombinedRule(PriorityRule<A>... rules) {
        super(COMBINED);
        this.rules = newArrayListWithExpectedSize(rules.length);
        for (PriorityRule<A> rule : rules) {
            add(rule);
        }
    }

    /**
     * Add a priority rule to the compound rule. This will also set the sorter
     * to that of this combined rule.
     *
     * @param rule a new rule to use
     */
    public void add(PriorityRule<A> rule) {
        if (rule == null) {
            throw new IllegalArgumentException("Provided priority rule was"
                    + "null!");
        }
        rules.add(rule);
        rule.setSorter(createSorter(rules));
    }

    /**
     *
     * @param rules
     * @return
     */
    public LigandSorter<A> createSorter(List<PriorityRule<A>> rules) {
        return new InsertionSorter<>(rules, CONSTITUTIONAL); // restriction should be configurable
    }

    /**
     * Iteratively compares ligands using the given priority rules. Each rule is
     * applied exhaustively. If a difference is found for a rule the comparison
     * is returned along without the rule type.
     *
     * @return
     * @see LigandComparison
     * @see uk.ac.ebi.centres.PriorityRule#getType()
     * @see Descriptor.Type
     */
    @Override
    public int compare(Ligand<A> o1, Ligand<A> o2) {

        // Try using each rule. The rules will expand the search exhaustively
        // to all child ligands
        for (PriorityRule<A> rule : rules) {

            if (isHalted()) {
                return 0;
            }

            // compare expands exhaustively across the whole graph
            int value = rule.recursiveCompare(o1, o2);

            if (value != 0) {
                return value;
            }

        }

        return 0;

    }

    /**
     * Iteratively compares ligands using the given priority rules. Each rule is
     * applied exhaustively. If a difference is found for a rule the comparison
     * is returned along with the rule type
     * {@link uk.ac.ebi.centres.PriorityRule#getType()}
     *
     * @return
     * @see LigandComparison
     * @see uk.ac.ebi.centres.PriorityRule#getType()
     * @see Descriptor.Type
     */
    @Override
    public Comparison compareLigands(Ligand<A> o1, Ligand<A> o2) {

        // Try using each rule. The rules will expand the search exhaustively
        // to all child ligands
        for (PriorityRule<A> rule : rules) {

            if (isHalted()) {
                return new LigandComparison(0, NON_STEREOGENIC);
            }

            // compare expands exhaustively across the whole graph
            int value = rule.recursiveCompare(o1, o2);

            if (value != 0) {
                return new LigandComparison(value, rule.getType());
            }

        }

        // can't really give a rule type here...
        return new LigandComparison(0, NON_STEREOGENIC);
    }

    @Override
    public void setHalt(boolean halt) {
        for (PriorityRule<A> rule : rules) {
            rule.setHalt(halt);
        }
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder("Combined rule:");
        for (PriorityRule<A> rule : rules) {
            builder.append(rule.toString()).append(", ");
        }
        return builder.toString();
    }
}
