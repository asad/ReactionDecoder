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

import static com.google.common.collect.Maps.newHashMapWithExpectedSize;
import java.util.Map;

import uk.ac.ebi.centres.Descriptor;
import static uk.ac.ebi.centres.Descriptor.Type.ASYMMETRIC;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.priority.AbstractPriorityRule;
import uk.ac.ebi.centres.priority.access.DescriptorAccessor;

/**
 * A configurable descriptor rule that allows ranking of ligands based on their
 * descriptors. The type of descriptor can be provided by a {@link
 * DescriptorAccessor}. The rule type will default to {@link
 * Descriptor.Type#ASYMMETRIC} but can also be configured for R>S rules. The
 * ranking is achieved by assigning each descriptor a value 1 .. n. If a given
 * descriptor is not found in the ranking is assign rank 0.
 *
 * @author John May
 * @param <A>
 */
public class DescriptorRule<A> extends AbstractPriorityRule<A> {

    private DescriptorAccessor<A> accessor;
    private Map<Descriptor, Integer> ranking;

    /**
     * Create an {@link Descriptor.Type#ASYMMETRIC} rule with a provided
     * accessor and given ligand order. Ligand order preceedes that the higher
     * index number the higher the priority.
     *
     * @param accessor a {@link DescriptorAccessor} for a descriptor label
     * @param ordering
     * @param descriptors ranking of descriptors low .. high priority
     */
    public DescriptorRule(DescriptorAccessor<A> accessor,
            Type ordering,
            Descriptor... descriptors) {
        this(ASYMMETRIC, ordering, accessor, descriptors);
    }

    /**
     * Create an rule with a provided rule type, accessor and given ligand
     * order. Ligand order preceedes that the higher index number the higher the
     * priority.
     *
     * @param type the type of priority rule
     * @param ordering
     * @param accessor a {@link DescriptorAccessor} for a descriptor label
     * @param descriptors ranking of descriptors low .. high priority
     */
    public DescriptorRule(Descriptor.Type type,
            Type ordering,
            DescriptorAccessor<A> accessor,
            Descriptor... descriptors) {
        super(type, ordering);
        this.accessor = accessor;

        ranking = newHashMapWithExpectedSize(descriptors.length);

        for (int i = 0; i < descriptors.length; i++) {
            ranking.put(descriptors[i], i + 1);
        }

    }

    /**
     * Access the rank using the accessor and the map.
     *
     * @param ligand the ligand which to access the rank of it's descriptor
     *
     * @return an integer ranking (higher number=higher priority), will default
     * to 0.
     */
    private int getRank(Ligand<A> ligand) {
        Descriptor descriptor = accessor.getDescriptor(ligand);
        Integer rank = ranking.get(descriptor);
        return rank == null ? 0 : rank;
    }

    /**
     * Compares ligands on the rank of their descriptors.
     *
     * @inheritDoc
     */
    @Override
    public int compare(Ligand<A> o1, Ligand<A> o2) {
        return getRank(o1) - getRank(o2);
    }
}
