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
package com.bioinceptionlabs.centres.priority.descriptor;

import static com.bioinceptionlabs.centres.Descriptor.Type.PSEUDO_ASYMMETRIC;
import static com.bioinceptionlabs.centres.PriorityRule.Type.TOPOGRAPHICAL;
import static com.bioinceptionlabs.centres.descriptor.Tetrahedral.R;
import static com.bioinceptionlabs.centres.descriptor.Tetrahedral.S;
import com.bioinceptionlabs.centres.priority.access.DescriptorAccessor;

/**
 * A rule with prioritises ligands in R configuration over those in S
 * configuration. This rule is pseudo-asymmetric
 *
 * @author John May
 * @param <A>
 */
public class RSRule<A> extends DescriptorRule<A> {

    /**
     *
     * @param accessor
     */
    public RSRule(DescriptorAccessor<A> accessor) {
        super(PSEUDO_ASYMMETRIC, TOPOGRAPHICAL,
                accessor, S, R);
    }
}
