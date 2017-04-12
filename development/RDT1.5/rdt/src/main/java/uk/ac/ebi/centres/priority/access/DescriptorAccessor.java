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
package uk.ac.ebi.centres.priority.access;

import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.descriptor.General;

/**
 * Provides access to a descriptor from the provided ligand. The accessor can be injected into priority rules to alter
 * their behaviour by accessing different descriptors (e.g. main, auxiliary or bond descriptor).
 *
 * @param <A> the class type of atom the ligand holds
 *
 * @author John May
 */
public interface DescriptorAccessor<A> {

    /**
     * Access the descriptor for the provided ligand. If no descriptor is set for the ligand this method should return
     * {@link General#NONE}
     *
     * @param ligand the ligand to get the descriptor for
     *
     * @return the descriptor for the ligand
     *
     * @see Descriptor
     * @see uk.ac.ebi.centres.descriptor.Tetrahedral
     * @see uk.ac.ebi.centres.descriptor.General
     * @see uk.ac.ebi.centres.descriptor.Planar
     * @see uk.ac.ebi.centres.descriptor.Trigonal
     */
    public Descriptor getDescriptor(Ligand<A> ligand);
}
