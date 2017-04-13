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
package uk.ac.ebi.centres;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;

/**
 * Holds some properties that are determined when sorting/prioritising ligands.
 *
 * @author John May
 */
public class Priority {
    private static final Logger LOG = getLogger(Priority.class.getName());

    private Boolean unique;
    private Descriptor.Type type;
    private Set<Set<Integer>> duplicates;

    /**
     *
     * @param unique
     * @param type
     */
    public Priority(Boolean unique, Descriptor.Type type) {
        this.unique = unique;
        this.type = type;
    }

    /**
     *
     * @param unique
     * @param type
     * @param duplicates
     */
    public Priority(Boolean unique, Descriptor.Type type, Set<Set<Integer>> duplicates) {
        this.unique = unique;
        this.type = type;
        this.duplicates = duplicates;
    }

    /**
     * Indicates whether the ligands were unique (i.e. could be ordered)
     *
     * @return whether the ligands were unique
     */
    public Boolean isUnique() {
        return unique;
    }

    /**
     * Indicates the descriptor type used to. This allows methods that represent pseudo-asymmetric molecules to indicate
     * that the centre is pseudo-asymmetric.
     *
     * @return The type of the descriptor that should be assigned
     */
    public Descriptor.Type getType() {
        return type;
    }

    /**
     *
     * @param <A>
     * @param ligands
     * @return
     */
    public <A> List<List<Ligand<A>>> createBins(List<Ligand<A>> ligands) {
        if (duplicates == null) {
            throw new IllegalArgumentException("No duplicates stored at time of sort!");
        }

        List<List<Ligand<A>>> bins = new ArrayList<>(ligands.size());

        // now need to place in bins
        for (Ligand<A> ligand : ligands) {
            List<Ligand<A>> bin = new ArrayList<>();
            bin.add(ligands.get(0));
            bins.add(bin);
        }

        Set<Integer> removed = new HashSet<>();
        // and compact (could be doing something wrong
        for (Set<Integer> pair : duplicates) {
            Iterator<Integer> it = pair.iterator();
            int i = it.next();
            int j = it.next();
            if (!removed.contains(i) || !removed.contains(j)) {
                bins.get(i).addAll(bins.get(j));
                removed.add(j);
            }
        }
        for (Integer remove : removed) {
            bins.remove(remove);
        }

        return bins;

    }
}
