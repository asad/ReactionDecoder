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
package uk.ac.ebi.centres;

import java.util.List;

/**
 * An injectable sorter for ligands.
 *
 * @author John May
 * @param <A>
 */
public interface LigandSorter<A> {

    /**
     * Sorts the provided ligands and indicates if all the ligands are different (i.e. unique). The method is named
     * prioritise to emphasise that ligands are sorting in descending order (i.e. ranked highest to lowest). The default
     * {@link java.util.Collections#sort(java.util.List)} performs an ascending order sort.
     *
     * @param ligands the ligands that will be sorted
     *
     * @return whether the ligands are all different and which type of rule was used
     */
    public Priority prioritise(List<Ligand<A>> ligands);

    /**
     *
     * @param sorted
     * @return
     */
    public List<List<Ligand<A>>> getGroups(List<Ligand<A>> sorted);
}
