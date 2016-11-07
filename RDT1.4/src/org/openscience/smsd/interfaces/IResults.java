/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.interfaces;

import java.util.List;
import org.openscience.smsd.AtomAtomMapping;

/**
 * Interface that holds basic core interface for all MCS algorithm.  
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public interface IResults {

    /**
     * Returns all plausible mappings between query and target molecules. Each map in the list has atom-atom equivalence
     * of the mappings between query and target molecule i.e. map.getKey() for the query and map.getValue() for the
     * target molecule
     *
     * @return All possible MCS atom Mappings
     */
    public abstract List<AtomAtomMapping> getAllAtomMapping();

    /**
     * Returns one of the best matches with atoms mapped.
     *
     * @return Best Atom Mapping
     */
    public abstract AtomAtomMapping getFirstAtomMapping();
}
