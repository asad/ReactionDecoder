/* 
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.util.Comparator;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;

 /*
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Map2ValueComparator implements Comparator<Map<IAtom, IAtom>> {

    private final SortOrder sortOrder;

    public Map2ValueComparator(SortOrder sortOrder) {
        this.sortOrder = sortOrder;
    }

    /**
     *
     * @param object1
     * @param object2
     * @return
     */
    @Override
    public int compare(Map<IAtom, IAtom> object1, Map<IAtom, IAtom> object2) {
        int size1 = object1.size();
        int size2 = object2.size();
        int compare = Integer.signum(new Integer(size1).compareTo(size2));

        if (sortOrder == SortOrder.ASCENDING) {
            return compare;
        } else {
            return compare * (-1);
        }
        //return size2 - size1;  assumes you want biggest to smallest;
    }
}
