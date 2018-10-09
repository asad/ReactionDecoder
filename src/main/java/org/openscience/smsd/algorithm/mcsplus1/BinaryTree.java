/* Copyright (C) 2009-2018  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.algorithm.mcsplus1;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * the second part of the program extents the mapping by the McGregor algorithm
 * in case that not all atoms of molecule A and molecule B are mapped by the
 * clique approach
 *
 */
public class BinaryTree {

    /**
     * Creates a new instance of BinaryTree
     */
    public BinaryTree() {
        this.equal = null;
        this.not_equal = null;
    }
    public BinaryTree equal;
    public BinaryTree not_equal;
    int value;

    public void setValue(int val) {
        this.value = val;
    }

    public int getValue() {
        return this.value;
    }

    /* 
     * Remove the nodes recursively 
     *
     *
     * @param cur_struc
     * @return
     */
    protected static int remove_tree_structure(BinaryTree cur_struc) {

        BinaryTree equal_struc = cur_struc.equal;
        BinaryTree not_equal_struc = cur_struc.not_equal;
//      delete(cur_struc);//TODO by ASAD in java here is automatic pointer deleting
        if (equal_struc != null) {
            remove_tree_structure(equal_struc);
        }
        if (not_equal_struc != null) {
            remove_tree_structure(not_equal_struc);
        }

        return 0;
    }
}
