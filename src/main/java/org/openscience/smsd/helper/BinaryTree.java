/**
 *
 * Copyright (C) 2009-2020 Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.helper;

/**
 * Class to construct a Binary tree for McGregor search.
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class BinaryTree {

    /**
     * @param value the value to set
     */
    public void setValue(int value) {
        this.value = value;
    }

    /**
     * Creates a new instance of BinaryTree.
     *
     * @param value node value
     */
    public BinaryTree(int value) {
        this.value = value;
        this.equal = null;
        this.notEqual = null;
    }

    public BinaryTree() {
        this.value = -1;
        this.equal = null;
        this.notEqual = null;
    }
    /**
     * declaration
     */
    private BinaryTree equal;
    private BinaryTree notEqual;
    private int value;

    /**
     * Return value of the node
     *
     * @return get the value of the current node
     */
    public synchronized int getValue() {
        return this.value;
    }

    /**
     * Returns equal node
     *
     * @return the equal
     */
    public synchronized BinaryTree getEqual() {
        return equal;
    }

    /**
     * Set equal node
     *
     * @param tree the equal to set
     */
    public synchronized void setEqual(BinaryTree tree) {
        this.equal = tree;
    }

    /**
     * Returns not equal node
     *
     * @return the notEqual
     */
    public synchronized BinaryTree getNotEqual() {
        return notEqual;
    }

    /**
     * Set not equal node
     *
     * @param tree the tree to set
     */
    public synchronized void setNotEqual(BinaryTree tree) {
        this.notEqual = tree;
    }

    /* 
     * Remove the nodes recursively 
     *
     *
     * @param cur_struc
     * @return
     */
    public static int remove_tree_structure(BinaryTree cur_struc) {

        BinaryTree equal_struc = cur_struc.equal;
        BinaryTree not_equal_struc = cur_struc.notEqual;
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
