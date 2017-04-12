/*
 *
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *                          Gilleain Torrance <gilleain.torrance@gmail.com>
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
 * 
 * 
 ** Copyright (C) 2009-2015 Kyle Lutz <kyle.r.lutz@gmail.com>
 **
 ** This file is part of chemkit. For more information see
 ** <http://www.chemkit.org>.
 **
 ** chemkit is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU Lesser General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** chemkit is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU Lesser General Public License for more details.
 **
 ** You should have received a copy of the GNU Lesser General Public License
 ** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
 **
 ******************************************************************************/
package org.openscience.smsd.algorithm.vflib.substructure;

/**
 * Holds source and target atoms
 * 
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
/**
 * @param <T>
 * @param <S>
 */
public class Pair<T, S> {

    private T source;
    private S target;

    public Pair(T a, S b) {
        this.source = a;
        this.target = b;
    }

    @Override
    public synchronized String toString() {
        return "(" + getSourceAtom() + ", " + getTargetAtom() + ")";
    }

    /**
     * @return the source
     */
    public synchronized T getSourceAtom() {
        return source;
    }

    /**
     * @param source the source to set
     */
    public synchronized void setSourceAtom(T first) {
        this.source = first;
    }

    /**
     * @return the target
     */
    public synchronized S getTargetAtom() {
        return target;
    }

    /**
     * @param target the target to set
     */
    public synchronized void setTargetAtom(S second) {
        this.target = second;
    }
}
