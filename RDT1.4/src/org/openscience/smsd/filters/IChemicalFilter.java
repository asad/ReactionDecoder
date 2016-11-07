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
package org.openscience.smsd.filters;

import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.smsd.AtomAtomMapping;

/**
 * A filter on SMSD results.
 *
 * @param <T>
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * @author maclean
 * 
 */
public interface IChemicalFilter<T> {

    /**
     * Calculates a score for each MCS, and sorts the results on that score,
     * returning the best.
     *
     * @param allAtomMCS
     * @param selectionMap
     * @return
     * @throws CDKException
     */
    public T sortResults(
            Map<Integer, AtomAtomMapping> allAtomMCS,
            Map<Integer, T> selectionMap) throws CDKException;

    public List<T> getScores();

    public void clearScores();

    public void addScore(int counter, T value);

    public void fillMap(Map<Integer, T> map);
}
