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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.interfaces;

import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;

/**
 * Interface for all MCS/Substructure algorithms.
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public interface IAtomMapping {

    /**
     * initialize query and target molecules.
     *
     * @param stereoFilter set true to rank the solutions as per stereo matches
     * @param fragmentFilter set true to return matches with minimum fragments
     * @param energyFilter set true to return matches with minimum bond changes
     * based on the bond breaking energy
     */
    public abstract void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter);

    /** 
     * Returns summation energy score of the disorder if the MCS is removed
     * from the target and query graph. Amongst the solutions, a solution
     * with lowest energy score is preferred.
     * 
     * @param Key Index of the mapping solution
     * @return Total bond breaking energy required to remove the mapped part
     */
    public abstract Double getEnergyScore(int Key);

    /** 
     * Returns number of fragment generated in the solution space,
     * if the MCS is removed from the target and query graph.
     * Amongst the solutions, a solution with lowest fragment size
     * is preferred.
     *
     * @param Key Index of the mapping solution
     * @return Fragment count(s) generated after removing the mapped parts
     */
    public abstract Integer getFragmentSize(int Key);

    /** 
     * Returns a number which denotes the quality of the mcs.
     * A solution with highest stereo score is preferred over other
     * scores.
     * @param Key Index of the mapping solution
     * @return true if no stereo mismatch occurs
     * else false if stereo mismatch occurs
     */
    public abstract Integer getStereoScore(int Key);

    /**
     * Returns all plausible mappings between query and target molecules
     * Each map in the list has atom-atom equivalence of the mappings
     * between query and target molecule i.e. map.getKey() for the query
     * and map.getValue() for the target molecule.
     * @return All possible MCS atom Mappings
     */
    public abstract List<AtomAtomMapping> getAllAtomMapping();

    /**
     * Returns one of the best matches with atoms mapped.
     * @return Best Atom Mapping
     */
    public abstract AtomAtomMapping getFirstAtomMapping();

    /** 
     * Returns Tanimoto similarity between query and target molecules
     * (Score is between 0-min and 1-max).
     *
     * @return Tanimoto Similarity between 0 and 1
     */
    public abstract double getTanimotoSimilarity();

    /** 
     * Returns Euclidean Distance between query and target molecule.
     * @return Euclidean Distance (lower the score, better the match)
     *
     */
    public abstract double getEuclideanDistance();

    /**
     *
     * Returns true if mols have different stereo
     * chemistry else false if no stereo mismatch.
     * 
     * @return true if mols have different stereo
     * chemistry else false if no stereo mismatch.
     * true if stereo mismatch occurs
     * else true if stereo mismatch occurs.
     */
    public abstract boolean isStereoMisMatch();

    /**
     * Returns total mapping count between query and target molecules
     * Each map in the list has atom-atom equivalence of the mappings
     * between query and target molecule i.e. map.getKey() for the query
     * and map.getValue() for the target molecule.
     * @return All possible MCS atom Mappings
     */
    public abstract int getMappingCount();
    
    /**
     * Convert to String
     * @return
     */
    @Override
    public abstract String toString();
}
