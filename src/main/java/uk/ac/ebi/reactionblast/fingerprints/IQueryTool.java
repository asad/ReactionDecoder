/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.fingerprints;

import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public interface IQueryTool {

    /**
     * Returns the number of times the pattern was found in the target molecule.
     * <p/>
     * This function should be called after
     * {@link #matches(org.openscience.cdk.interfaces.IAtomContainer)}. If not,
     * the results may be undefined.
     *
     * @return The number of times the pattern was found in the target molecule
     */
    int countMatches();

    /**
     * Get the atoms in the target molecule that match the query pattern.
     * <p/>
     * Since there may be multiple matches, the return value is a List of List
     * objects. Each List object contains the indices of the atoms in the target
     * molecule, that match the query pattern
     *
     * @return A List of List of atom indices in the target molecule
     */
    List<List<Integer>> getMatchingAtoms();

    /**
     * Returns the current SMARTS pattern being used.
     *
     * @return The SMARTS pattern
     */
    String getSmarts();

    /**
     * Get the atoms in the target molecule that match the query pattern.
     * <p/>
     * Since there may be multiple matches, the return value is a List of List
     * objects. Each List object contains the unique set of indices of the atoms
     * in the target molecule, that match the query pattern
     *
     * @return A List of List of atom indices in the target molecule
     */
    List<List<Integer>> getUniqueMatchingAtoms();

    /**
     * Perform a SMARTS match and check whether the query is present in the
     * target molecule.
     * <p/>
     * This function simply checks whether the query pattern matches the
     * specified molecule. However the function will also, internally, save the
     * mapping of query atoms to the target molecule
     * <p/>
     * <b>Note</b>: This method performs a simple caching scheme, by comparing
     * the current molecule to the previous molecule by reference. If you
     * repeatedly match different SMARTS on the same molecule, this method will
     * avoid initializing ( ring perception, aromaticity etc.) the molecule each
     * time. If however, you modify the molecule between such multiple matchings
     * you should use the other form of this method to force initialization.
     *
     * @param atomContainer The target moleculoe
     * @return true if the pattern is found in the target molecule, false
     * otherwise
     * @throws CDKException if there is an error in ring, aromaticity or
     * isomorphism perception
     * @see #getMatchingAtoms()
     * @see #countMatches()
     * @see #matches(org.openscience.cdk.interfaces.IAtomContainer, boolean)
     */
    boolean matches(IAtomContainer atomContainer) throws CDKException;

    /**
     * Perform a SMARTS match and check whether the query is present in the
     * target molecule.
     * <p/>
     * This function simply checks whether the query pattern matches the
     * specified molecule. However the function will also, internally, save the
     * mapping of query atoms to the target molecule
     *
     * @param atomContainer The target moleculoe
     * @param forceInitialization If true, then the molecule is initialized
     * (ring perception, aromaticity etc). If false, the molecule is only
     * initialized if it is different (in terms of object reference) than one
     * supplied in a previous call to this method.
     * @return true if the pattern is found in the target molecule, false
     * otherwise
     * @throws CDKException if there is an error in ring, aromaticity or
     * isomorphism perception
     * @see #getMatchingAtoms()
     * @see #countMatches()
     * @see #matches(org.openscience.cdk.interfaces.IAtomContainer)
     */
    boolean matches(IAtomContainer atomContainer, boolean forceInitialization) throws CDKException;

    /**
     * Set the maximum size of the query cache.
     *
     * @param maxEntries The maximum number of entries
     */
    void setQueryCacheSize(int maxEntries);

    /**
     * Set a new SMARTS pattern.
     *
     * @param smarts The new SMARTS pattern
     * @throws CDKException if there is an error in parsing the pattern
     */
    void setSmarts(String smarts) throws CDKException;

}
