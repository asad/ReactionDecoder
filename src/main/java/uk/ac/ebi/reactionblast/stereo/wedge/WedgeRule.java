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
package uk.ac.ebi.reactionblast.stereo.wedge;

import java.util.SortedMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;

/**
 *
 * @author asad
 */
public abstract class WedgeRule {

    private int matchPoint;

    /**
     * A match between a particular bond stereo list and the pattern in a rule
     * may be cyclicly permuted - this method returns the permutation.
     *
     * @return a permutation as an int array
     */
    public int[] getMatchPermutation() {
        int length = getPattern().length;
        int[] permutation = new int[length];
        if (matchPoint != -1) {
            int pi = matchPoint;
            for (int i = 0; i < length; i++) {
                permutation[i] = pi;
                pi++;
                if (pi >= length) {
                    pi = 0;
                }
            }
        } else {
            // Didn't match - should have called this method!
        }
        return permutation;
    }

    /**
     * Match a list of {@link IBond.Stereo}s to the pattern list of stereos,
     * allowing for circular permutation of the list. As a side-effect, it
     * stores the point in the pattern where the match was made.
     *
     * @param stereoList the array of stereo constants to match
     * @return true if this Rule's pattern matches
     */
    public boolean matches(IBond.Stereo[] stereoList) {
        IBond.Stereo[] pattern = getPattern();
        if (stereoList.length != pattern.length) {
            return false;
        }
        int length = pattern.length;

        int patternIndex = 0;
        int matchIndex = 0;

        // reset the match point
        matchPoint = -1;
        while (patternIndex < length && matchIndex < 2 * length) {
            IBond.Stereo patternStereo = pattern[patternIndex];
            IBond.Stereo matchStereo;

            // XXX could fail faster here : pI - l > l?
            if (matchIndex < length) {
                matchStereo = stereoList[matchIndex];
            } else {
                matchStereo = stereoList[matchIndex - length];
            }

            if (patternStereo == matchStereo) {
                patternIndex++;
                matchIndex++;
            } else if (patternIndex == 0) {
                // look for a new start
                matchIndex++;
            } else {
                // start again
                patternIndex = 0;
            }
        }

        // less than a whole pattern was matched
        if (patternIndex < pattern.length) {
            return false;
        }

        // store the point where the match started
        matchPoint = matchIndex - length;

        return true;
    }

    /**
     *
     * @return
     */
    public abstract IBond.Stereo[] getPattern();

    /**
     *
     * @param centralAtom
     * @param atomContainer
     * @param angleMap
     * @return
     */
    public abstract IStereoElement execute(
            IAtom centralAtom,
            IAtomContainer atomContainer,
            SortedMap<Double, IBond> angleMap);
}
