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

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IStereoElement;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.stereo.compare.StereoCenterAnalyser.hasPotentialStereoCenter;
import static uk.ac.ebi.reactionblast.stereo.tools.Chirality2DTool.getChirality2D;
import static uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult.CHIRAL_R;
import static uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult.CHIRAL_S;
import static uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult.ERROR;
import static uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult.MISSING;
import static uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult.NONE;

/**
 * Analyse the stereo wedges around an atom, to determine if they are correct.
 *
 * @author maclean
 *
 */
public class WedgeStereoAnalyser {

    /**
     *
     * @param atom
     * @param atomContainer
     * @param lifter
     * @return
     */
    public static WedgeStereoAnalysisResult getResult(IAtom atom, IAtomContainer atomContainer, WedgeStereoLifter lifter) {
        boolean isPotentialStereoCenter
                = hasPotentialStereoCenter(atom, atomContainer);
        IStereoElement element = lifter.lift(atom, atomContainer);
        return getResult(atomContainer, isPotentialStereoCenter, element);
    }

    private static WedgeStereoAnalysisResult getResult(IAtomContainer atomContainer, boolean isPotentialStereoCenter, IStereoElement stereoElement) {
        if (isPotentialStereoCenter) {
            if (stereoElement == null) {
                return MISSING;
            } else {
                IStereoAndConformation chirality = getChirality2D(stereoElement, atomContainer);
                WedgeStereoAnalysisResult result = convertCipToResult(chirality);
                if (result == NONE) {
                    // should have R or S!
                    return ERROR;
                } else {
                    return result;
                }
            }
        } else if (stereoElement == null) {
            return NONE;
        } else {
            return ERROR;
        }
    }

    private static WedgeStereoAnalysisResult convertCipToResult(IStereoAndConformation cipChirality) {
        switch (cipChirality) {
            case NONE:
                return NONE;
            case R:
                return CHIRAL_R;
            case S:
                return CHIRAL_S;
            default:
                return NONE;
        }
    }

    private WedgeStereoAnalyser() {
    }
}
