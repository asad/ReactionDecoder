/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IStereoElement;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import uk.ac.ebi.reactionblast.stereo.compare.StereoCenterAnalyser;
import uk.ac.ebi.reactionblast.stereo.tools.Chirality2DTool;

/**
 * Analyse the stereo wedges around an atom, to determine if they are correct.
 *
 * @author maclean
 *
 */
public class WedgeStereoAnalyser {

    private static final Logger LOG = Logger.getLogger(WedgeStereoAnalyser.class.getName());

    public static WedgeStereoAnalysisResult getResult(IAtom atom, IAtomContainer atomContainer, WedgeStereoLifter lifter) {
        boolean isPotentialStereoCenter
                = StereoCenterAnalyser.hasPotentialStereoCenter(atom, atomContainer);
        IStereoElement element = lifter.lift(atom, atomContainer);
        return getResult(atomContainer, isPotentialStereoCenter, element);
    }

    private static WedgeStereoAnalysisResult getResult(IAtomContainer atomContainer, boolean isPotentialStereoCenter, IStereoElement stereoElement) {
        if (isPotentialStereoCenter) {
            if (stereoElement == null) {
                return WedgeStereoAnalysisResult.MISSING;
            } else {
                IStereoAndConformation chirality = Chirality2DTool.getChirality2D(stereoElement, atomContainer);
                WedgeStereoAnalysisResult result = convertCipToResult(chirality);
                if (result == WedgeStereoAnalysisResult.NONE) {
                    // should have R or S!
                    return WedgeStereoAnalysisResult.ERROR;
                } else {
                    return result;
                }
            }
        } else {
            if (stereoElement == null) {
                return WedgeStereoAnalysisResult.NONE;
            } else {
                return WedgeStereoAnalysisResult.ERROR;
            }
        }
    }

    private static WedgeStereoAnalysisResult convertCipToResult(IStereoAndConformation cipChirality) {
        switch (cipChirality) {
            case NONE:
                return WedgeStereoAnalysisResult.NONE;
            case R:
                return WedgeStereoAnalysisResult.CHIRAL_R;
            case S:
                return WedgeStereoAnalysisResult.CHIRAL_S;
            default:
                return WedgeStereoAnalysisResult.NONE;
        }
    }
}
