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
package uk.ac.ebi.reactionblast.stereo.compare;

import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.stereo.wedge.WedgeStereoAnalysisResult;

/**
 * Tool for comparing chiralities.
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 *
 */
public class WedgeStereoComparisonResult {
    private static final Logger LOG = Logger.getLogger(WedgeStereoComparisonResult.class.getName());

    private final IAtom firstAtom;
    private final IAtomContainer firstAtomContainer;
    private final WedgeStereoAnalysisResult resultForFirst;
    private final IAtom secondAtom;
    private final IAtomContainer secondAtomContainer;
    private final WedgeStereoAnalysisResult resultForSecond;

    public WedgeStereoComparisonResult(
            IAtom firstAtom, IAtomContainer firstAtomContainer, WedgeStereoAnalysisResult resultForFirst,
            IAtom secondAtom, IAtomContainer secondAtomContainer, WedgeStereoAnalysisResult resultForSecond) {
        this.firstAtom = firstAtom;
        this.firstAtomContainer = firstAtomContainer;
        this.resultForFirst = resultForFirst;
        this.secondAtom = secondAtom;
        this.secondAtomContainer = secondAtomContainer;
        this.resultForSecond = resultForSecond;
    }

    public int getIndexOfFirst() {
        return firstAtomContainer.getAtomNumber(firstAtom);
    }

    public WedgeStereoAnalysisResult getResultForFirst() {
        return resultForFirst;
    }

    public int getIndexOfSecond() {
        return secondAtomContainer.getAtomNumber(secondAtom);
    }

    public WedgeStereoAnalysisResult getResultForSecond() {
        return resultForSecond;
    }

    public IAtom getFirst() {
        return firstAtom;
    }

    public IAtomContainer getFirstContainer() {
        return firstAtomContainer;
    }

    public IAtom getSecond() {
        return secondAtom;
    }

    public IAtomContainer getSecondContainer() {
        return secondAtomContainer;
    }

    public String toString() {
        String firstID;
        if (firstAtom.getID() == null) {
            firstID = String.valueOf(firstAtomContainer.getAtomNumber(firstAtom));
        } else {
            firstID = firstAtom.getID();
        }
        String secondID;
        if (secondAtom.getID() == null) {
            secondID = String.valueOf(secondAtomContainer.getAtomNumber(secondAtom));
        } else {
            secondID = secondAtom.getID();
        }
        return firstID + "\t:\t" + resultForFirst + "\t\t"
                + secondID + "\t:\t" + resultForSecond;
    }
}
