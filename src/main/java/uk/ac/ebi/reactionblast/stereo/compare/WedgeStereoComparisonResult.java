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
package uk.ac.ebi.reactionblast.stereo.compare;

import static java.lang.String.valueOf;


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

    private final IAtom firstAtom;
    private final IAtomContainer firstAtomContainer;
    private final WedgeStereoAnalysisResult resultForFirst;
    private final IAtom secondAtom;
    private final IAtomContainer secondAtomContainer;
    private final WedgeStereoAnalysisResult resultForSecond;

    /**
     *
     * @param firstAtom
     * @param firstAtomContainer
     * @param resultForFirst
     * @param secondAtom
     * @param secondAtomContainer
     * @param resultForSecond
     */
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

    /**
     *
     * @return
     */
    public int getIndexOfFirst() {
        return firstAtomContainer.indexOf(firstAtom);
    }

    /**
     *
     * @return
     */
    public WedgeStereoAnalysisResult getResultForFirst() {
        return resultForFirst;
    }

    /**
     *
     * @return
     */
    public int getIndexOfSecond() {
        return secondAtomContainer.indexOf(secondAtom);
    }

    /**
     *
     * @return
     */
    public WedgeStereoAnalysisResult getResultForSecond() {
        return resultForSecond;
    }

    /**
     *
     * @return
     */
    public IAtom getFirst() {
        return firstAtom;
    }

    /**
     *
     * @return
     */
    public IAtomContainer getFirstContainer() {
        return firstAtomContainer;
    }

    /**
     *
     * @return
     */
    public IAtom getSecond() {
        return secondAtom;
    }

    /**
     *
     * @return
     */
    public IAtomContainer getSecondContainer() {
        return secondAtomContainer;
    }

    @Override
    public String toString() {
        String firstID;
        if (firstAtom.getID() == null) {
            firstID = valueOf(firstAtomContainer.indexOf(firstAtom));
        } else {
            firstID = firstAtom.getID();
        }
        String secondID;
        if (secondAtom.getID() == null) {
            secondID = valueOf(secondAtomContainer.indexOf(secondAtom));
        } else {
            secondID = secondAtom.getID();
        }
        return firstID + "\t:\t" + resultForFirst + "\t\t"
                + secondID + "\t:\t" + resultForSecond;
    }
}
