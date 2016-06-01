/* Copyright (C) 2011  Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class SuperAtoms {
    private static final Logger LOG = getLogger(SuperAtoms.class.getName());

    private String sgroupType;
    private IPseudoAtom pseudoAtom;
    private int index;
    private IAtom prevAtom;
    private int sgroupIndex;
    private IBond crossingBond;

    /**
     * @return the sgroupType
     */
    public String getSgroupType() {
        return sgroupType;
    }

    /**
     * @param sgroupType the sgroupType to set
     */
    public void setSgroupType(String sgroupType) {
        this.sgroupType = sgroupType;
    }

    /**
     * @return the pseudoAtom
     */
    public IPseudoAtom getPseudoAtom() {
        return pseudoAtom;
    }

    /**
     * @param pseudoAtom the pseudoAtom to set
     */
    public void setPseudoAtom(IPseudoAtom pseudoAtom) {
        this.pseudoAtom = pseudoAtom;
    }

    /**
     * @return the index
     */
    public int getIndex() {
        return index;
    }

    /**
     * @param index the index to set
     */
    public void setIndex(int index) {
        this.index = index;
    }

    /**
     * @return the prevAtom
     */
    public IAtom getPrevAtom() {
        return prevAtom;
    }

    /**
     * @param prevAtom the prevAtom to set
     */
    public void setPrevAtom(IAtom prevAtom) {
        this.prevAtom = prevAtom;
    }

    /**
     * @return the sgroupIndex
     */
    public int getSgroupIndex() {
        return sgroupIndex;
    }

    /**
     * @param sgroupIndex the sgroupIndex to set
     */
    public void setSgroupIndex(int sgroupIndex) {
        this.sgroupIndex = sgroupIndex;
    }

    /**
     * @return the crossingBond
     */
    public IBond getCrossingBond() {
        return crossingBond;
    }

    /**
     * @param crossingBond the crossingBond to set
     */
    public void setCrossingBond(IBond crossingBond) {
        this.crossingBond = crossingBond;
    }
}
