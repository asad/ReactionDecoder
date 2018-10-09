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
package uk.ac.ebi.reactionblast.mechanism.interfaces;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import uk.ac.ebi.reactionblast.mechanism.BEMatrix;
import uk.ac.ebi.reactionblast.mechanism.RMatrix;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomStereoChangeInformation;
import uk.ac.ebi.reactionblast.mechanism.helper.BondChange;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public interface IChangeCalculator {

    /**
     *
     * @return
     */
    public BEMatrix getEductBEMatrix();

    /**
     *
     * @return
     */
    public BEMatrix getProductBEMatrix();

    /**
     *
     * @return
     */
    public RMatrix getRMatrix();

    /**
     *
     */
    public void printBMatrix();

    /**
     *
     */
    public void printEMatrix();

    /**
     *
     */
    public void printRMatrix();

    /**
     *
     * @param outputFile
     */
    public void writeBMatrix(File outputFile);

    /**
     *
     * @param outputFile
     */
    public void writeEMatrix(File outputFile);

    /**
     *
     * @param outputFile
     */
    public void writeRMatrix(File outputFile);

    /**
     *
     * @return
     */
    public boolean hasRMatrix();

    /**
     *
     * @return
     */
    public Map<IAtom, IAtom> getMappingMap();

    /**
     *
     * @return
     */
    public List<BondChange> getBondChangeList();

    /**
     *
     * @return
     */
    public Collection<IAtom> getReactionCenterSet();

    /**
     *
     * @return
     */
    public List<AtomStereoChangeInformation> getStereoChangeList();

    /**
     *
     * @return
     */
    public Iterable<AtomStereoChangeInformation> getConformationChangeList();

    /**
     *
     * @return
     */
    public AtomAtomMappingContainer getMappingContainer();
}
