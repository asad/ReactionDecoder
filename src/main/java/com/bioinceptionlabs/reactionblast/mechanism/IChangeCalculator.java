/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;

/**
 * Package-private interface for change calculators.
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
interface IChangeCalculator {

    BEMatrix getEductBEMatrix();
    BEMatrix getProductBEMatrix();
    RMatrix getRMatrix();
    void printBMatrix();
    void printEMatrix();
    void printRMatrix();
    void writeBMatrix(File outputFile);
    void writeEMatrix(File outputFile);
    void writeRMatrix(File outputFile);
    boolean hasRMatrix();
    Map<IAtom, IAtom> getMappingMap();
    List<MechanismHelpers.BondChange> getBondChangeList();
    Collection<IAtom> getReactionCenterSet();
    List<MechanismHelpers.AtomStereoChangeInformation> getStereoChangeList();
    Iterable<MechanismHelpers.AtomStereoChangeInformation> getConformationChangeList();
    MechanismHelpers.AtomAtomMappingContainer getMappingContainer();
}
