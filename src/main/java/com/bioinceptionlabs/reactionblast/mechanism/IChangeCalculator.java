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
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface IChangeCalculator {

    public BEMatrix getEductBEMatrix();
    public BEMatrix getProductBEMatrix();
    public RMatrix getRMatrix();
    public void printBMatrix();
    public void printEMatrix();
    public void printRMatrix();
    public void writeBMatrix(File outputFile);
    public void writeEMatrix(File outputFile);
    public void writeRMatrix(File outputFile);
    public boolean hasRMatrix();
    public Map<IAtom, IAtom> getMappingMap();
    public List<MechanismHelpers.BondChange> getBondChangeList();
    public Collection<IAtom> getReactionCenterSet();
    public List<MechanismHelpers.AtomStereoChangeInformation> getStereoChangeList();
    public Iterable<MechanismHelpers.AtomStereoChangeInformation> getConformationChangeList();
    public MechanismHelpers.AtomAtomMappingContainer getMappingContainer();
}
