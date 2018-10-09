/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.interfaces;

import java.io.IOException;
import java.util.List;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public interface IGraphTheoryMatrix {

    /**
     * clears all the containers in this object
     *
     * @throws java.io.IOException
     */
    void Clear() throws IOException;

    /**
     *
     * @return
     */
    int getDelta();

    /**
     *
     * @return reactant molecule and their index (key) in a Map
     */
    List<String> getEductCounter();

    /**
     * @return the matrixHolder
     */
    Holder getMatrixHolder();

    /**
     *
     * @return product molecule and their index (key) in a Map
     */
    List<String> getProductCounter();

    /**
     *
     * @return
     */
    MoleculeMoleculeMapping getReactionMolMapping();

    /**
     *
     * @param reactionMolMapping
     */
    void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping);
}
