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
import java.util.Collection;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;
import uk.ac.ebi.reactionblast.mapping.graph.MCSSolution;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public interface IGameTheory {

    /**
     * @return the delta
     */
    public int getDelta();

    /**
     * @return the reactionMolMapping
     */
    public abstract MoleculeMoleculeMapping getReactionMolMapping();

    /**
     *
     * @return @throws IOException
     */
    public String getSuffix() throws IOException;

    /**
     * @param reactionMolMapping the reactionMolMapping to set
     */
    public abstract void setReactionMolMapping(MoleculeMoleculeMapping reactionMolMapping);

    /**
     *
     * @param mh
     * @param removeHydrogen
     * @throws Exception
     */
    public void UpdateMatrix(Holder mh, boolean removeHydrogen) throws Exception;

    /**
     *
     * @param mcsSolutions
     * @param mh
     * @param removeHydrogen
     * @throws Exception
     */
    public void UpdateMatrix(Collection<MCSSolution> mcsSolutions, Holder mh, boolean removeHydrogen) throws Exception;
}
