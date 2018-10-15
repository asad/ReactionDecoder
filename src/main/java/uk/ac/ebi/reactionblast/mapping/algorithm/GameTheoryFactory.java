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
package uk.ac.ebi.reactionblast.mapping.algorithm;

import java.io.Serializable;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.mapping.interfaces.IGameTheory;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;

/**
 * This class initiates the algorithm
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class GameTheoryFactory implements Serializable {

    private static final long serialVersionUID = 01567272317571L;

    /**
     *
     * @param theory
     * @param reaction
     * @param removeHydrogen
     * @param educts
     * @param products
     * @param rpsh
     * @return
     * @throws Exception
     */
    public static synchronized IGameTheory make(IMappingAlgorithm theory, IReaction reaction, boolean removeHydrogen, Map<Integer, IAtomContainer> educts, Map<Integer, IAtomContainer> products, GameTheoryMatrix rpsh) throws Exception {
        switch (theory) {
            case MIXTURE:
                return new GameTheoryMixture(
                        reaction, removeHydrogen,
                        educts, products,
                        rpsh);
            case MIN:
                return new GameTheoryMin(
                        reaction, removeHydrogen,
                        educts, products,
                        rpsh);
            case MAX:
                return new GameTheoryMax(
                        reaction, removeHydrogen,
                        educts, products,
                        rpsh);
            case RINGS:
                return new GameTheoryRings(
                        reaction, removeHydrogen,
                        educts, products,
                        rpsh);
            default:
                return null;
        }
    }

    private GameTheoryFactory() {
    }
}
