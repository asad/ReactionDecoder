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
package uk.ac.ebi.reactionblast.mapping.helper;

import java.io.IOException;
import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public abstract class AbstractReactor extends MappingHandler {

    /**
     *
     */
    public AbstractReactor() {
    }

    /**
     *
     * @return bonds of reactants
     */
    public abstract List<IBond> getEductBonds();

    /**
     *
     * @param i Index Ith position
     * @return Stoichiometry weight of the product molecule at ith Position
     *
     */
    public abstract Double getExpandedProductStoichiometry(int i);

    /**
     *
     * @return products expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
    public abstract IAtomContainerSet getExpandedProducts() throws IOException;

    /**
     *
     * @param i Index Ith position
     * @return Stoichiometry weight of the reactant molecule at ith Position
     */
    public abstract Double getExpandedReactantStoichiometry(int i);

    /**
     *
     * @return reactants expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
    public abstract IAtomContainerSet getExpandedReactants() throws IOException;

    /**
     *
     * @return
     */
    public abstract int getMappingCount();

    /**
     * @return bonds of products
     */
    public abstract List<IBond> getProductBonds();

    /**
     *
     * @return true if its a balanced reaction else false
     * @throws IOException
     *
     */
    public abstract boolean getReactionBalanceFlag() throws IOException;

    /**
     *
     * @return @throws java.io.IOException
     */
    public abstract boolean getReactionBalanceFlagWithChargeBalance() throws IOException;

    /**
     *
     * @return true if its a balanced reaction else false Note: This does not
     * consider whether Hydrogens are balanced or not
     *
     *
     */
    public abstract boolean getReactionBalanceFlagWithoutHydrogen();

    /**
     *
     * @return IReaction object with unique atom labeling
     * @throws Exception
     */
    public abstract IReaction getReactionWithAtomAtomMapping() throws Exception;
}
