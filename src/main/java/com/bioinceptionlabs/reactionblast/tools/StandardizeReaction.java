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
package com.bioinceptionlabs.reactionblast.tools;

import static java.lang.System.currentTimeMillis;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import com.bioinceptionlabs.reactionblast.mapping.container.CDKReactionBuilder;
import static com.bioinceptionlabs.reactionblast.mapping.helper.MappingHandler.cleanMapping;

/**
 * Standardizes reaction objects for atom-atom mapping.
 * Validates atom balance and prepares reaction containers.
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class StandardizeReaction {

    private static final ILoggingTool LOGGER = createLoggingTool(StandardizeReaction.class);

    /**
     * Standardize a reaction: clean mappings, validate balance, build containers.
     *
     * @param reaction the input reaction
     * @return New Standardized reaction Object
     * @throws Exception if standardization fails
     */
    public IReaction standardize(IReaction reaction) throws Exception {
        String reactionID = reaction.getID();
        cleanMapping(reaction);

        if (reactionID == null) {
            reactionID = Long.toString(currentTimeMillis());
            reaction.setID(reactionID);
        }

        // Validate atom balance (warn but don't fail — some reactions are intentionally unbalanced)
        checkAtomBalance(reaction);

        CDKReactionBuilder rBuilder = new CDKReactionBuilder();
        return rBuilder.standardize(reaction);
    }

    /**
     * Check if a reaction is atom-balanced. Logs a warning if not.
     * Does not throw — unbalanced reactions are handled gracefully.
     *
     * @param reaction the reaction to check
     */
    private void checkAtomBalance(IReaction reaction) {
        Map<String, Integer> reactantAtoms = countAtoms(reaction.getReactants());
        Map<String, Integer> productAtoms = countAtoms(reaction.getProducts());

        if (!reactantAtoms.equals(productAtoms)) {
            LOGGER.warn("Reaction " + reaction.getID() + " may be unbalanced: "
                    + "reactants=" + reactantAtoms + " products=" + productAtoms);
        }
    }

    private Map<String, Integer> countAtoms(org.openscience.cdk.interfaces.IAtomContainerSet molSet) {
        Map<String, Integer> counts = new LinkedHashMap<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                counts.merge(atom.getSymbol(), 1, Integer::sum);
            }
        }
        return counts;
    }
}
