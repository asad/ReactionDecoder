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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;

import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.CDKReactionBuilder;
import static com.bioinceptionlabs.reactionblast.mapping.Reactor.MappingHandler.cleanMapping;

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

        // Filter reagents/solvents before mapping (improves accuracy for multi-component reactions)
        reaction = filterReagents(reaction);

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

    /**
     * Filter out reagents and solvents from reactants that don't participate
     * in the actual bond-changing reaction. Uses Tanimoto fingerprint similarity
     * to identify reactant molecules that have no corresponding product.
     *
     * A reactant is classified as a reagent/solvent if:
     * 1. Its max Tanimoto similarity to any product is < 0.3 (no product resembles it)
     * 2. AND it has no atoms that appear in the product atom balance
     *
     * This is conservative — it only removes molecules that clearly don't
     * participate. If in doubt, the molecule is kept as a reactant.
     *
     * @param reaction the reaction to filter
     * @return filtered reaction with reagents moved to agents
     */
    public IReaction filterReagents(IReaction reaction) {
        if (reaction.getReactantCount() <= 1) {
            return reaction; // nothing to filter
        }

        try {
            CircularFingerprinter fp = new CircularFingerprinter();
            IAtomContainerSet products = reaction.getProducts();

            // Pre-compute product fingerprints
            List<org.openscience.cdk.fingerprint.IBitFingerprint> productFPs = new ArrayList<>();
            for (IAtomContainer prod : products.atomContainers()) {
                try {
                    productFPs.add(fp.getBitFingerprint(prod));
                } catch (Exception e) {
                    productFPs.add(null);
                }
            }

            // Collect product atom types for balance check
            Map<String, Integer> productAtomCounts = countAtoms(products);

            List<IAtomContainer> keptReactants = new ArrayList<>();
            List<IAtomContainer> reagents = new ArrayList<>();

            for (IAtomContainer reactant : reaction.getReactants().atomContainers()) {
                boolean isReagent = false;
                try {
                    org.openscience.cdk.fingerprint.IBitFingerprint reactantFP = fp.getBitFingerprint(reactant);

                    // Find max similarity to any product
                    double maxSim = 0.0;
                    for (org.openscience.cdk.fingerprint.IBitFingerprint prodFP : productFPs) {
                        if (prodFP != null) {
                            double sim = Tanimoto.calculate(reactantFP, prodFP);
                            maxSim = Math.max(maxSim, sim);
                        }
                    }

                    // If no product resembles this reactant, it's likely a reagent
                    if (maxSim < 0.3 && reactant.getAtomCount() > 0) {
                        // Double-check: does any atom type in this molecule appear
                        // exclusively in the reactant (i.e., removed in products)?
                        // If so, it might still be a real reactant (leaving group)
                        boolean hasUniqueContribution = false;
                        Map<String, Integer> reactantAtomCounts = new LinkedHashMap<>();
                        for (IAtom atom : reactant.atoms()) {
                            reactantAtomCounts.merge(atom.getSymbol(), 1, Integer::sum);
                        }
                        // Check if this reactant contributes atoms not in products
                        for (Map.Entry<String, Integer> entry : reactantAtomCounts.entrySet()) {
                            if (!productAtomCounts.containsKey(entry.getKey())) {
                                hasUniqueContribution = true;
                                break;
                            }
                        }

                        // Only filter if: low similarity AND no unique atom contribution
                        // AND molecule is small (≤ 10 heavy atoms) — large molecules
                        // are more likely to be real reactants
                        int heavyAtomCount = 0;
                        for (IAtom atom : reactant.atoms()) {
                            if (!"H".equals(atom.getSymbol())) heavyAtomCount++;
                        }
                        if (!hasUniqueContribution && heavyAtomCount <= 10) {
                            isReagent = true;
                            LOGGER.debug("Filtered reagent/solvent: " + reactant.getID()
                                    + " (Tanimoto=" + String.format("%.2f", maxSim)
                                    + ", atoms=" + heavyAtomCount + ")");
                        }
                    }
                } catch (Exception e) {
                    // If fingerprinting fails, keep the molecule as reactant
                    LOGGER.debug("Fingerprint failed for " + reactant.getID() + ": " + e.getMessage());
                }

                if (isReagent) {
                    reagents.add(reactant);
                } else {
                    keptReactants.add(reactant);
                }
            }

            // Only filter if we'd keep at least 1 reactant
            if (keptReactants.isEmpty() || reagents.isEmpty()) {
                return reaction; // nothing filtered or would remove all
            }

            // Build filtered reaction
            IReaction filtered = reaction.getBuilder().newInstance(IReaction.class);
            filtered.setID(reaction.getID());
            filtered.setDirection(reaction.getDirection());
            for (IAtomContainer r : keptReactants) {
                filtered.addReactant(r);
            }
            for (IAtomContainer p : products.atomContainers()) {
                filtered.addProduct(p);
            }
            for (IAtomContainer agent : reagents) {
                filtered.addAgent(agent);
            }
            // Copy existing agents
            if (reaction.getAgents() != null) {
                for (IAtomContainer agent : reaction.getAgents().atomContainers()) {
                    filtered.addAgent(agent);
                }
            }

            LOGGER.debug("Filtered " + reagents.size() + " reagent(s) from "
                    + reaction.getReactantCount() + " reactants → "
                    + keptReactants.size() + " reactants");
            return filtered;

        } catch (Exception e) {
            LOGGER.debug("Reagent filtering failed: " + e.getMessage());
            return reaction; // return unfiltered on error
        }
    }

    private Map<String, Integer> countAtoms(IAtomContainerSet molSet) {
        Map<String, Integer> counts = new LinkedHashMap<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                counts.merge(atom.getSymbol(), 1, Integer::sum);
            }
        }
        return counts;
    }
}
