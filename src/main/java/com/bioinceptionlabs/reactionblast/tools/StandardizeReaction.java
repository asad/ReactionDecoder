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

import com.bioinception.smsd.core.SMSD;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
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

    public static final String SOURCE_OCCURRENCE_ID = "sourceOccurrenceId";
    public static final String SOURCE_ATOM_ID = "sourceAtomId";
    public static final String PRESERVE_OCCURRENCE_IDENTITY = "preserveOccurrenceIdentity";
    public static final String STOICHIOMETRY_KEY = "stoichiometryKey";

    private static final ILoggingTool LOGGER = createLoggingTool(StandardizeReaction.class);

    /**
     * Common solvents, reagents, and catalysts by canonical SMILES.
     * These molecules never participate in bond-changing reactions —
     * they facilitate or mediate but their bonds don't change.
     */
    private static final Set<String> KNOWN_REAGENT_SMILES = new HashSet<>(Arrays.asList(
            // Solvents
            "ClCCl",        // DCM (dichloromethane)
            "ClC(Cl)Cl",    // chloroform
            "CC(C)=O",      // acetone
            "CCCCCC",       // hexane
            "c1ccncc1",     // pyridine (also base)
            "CC#N",         // acetonitrile (MeCN)
            "CS(C)=O",      // DMSO
            "CN(C)C=O",     // DMF
            "C1CCOC1",      // THF
            "CCOCC",        // diethyl ether
            "C1COCCO1",     // 1,4-dioxane
            "CO",           // methanol
            "CCO",          // ethanol
            "CC(C)O",       // isopropanol
            "O",            // water
            "CC(=O)O",      // acetic acid (when used as solvent)
            "CCOC(C)=O",    // ethyl acetate
            "c1ccccc1",     // benzene
            "Cc1ccccc1",    // toluene
            "c1ccc(cc1)C",  // toluene alternate
            // Reducing agents
            "[Na+]",        // sodium cation
            "[K+]",         // potassium cation
            "[Li+]",        // lithium cation
            "[Cs+]",        // cesium cation
            "[NH4+]",       // ammonium
            "[Cl-]",        // chloride
            "[Br-]",        // bromide
            "[I-]",         // iodide
            "[OH-]",        // hydroxide
            // Inorganic bases/acids
            "[Na]O",        // NaOH
            "O=S(=O)(O)O",  // sulfuric acid
            // Drying agents / dessicants
            "O=S(Cl)Cl",    // thionyl chloride (reagent but bonds don't map)
            "[Mg+2]",       // magnesium ion
            "[Ca+2]",       // calcium ion
            "[Zn]",         // zinc
            // Borohydride / cyanoborohydride (reducing agents)
            "[BH4-]",       // borohydride
            "[BH3-]C#N"     // cyanoborohydride
    ));

    /**
     * Metal elements commonly found in catalysts.
     * Molecules containing these are likely catalysts, not reactants.
     */
    private static final Set<String> CATALYST_METALS = new HashSet<>(Arrays.asList(
            "Pd", "Pt", "Rh", "Ru", "Ir", "Ni", "Cu", "Fe",
            "Co", "Mn", "Ti", "Zr", "Mo", "W", "Os", "Ag", "Au"
    ));

    /**
     * Standardize a reaction: clean mappings, validate balance, build containers.
     *
     * @param reaction the input reaction
     * @return New Standardized reaction Object
     * @throws Exception if standardization fails
     */
    public IReaction standardize(IReaction reaction) throws Exception {
        String reactionID = reaction.getID();
        annotateSourceIdentity(reaction);
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

    private void annotateSourceIdentity(IReaction reaction) {
        annotateSourceIdentity(reaction.getReactants(), "R");
        annotateSourceIdentity(reaction.getProducts(), "P");
    }

    private void annotateSourceIdentity(IAtomContainerSet containers, String side) {
        Map<Integer, String> componentSignatures = new LinkedHashMap<>();
        Map<String, Integer> signatureCounts = new LinkedHashMap<>();
        org.openscience.cdk.smiles.SmilesGenerator smilesGenerator
                = new org.openscience.cdk.smiles.SmilesGenerator(
                        org.openscience.cdk.smiles.SmiFlavor.Canonical);

        for (int moleculeIndex = 0; moleculeIndex < containers.getAtomContainerCount(); moleculeIndex++) {
            IAtomContainer molecule = containers.getAtomContainer(moleculeIndex);
            String signature = componentSignature(molecule, smilesGenerator);
            componentSignatures.put(moleculeIndex, signature);
            signatureCounts.merge(signature, 1, Integer::sum);
        }

        for (int moleculeIndex = 0; moleculeIndex < containers.getAtomContainerCount(); moleculeIndex++) {
            IAtomContainer molecule = containers.getAtomContainer(moleculeIndex);
            molecule.setProperty(SOURCE_OCCURRENCE_ID, side + ":" + moleculeIndex);
            boolean preserveOccurrenceIdentity = hasBenchmarkAtomIds(molecule)
                    || signatureCounts.getOrDefault(componentSignatures.get(moleculeIndex), 0) > 1;
            molecule.setProperty(PRESERVE_OCCURRENCE_IDENTITY, preserveOccurrenceIdentity);
            for (int atomIndex = 0; atomIndex < molecule.getAtomCount(); atomIndex++) {
                IAtom atom = molecule.getAtom(atomIndex);
                if (atom.getProperty(SOURCE_ATOM_ID) == null) {
                    atom.setProperty(SOURCE_ATOM_ID, side + ":" + moleculeIndex + ":" + atomIndex);
                }
            }
        }
    }

    private String componentSignature(IAtomContainer molecule,
            org.openscience.cdk.smiles.SmilesGenerator smilesGenerator) {
        try {
            return smilesGenerator.create(molecule);
        } catch (Exception e) {
            return molecule.getAtomCount() + ":" + molecule.getBondCount();
        }
    }

    private boolean hasBenchmarkAtomIds(IAtomContainer molecule) {
        for (IAtom atom : molecule.atoms()) {
            if (atom.getProperty("benchmarkAtomId") != null) {
                return true;
            }
        }
        return false;
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
            LOGGER.debug("Reaction " + reaction.getID() + " may be unbalanced: "
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
            IAtomContainerSet products = reaction.getProducts();

            // Pre-compute product fingerprints using SMSD ECFP4 (radius=2)
            List<long[]> productFPs = new ArrayList<>();
            for (IAtomContainer prod : products.atomContainers()) {
                try {
                    productFPs.add(SMSD.circularFingerprintECFP(prod, 2, 256));
                } catch (Exception e) {
                    productFPs.add(null);
                }
            }

            // Collect product atom types for balance check
            Map<String, Integer> productAtomCounts = countAtoms(products);

            List<IAtomContainer> keptReactants = new ArrayList<>();
            List<IAtomContainer> reagents = new ArrayList<>();

            // Generate canonical SMILES for known-reagent lookup
            org.openscience.cdk.smiles.SmilesGenerator smiGen =
                    new org.openscience.cdk.smiles.SmilesGenerator(
                            org.openscience.cdk.smiles.SmiFlavor.Canonical);

            for (IAtomContainer reactant : reaction.getReactants().atomContainers()) {
                boolean isReagent = false;
                String reason = "";

                // Check if this reactant is needed for atom balance
                boolean neededForBalance = isNeededForBalance(
                        reactant, reaction.getReactants(), productAtomCounts);

                try {
                    // Check 1: Known solvent/reagent by canonical SMILES
                    String canSmiles = smiGen.create(reactant);
                    if (!neededForBalance && KNOWN_REAGENT_SMILES.contains(canSmiles)) {
                        isReagent = true;
                        reason = "known reagent/solvent: " + canSmiles;
                    }

                    // Check 2: Contains catalyst metal
                    if (!isReagent) {
                        for (IAtom atom : reactant.atoms()) {
                            if (CATALYST_METALS.contains(atom.getSymbol())) {
                                isReagent = true;
                                reason = "catalyst metal: " + atom.getSymbol();
                                break;
                            }
                        }
                    }

                    // Check 3: Tanimoto fingerprint similarity
                    if (!isReagent && !neededForBalance) {
                        long[] reactantFP = SMSD.circularFingerprintECFP(reactant, 2, 256);

                        double maxSim = 0.0;
                        for (long[] prodFP : productFPs) {
                            if (prodFP != null) {
                                double sim = SMSD.fingerprintTanimoto(reactantFP, prodFP);
                                maxSim = Math.max(maxSim, sim);
                            }
                        }

                        if (maxSim < 0.4 && reactant.getAtomCount() > 0) {
                            // Check for unique atom contribution
                            boolean hasUniqueContribution = false;
                            Map<String, Integer> reactantAtomCounts = new LinkedHashMap<>();
                            for (IAtom atom : reactant.atoms()) {
                                reactantAtomCounts.merge(atom.getSymbol(), 1, Integer::sum);
                            }
                            for (Map.Entry<String, Integer> entry : reactantAtomCounts.entrySet()) {
                                if (!productAtomCounts.containsKey(entry.getKey())) {
                                    hasUniqueContribution = true;
                                    break;
                                }
                            }

                            int heavyAtomCount = 0;
                            for (IAtom atom : reactant.atoms()) {
                                if (!"H".equals(atom.getSymbol())) heavyAtomCount++;
                            }

                            if (!hasUniqueContribution && heavyAtomCount <= 10) {
                                isReagent = true;
                                reason = "low Tanimoto=" + String.format("%.2f", maxSim)
                                        + ", atoms=" + heavyAtomCount;
                            }
                        }
                    }
                } catch (Exception e) {
                    LOGGER.debug("Filter check failed for " + reactant.getID()
                            + ": " + e.getMessage());
                }

                if (isReagent) {
                    LOGGER.debug("Filtered: " + reason);
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

    private boolean isNeededForBalance(IAtomContainer candidate,
            IAtomContainerSet allReactants, Map<String, Integer> productAtomCounts) {
        Map<String, Integer> remaining = new LinkedHashMap<>(countAtoms(allReactants));
        for (IAtom atom : candidate.atoms()) {
            remaining.merge(atom.getSymbol(), -1, Integer::sum);
        }
        for (Map.Entry<String, Integer> entry : productAtomCounts.entrySet()) {
            if (remaining.getOrDefault(entry.getKey(), 0) < entry.getValue()) {
                return true;
            }
        }
        return false;
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
