/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping.algorithm;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import static java.lang.System.getProperty;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import static org.openscience.cdk.CDKConstants.RING_CONNECTIONS;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import static org.openscience.cdk.graph.Cycles.all;
import static org.openscience.cdk.graph.Cycles.or;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.MoleculeInitializer;
import org.openscience.smsd.ExtAtomContainerManipulator;
import static org.openscience.smsd.ExtAtomContainerManipulator.Utility.findSubgraph;
import static org.openscience.smsd.ExtAtomContainerManipulator.Utility.isMatch;
import static com.bioinceptionlabs.reactionblast.mapping.algorithm.GameTheoryEngine.GameTheoryFactory.make;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.MoleculeMoleculeMapping;
import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MAX;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MIN;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.MIXTURE;
import static com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm.RINGS;

/**
 * Combines CalculationProcess and CaseHandler (merged from CaseHandler.java).
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CalculationProcess implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static ILoggingTool LOGGER
            = createLoggingTool(CalculationProcess.class);
    private static final long serialVersionUID = 0x4a0bba049L;
    private final boolean removeHydrogen;
    private int delta = 0;
    private MoleculeMoleculeMapping reactionBlastMolMapping;
    private final IMappingAlgorithm algorithm;

    // --- Fields from CaseHandler ---
    private final Map<IRingSet, IAtomContainer> ringContainerCountR;
    private final Map<IRingSet, IAtomContainer> ringContainerCountP;
    private final IRingSet sssrEduct;
    private final IRingSet sssrProduct;
    protected final IReaction reaction;

    /**
     *
     * @param removeHydrogen
     * @param reaction
     * @param algorithm
     * @throws org.openscience.cdk.exception.Intractable
     */
    public CalculationProcess(
            boolean removeHydrogen,
            IReaction reaction,
            IMappingAlgorithm algorithm) throws Intractable {

        // --- CaseHandler initialization (merged) ---
        LOGGER.debug("=====CaseHandler====");
        this.reaction = reaction;
        ringContainerCountR = getRingContainerCount(reaction.getReactants());
        ringContainerCountP = getRingContainerCount(reaction.getProducts());

        if (ringContainerCountR.size() == 1 && ringContainerCountP.size() == 1) {
            IAtomContainer educt = ringContainerCountR.values().iterator().next();
            IAtomContainer product = ringContainerCountP.values().iterator().next();
            try {
                initializeMolecule(educt);
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
            try {
                initializeMolecule(product);
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }

            CycleFinder cf = Cycles.essential();
            Cycles cycles = cf.find(educt);
            sssrEduct = cycles.toRingSet();
            cycles = cf.find(product);
            sssrProduct = cycles.toRingSet();
            chipBondInTheRing(educt, product);
        } else {
            sssrEduct = null;
            sssrProduct = null;
            IAtomContainerSet reactants = reaction.getReactants();
            IAtomContainerSet products = reaction.getProducts();

            if (reactants.getAtomContainerCount() == 1 && products.getAtomContainerCount() == 1) {
                try {
                    IAtomContainer educt = reactants.atomContainers().iterator().next();
                    IAtomContainer product = products.atomContainers().iterator().next();
                    chipPhophateInSingleReactantProductNotInRing(educt, product);
                } catch (CDKException | IOException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
        for (IAtomContainer ac1 : reaction.getReactants().atomContainers()) {
            for (IAtomContainer ac2 : reaction.getProducts().atomContainers()) {
                try {
                    deleteBonds(ac1, ac2);
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
        LOGGER.debug("=====DONE CaseHandler====");

        // --- CalculationProcess initialization ---
        LOGGER.debug("=====CalculationProcess====");
        this.removeHydrogen = removeHydrogen;
        LOGGER.debug(NEW_LINE + "|++++++++++++++++++++++++++++|");
        LOGGER.debug("Performing Atom-Atom Mapping ....... " + reaction.getID() + " .......");
        LOGGER.debug(NEW_LINE + "|++++++++++++++++++++++++++++|");
        this.algorithm = algorithm;
        run();
        LOGGER.debug("=====Done CalculationProcess====");
    }

    private void run() {
        switch (algorithm) {
            case MIN:
                LOGGER.debug("Processing Reaction for Local Minimum: ");
                delta = (int) calRelation(reaction, MIN);
                break;
            case MAX:
                LOGGER.debug("Processing Reaction for Global Minimum: ");
                delta = (int) calRelation(reaction, MAX);
                break;
            case MIXTURE:
                LOGGER.debug("Processing Reaction for Max-Mixture Model: ");
                delta = (int) calRelation(reaction, MIXTURE);
                break;
            case RINGS:
                LOGGER.debug("Processing Reaction for Ring Model: ");
                delta = (int) calRelation(reaction, RINGS);
                break;
        }
    }

    /**
     *
     * @return
     */
    public IReaction getMappedReaction() {
        return reaction;
    }

    private double calRelation(IReaction reaction, IMappingAlgorithm theory) {
        try {
            Map<Integer, IAtomContainer> educts = new TreeMap<>();
            for (int i = 0; i < reaction.getReactantCount(); i++) {
                educts.put(i, reaction.getReactants().getAtomContainer(i));
            }

            Map<Integer, IAtomContainer> products = new TreeMap<>();
            for (int i = 0; i < reaction.getProductCount(); i++) {
                products.put(i, reaction.getProducts().getAtomContainer(i));
            }

            GameTheoryEngine.GameTheoryMatrix EDSH
                    = new GameTheoryEngine.GameTheoryMatrix(theory, reaction, removeHydrogen);

            LOGGER.debug("=====AGORITHM====" + theory);
            IGameTheory gameTheory = make(theory,
                    reaction,
                    removeHydrogen,
                    educts,
                    products,
                    EDSH);

            LOGGER.debug("=====DONE AGORITHM====" + theory);
            this.reactionBlastMolMapping = gameTheory.getReactionMolMapping();
            EDSH.Clear();

            return gameTheory.getDelta();
        } catch (Exception e) {
            LOGGER.error(e);
            return -1;
        }
    }

    /**
     * @return the delta
     */
    public int getDelta() {
        return delta;
    }

    /**
     * @return the reactionBlastMolMapping
     */
    public MoleculeMoleculeMapping getReactionBlastMolMapping() {
        return reactionBlastMolMapping;
    }

    // ---- Methods merged from CaseHandler.java ----

    private static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        MoleculeInitializer.initializeMolecule(atomContainer);
        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
    }

    private Map<IRingSet, IAtomContainer> getRingContainerCount(IAtomContainerSet acs) {
        CycleFinder cycles = or(all(), all());
        Map<IRingSet, IAtomContainer> ringSet = new HashMap<>();
        for (IAtomContainer ac : acs.atomContainers()) {
            IRingSet basicRings;
            try {
                Cycles rings = cycles.find(ac);
                basicRings = rings.toRingSet();
                if (!basicRings.isEmpty()) {
                    basicRings.setID(ac.getID());
                    ringSet.put(basicRings, ac);
                }
            } catch (Intractable ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        return ringSet;
    }

    private boolean chipBondInTheRing(IAtomContainer educt, IAtomContainer product) {
        if (sssrEduct.getAtomContainerCount() == 1
                && sssrProduct.getAtomContainerCount() == 1) {
            IAtomContainer ringE = sssrEduct.atomContainers().iterator().next();
            IAtomContainer ringP = sssrProduct.atomContainers().iterator().next();

            if (ringE.getAtomCount() == 6 && ringP.getAtomCount() == 5) {
                return findAndChipBond(ringE, educt);
            } else if (ringE.getAtomCount() == 5 && ringP.getAtomCount() == 6) {
                return findAndChipBond(ringP, product);
            }
        } else if (sssrEduct.getAtomContainerCount() == 2
                && sssrEduct.getAtomContainerCount() == sssrProduct.getAtomContainerCount()) {
            return findAndChipBondBetweenRings(educt) | findAndChipBondBetweenRings(product);
        }
        return false;
    }

    private boolean chipPhophateInSingleReactantProductNotInRing(IAtomContainer educt, IAtomContainer product) throws CDKException, IOException {
        if (ringContainerCountR.isEmpty() && ringContainerCountP.isEmpty()) {
            String phosphateSMILES = "OP(O)(O)=O";
            SmartsPattern smartsPhosphate;

            LOGGER.debug("String phosphateSMILES = \"OP(O)(O)=O\";");

            smartsPhosphate = SmartsPattern.create(phosphateSMILES, SilentChemObjectBuilder.getInstance());
            boolean matchesE = smartsPhosphate.matches(educt);
            boolean matchesP = smartsPhosphate.matches(product);

            LOGGER.debug("String smartsPhosphate Q " + matchesE);
            LOGGER.debug("String smartsPhosphate T " + matchesP);
            if (matchesE && matchesP) {
                boolean findAndChipBondPhophate1 = findAndChipBondPhophate(educt);
                boolean findAndChipBondPhophate2 = findAndChipBondPhophate(product);
                return findAndChipBondPhophate1 && findAndChipBondPhophate2;
            }
        }
        return false;
    }

    private boolean findAndChipBond(IAtomContainer container, IAtomContainer referenceContainer) {
        boolean flag = false;
        if (container != null) {
            for (IBond bond : container.bonds()) {
                if ((bond.getAtom(0).getSymbol().equalsIgnoreCase("O")
                        && bond.getAtom(1).getSymbol().equalsIgnoreCase("C"))
                        || (bond.getAtom(0).getSymbol().equalsIgnoreCase("C")
                        && bond.getAtom(1).getSymbol().equalsIgnoreCase("O"))) {
                    if (!bond.getAtom(0).getFlag(ISAROMATIC)
                            && !bond.getAtom(1).getFlag(ISAROMATIC)) {
                        if (referenceContainer.contains(bond)) {
                            IAtom atom = bond.getAtom(0).getSymbol().equalsIgnoreCase("C") ? bond.getAtom(0) : bond.getAtom(1);
                            List<IBond> neighbourhoodBonds = referenceContainer.getConnectedBondsList(atom);
                            flag = false;
                            for (IBond neighbourhoodBond : neighbourhoodBonds) {
                                if (neighbourhoodBond.contains(atom) && !neighbourhoodBond.getFlag(ISINRING)) {
                                    if ((neighbourhoodBond.getAtom(0).getSymbol().equalsIgnoreCase("O")
                                            && neighbourhoodBond.getAtom(1).getSymbol().equalsIgnoreCase("C"))
                                            || (neighbourhoodBond.getAtom(0).getSymbol().equalsIgnoreCase("C")
                                            && neighbourhoodBond.getAtom(1).getSymbol().equalsIgnoreCase("O"))) {
                                        flag = true;
                                    }
                                }
                            }
                            if (flag) {
                                referenceContainer.removeBond(bond);
                                break;
                            }
                        }
                    }
                }
            }
        }
        return flag;
    }

    private boolean findAndChipBondPhophate(IAtomContainer container) {
        boolean flag = false;
        for (IBond bond : container.bonds()) {
            IAtom atomE = bond.getAtom(0);
            IAtom atomP = bond.getAtom(1);
            if ((atomE.getSymbol().equals("O") && atomP.getSymbol().equals("P"))
                    || (atomE.getSymbol().equals("P") && atomP.getSymbol().equals("O"))) {
                if (bond.getOrder().equals(SINGLE)) {
                    IAtom oxygen = atomE.getSymbol().equals("O") ? atomE : atomP;
                    List<IBond> neighbourBonds = container.getConnectedBondsList(oxygen);
                    if (neighbourBonds.size() == 2) {
                        neighbourBonds.stream().filter((b) -> (b.getAtom(0).getSymbol().equals("O")
                                || b.getAtom(0).getSymbol().equals("P"))).filter((b) -> (b.getAtom(1).getSymbol().equals("O")
                                || b.getAtom(1).getSymbol().equals("P"))).forEach((b) -> {
                            container.removeBond(b);
                        });
                        return true;
                    }
                }
            }
        }
        return flag;
    }

    private boolean findAndChipBondBetweenRings(IAtomContainer container) {
        LOGGER.debug("Find and Chip Bond Between Rings");
        List<IBond> bond_to_be_removed = new ArrayList<>();
        for (IAtom atom : container.atoms()) {
            if (atom.getSymbol().equals("O") && !atom.isAromatic()) {
                Integer ringConns = atom.getProperty(RING_CONNECTIONS);
                int number_of_rings = ringConns != null ? ringConns : 0;
                LOGGER.debug("number_of_rings " + number_of_rings);

                List<IBond> bonds = container.getConnectedBondsList(atom);
                if (bonds.size() == 2 && number_of_rings == 2) {
                    IBond bondToBeChipped = bonds.iterator().next();
                    bond_to_be_removed.add(bondToBeChipped);
                }
            }
        }
        bond_to_be_removed.forEach((bond) -> {
            container.removeBond(bond);
        });
        return !bond_to_be_removed.isEmpty();
    }

    private boolean deleteBonds(IAtomContainer s, IAtomContainer t) throws InvalidSmilesException, CDKException {
        boolean flag = false;
        flag = flag | case1(s, t);
        LOGGER.debug("Case 1: " + flag);
        return flag;
    }

    private boolean case1(IAtomContainer s, IAtomContainer t) throws InvalidSmilesException, CDKException {
        boolean flag = false;
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String moiety = "NC(=O)C1=CN(C=N1)C1OC(COP(O)(O)=O)C(O)C1O";
        IAtomContainer query = smilesParser.parseSmiles(moiety);

        if (isMatch(query, s, false) && isMatch(query, t, false)) {
            IAtomContainer ac1 = s;
            IAtomContainer ac2 = t;

            Map<IAtom, IAtom> subgraph1 = findSubgraph(query, ac1, true, true, true, true);
            Map<IAtom, IAtom> subgraph2 = findSubgraph(query, ac2, true, true, true, true);

            if (subgraph1 != null && subgraph2 != null
                    && subgraph1.isEmpty() && subgraph2.isEmpty()) {
                ac1 = t;
                ac2 = s;
                subgraph1 = findSubgraph(query, ac1, true, true, true, true);
                subgraph2 = findSubgraph(query, ac2, true, true, true, true);
            }

            if (subgraph1 != null && !subgraph1.isEmpty()) {
                Map<IAtom, IAtom> map = subgraph1;
                for (IAtom a : map.values()) {
                    for (IAtom b : map.values()) {
                        if (a != b
                                && a.getAtomicNumber() == 6
                                && b.getAtomicNumber() == 6) {
                            IBond bond = ac1.getBond(a, b);
                            ac1.removeBond(bond);
                            flag = true;
                        }
                    }
                }
            }

            if (subgraph2 != null && !subgraph2.isEmpty()) {
                Map<IAtom, IAtom> map = subgraph2;
                for (IAtom a : map.values()) {
                    for (IAtom b : map.values()) {
                        if (a != b
                                && a.getAtomicNumber() == 6
                                && b.getAtomicNumber() == 6) {
                            IBond bond = ac2.getBond(a, b);
                            ac2.removeBond(bond);
                            flag = true;
                        }
                    }
                }
            }
        }
        return flag;
    }
}
