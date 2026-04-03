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
package com.bioinceptionlabs.reactionblast.cdk;

import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.AtomBondMatcher;
import org.openscience.smsd.BaseMapping;
import org.openscience.smsd.MoleculeInitializer;

import com.bioinceptionlabs.reactionblast.mapping.ReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.mapping.SmsdReactionMappingEngine;
import com.bioinceptionlabs.reactionblast.model.AtomNode;
import com.bioinceptionlabs.reactionblast.model.BondEdge;
import com.bioinceptionlabs.reactionblast.model.ChemToolkit;
import com.bioinceptionlabs.reactionblast.model.MolecularGraph;
import com.bioinceptionlabs.reactionblast.model.ReactionGraph;

/**
 * CDK implementation of ChemToolkit. Provides parsing, serialization,
 * perception, and substructure matching using the Chemistry Development Kit.
 *
 * Usage:
 * <pre>
 *   ChemToolkit.register(new CDKToolkit());
 *   ReactionGraph rxn = ChemToolkit.get().parseReactionSmiles("CC>>CC");
 * </pre>
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class CDKToolkit implements ChemToolkit {

    private static final ReactionMappingEngine MAPPING_ENGINE
            = SmsdReactionMappingEngine.getInstance();

    private final SmilesParser smilesParser;
    private final SmilesGenerator canonicalSmilesGen;
    private final SmilesGenerator mappedSmilesGen;
    private final Aromaticity aromaticity;

    public CDKToolkit() {
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.canonicalSmilesGen = new SmilesGenerator(SmiFlavor.Canonical);
        this.mappedSmilesGen = new SmilesGenerator(
                SmiFlavor.Stereo | SmiFlavor.AtomAtomMap);
        this.aromaticity = new Aromaticity(ElectronDonation.piBonds(),
                Cycles.or(Cycles.all(), Cycles.or(Cycles.relevant(), Cycles.essential())));
    }

    @Override
    public ReactionGraph parseReactionSmiles(String smiles) {
        try {
            IReaction rxn = smilesParser.parseReactionSmiles(smiles);
            return new CDKReactionGraph(rxn);
        } catch (CDKException e) {
            throw new RuntimeException("Failed to parse reaction SMILES: " + smiles, e);
        }
    }

    @Override
    public MolecularGraph parseMoleculeSmiles(String smiles) {
        try {
            IAtomContainer mol = smilesParser.parseSmiles(smiles);
            return new CDKMolecularGraph(mol);
        } catch (CDKException e) {
            throw new RuntimeException("Failed to parse SMILES: " + smiles, e);
        }
    }

    @Override
    public String toSmiles(MolecularGraph mol) {
        try {
            return mappedSmilesGen.create(unwrap(mol));
        } catch (CDKException e) {
            throw new RuntimeException("Failed to generate SMILES", e);
        }
    }

    @Override
    public String toSmiles(ReactionGraph rxn) {
        try {
            return mappedSmilesGen.create(unwrapReaction(rxn));
        } catch (CDKException e) {
            throw new RuntimeException("Failed to generate reaction SMILES", e);
        }
    }

    @Override
    public String toCanonicalSmiles(MolecularGraph mol) {
        try {
            return canonicalSmilesGen.create(unwrap(mol));
        } catch (CDKException e) {
            throw new RuntimeException("Failed to generate canonical SMILES", e);
        }
    }

    @Override
    public void perceiveAtomTypes(MolecularGraph mol) {
        try {
            IAtomContainer ac = unwrap(mol);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        } catch (CDKException e) {
            throw new RuntimeException("Failed to perceive atom types", e);
        }
    }

    @Override
    public void perceiveAromaticity(MolecularGraph mol) {
        try {
            aromaticity.apply(unwrap(mol));
        } catch (CDKException e) {
            throw new RuntimeException("Failed to perceive aromaticity", e);
        }
    }

    @Override
    public void addImplicitHydrogens(MolecularGraph mol) {
        try {
            IAtomContainer ac = unwrap(mol);
            CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(ac);
        } catch (CDKException e) {
            throw new RuntimeException("Failed to add implicit hydrogens", e);
        }
    }

    @Override
    public boolean isSubstructure(MolecularGraph query, MolecularGraph target) {
        try {
            IAtomContainer q = unwrap(query);
            IAtomContainer t = unwrap(target);
            BaseMapping sub = MAPPING_ENGINE.findSubstructure(q, t,
                    AtomBondMatcher.atomMatcher(true, true),
                    AtomBondMatcher.bondMatcher(true, true), true);
            return sub.isSubgraph();
        } catch (CDKException e) {
            return false;
        }
    }

    @Override
    public Map<AtomNode, AtomNode> findMCS(MolecularGraph mol1, MolecularGraph mol2) {
        try {
            IAtomContainer ac1 = unwrap(mol1);
            IAtomContainer ac2 = unwrap(mol2);
            MoleculeInitializer.initializeMolecule(ac1);
            MoleculeInitializer.initializeMolecule(ac2);
            BaseMapping iso = MAPPING_ENGINE.findMcs(ac1, ac2,
                    org.openscience.smsd.BaseMapping.Algorithm.VFLibMCS,
                    AtomBondMatcher.atomMatcher(false, false),
                    AtomBondMatcher.bondMatcher(false, false));
            Map<AtomNode, AtomNode> result = new HashMap<>();
            AtomAtomMapping mapping = iso.getFirstAtomMapping();
            if (mapping != null) {
                for (Map.Entry<IAtom, IAtom> entry : mapping.getMappingsByAtoms().entrySet()) {
                    result.put(new CDKAtomNode(entry.getKey()), new CDKAtomNode(entry.getValue()));
                }
            }
            return result;
        } catch (CDKException e) {
            return new HashMap<>();
        }
    }

    @Override
    public MolecularGraph createMolecularGraph() {
        return new CDKMolecularGraph(new AtomContainer());
    }

    @Override
    public AtomNode createAtomNode(String symbol) {
        IAtom atom = SilentChemObjectBuilder.getInstance().newInstance(IAtom.class, symbol);
        return new CDKAtomNode(atom);
    }

    @Override
    public BondEdge createBondEdge(AtomNode source, AtomNode target, BondEdge.BondOrder order) {
        IAtom a1 = ((CDKAtomNode) source).getCDKAtom();
        IAtom a2 = ((CDKAtomNode) target).getCDKAtom();
        IBond.Order cdkOrder;
        switch (order) {
            case DOUBLE: cdkOrder = IBond.Order.DOUBLE; break;
            case TRIPLE: cdkOrder = IBond.Order.TRIPLE; break;
            case QUADRUPLE: cdkOrder = IBond.Order.QUADRUPLE; break;
            default: cdkOrder = IBond.Order.SINGLE; break;
        }
        return new CDKBondEdge(new Bond(a1, a2, cdkOrder));
    }

    @Override
    public ReactionGraph createReactionGraph() {
        return new CDKReactionGraph(new Reaction());
    }

    // ---- Helper: unwrap graph model back to CDK ----

    private static IAtomContainer unwrap(MolecularGraph mol) {
        if (mol instanceof CDKMolecularGraph) {
            return ((CDKMolecularGraph) mol).getCDKContainer();
        }
        throw new IllegalArgumentException("Expected CDKMolecularGraph, got " + mol.getClass().getName());
    }

    private static IReaction unwrapReaction(ReactionGraph rxn) {
        if (rxn instanceof CDKReactionGraph) {
            return ((CDKReactionGraph) rxn).getCDKReaction();
        }
        throw new IllegalArgumentException("Expected CDKReactionGraph, got " + rxn.getClass().getName());
    }
}
