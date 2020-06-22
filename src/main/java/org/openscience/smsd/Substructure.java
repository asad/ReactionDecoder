/* Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.smsd.algorithm.ventofoggia.VF2Substructure;

/**
 * This is an ultra fast method to report if query is a substructure for target
 * molecule. If this case is true then it returns only all mapping.
 *
 *
 * Please call MoleculeInitializer before calling substructure search
 *
 * if (super.isMatchRings()) { try {
 * MoleculeInitializer.initializeMolecule(super.getQuery());
 * MoleculeInitializer.initializeMolecule(super.getTarget()); } catch
 * (CDKException ex) { } }
 *
 * This is much faster than {@link
 * org.openscience.smsd.algorithm.vflib.substructure} class as it only reports
 * first match and backtracks.
 *
 * This class should only be used to report if a query graph is a substructure
 * of the target graph.
 *
 *  *
 * <p>
 * An example for <b>Substructure search</b>:</p> <font color="#003366">
 * <pre>
 * SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
 * IAtomContainer query = sp.parseSmiles("CC");
 * IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
 * Substructure smsd = new Substructure(query, target, true, false, true, true);
 * Assert.assertTrue(smsd.isSubgraph());
 * Assert.assertEquals(18, smsd.getAllAtomMapping().size());
 *
 * IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
 * smsd = new Substructure(queryContainer, target, true);
 * Assert.assertTrue(smsd.isSubgraph());
 *
 * </pre> </font>
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public final class Substructure extends BaseMapping {

    private final static boolean DEBUG = false;
    private int vfMappingSize = -1;
    private final ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(Substructure.class);

    /**
     * Constructor for VF Substructure Algorithm
     *
     * Please call before calling substructure search
     *
     * if (super.isMatchRings()) { try {
     * MoleculeInitializer.initializeMolecule(super.getQuery());
     * MoleculeInitializer.initializeMolecule(super.getTarget()); } catch
     * (CDKException ex) { } }
     *
     * @param query
     * @param target
     * @param am
     * @param bm
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IAtomContainer query,
            IAtomContainer target,
            AtomMatcher am,
            BondMatcher bm,
            boolean findAllSubgraph) throws CDKException {
        super(query, target, am, bm);
        super.setSubgraph(findSubgraphs(findAllSubgraph));
    }

    /**
     * Constructor for VF Substructure Algorithm
     *
     * @param query
     * @param target
     * @param findAllSubgraphFlag report all subgraphs
     * @throws CDKException
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllSubgraphFlag) throws CDKException {
        super(query, target, AtomMatcher.forQuery(), BondMatcher.forQuery());
        super.setSubgraph(findSubgraphs(findAllSubgraphFlag));
    }

    private synchronized boolean hasMap(AtomAtomMapping map, List<AtomAtomMapping> mapGlobal) {
        return mapGlobal.stream().anyMatch((test) -> (test.equals(map)));
    }

    /**
     * Returns true if query is a subgraph of target molecule
     *
     * @return
     * @throws CDKException
     */
    private synchronized boolean findSubgraphs(boolean findAllMatches) throws CDKException {
        boolean isSubgraph;

        if ((getTarget() == null) || (getQuery() == null)) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (getQuery().getAtomCount() > getTarget().getAtomCount()) {
            return false;
        }

        int rBondCount = getQuery().getBondCount();
        int pBondCount = getTarget().getBondCount();

        int rAtomCount = getQuery().getAtomCount();
        int pAtomCount = getTarget().getAtomCount();

        int expectedMaxGraphmatch = expectedMaxGraphmatch(getQuery(), getTarget());

        if (DEBUG) {
            System.out.println("Expected match: " + expectedMaxGraphmatch);
            System.out.println("rAtomCount " + rAtomCount + ", rBondCount " + rBondCount);
            System.out.println("pAtomCount " + pAtomCount + ", pBondCount " + pBondCount);
        }
        if (expectedMaxGraphmatch == 1 && rAtomCount <= pAtomCount) {
            isSubgraph = singleMapping();
        } else {
            List<AtomAtomMapping> mappingsVF2 = new ArrayList<>();
            VF2Substructure mapper;
            if (getQuery() instanceof IQueryAtomContainer) {
                mapper = new VF2Substructure((IQueryAtomContainer) getQuery(), getTarget(), findAllMatches, atomMatcher, bondMatcher);
            } else {
                if (DEBUG) {
                    System.out.println("calling VF2Sub");
                }
                mapper = new VF2Substructure(getQuery(), getTarget(), atomMatcher, bondMatcher, findAllMatches);
                if (DEBUG) {
                    System.out.println("done calling VF2Sub");
                }
            }
            isSubgraph = mapper.isSubgraph();
            List<AtomAtomMapping> atomMappings = mapper.getAllAtomMapping();
            if (isSubgraph) {
                mappingsVF2.addAll(atomMappings);
            } else {
                return false;
            }
            setVFMappings(mappingsVF2);
        }
        return isSubgraph;
    }

    private synchronized void setVFMappings(List<AtomAtomMapping> mappingsVF2) {
        int counter = 0;
        for (AtomAtomMapping solution : mappingsVF2) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(getQuery(), getTarget());
            if (solution.getCount() > vfMappingSize) {
                this.vfMappingSize = solution.getCount();
                counter = 0;
            }
            solution.getMappingsByAtoms().entrySet().stream().forEach((mapping) -> {
                IAtom qAtom;
                IAtom tAtom;

                qAtom = mapping.getKey();
                tAtom = mapping.getValue();

                if (qAtom != null && tAtom != null) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            });
            if (!atomatomMapping.isEmpty() && !hasMap(atomatomMapping, getMCSList())
                    && atomatomMapping.getCount() == vfMappingSize) {
                getMCSList().add(counter, atomatomMapping);
                counter++;
            }
        }
    }

    private synchronized boolean singleMapping() {
        SingleMappingHandler mcs;
        if (!(getQuery() instanceof IQueryAtomContainer) && !(getTarget() instanceof IQueryAtomContainer)) {
            mcs = new SingleMappingHandler(getQuery(), getTarget(), atomMatcher);
        } else {
            mcs = new SingleMappingHandler((IQueryAtomContainer) getQuery(), getTarget(), atomMatcher);
        }
        return mcs.getAllAtomMapping() != null && !mcs.getAllAtomMapping().isEmpty();
    }
}
