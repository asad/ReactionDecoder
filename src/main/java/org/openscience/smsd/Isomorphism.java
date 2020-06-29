/* 
 * Copyright (C) 2009-2020  Syed Asad Rahman <asad at ebi.ac.uk>
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

import java.io.Serializable;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.algorithm.rgraph.CDKMCSHandler;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.smsd.algorithm.ventofoggia.VF2MCS;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.interfaces.IResults;

/**
 * <p>
 * This class implements the Isomorphism- a multipurpose structure comparison
 * tool. It allows users to, i) find the maximal common substructure(s) (MCS);
 * ii) perform the mapping of a substructure in another structure, and; iii) map
 * two isomorphic structures.</p>
 *
 * <p>
 * It also comes with various published algorithms. The user is free to choose
 * his favorite algorithm to perform MCS or substructure search. For
 * example:</p> <OL> <lI>0: Default, <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3:
 * CDKMCS </OL>
 * <p>
 * It also has a set of robust chemical filters (i.e. bond energy, fragment
 * count, stereo & bond match) to sort the reported MCS solutions in a
 * chemically relevant manner. Each comparison can be made with or without using
 * the bond sensitive mode and with implicit or explicit hydrogens.</p>
 *
 * <p>
 * If you are using <font color="#FF0000">Isomorphism, please cite Rahman
 * <i>et.al. 2009</i></font> {
 *
 * @cdk.cite SMSD2009}. The Isomorphism algorithm is described in this paper.
 * </p>
 *
 * <p>
 * An example for <b>MCS search</b>:</p> <font color="#003366">  <pre>
 *
 *
 * SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
 * // Benzene
 * IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
 * // Napthalene
 * IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
 * //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
 * //Algorithm is VF2MCS
 * //Bond Sensitive is set True
 * //Ring Match is set True
 * Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
 * // set chemical filter true
 * comparison.setChemFilters(true, true, true);
 * //Get similarity score
 * System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 * Assert.assertEquals(0.6, comparison.getTanimotoSimilarity());
 * Assert.assertEquals(12, comparison.getAllAtomMapping().size());
 * // Print the mapping between molecules
 * System.out.println(" Mappings: ");
 * for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
 *      for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappingsByAtoms().entrySet()) {
 *          IAtom sourceAtom = mapping.getKey();
 *          IAtom targetAtom = mapping.getValue();
 *          System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
 *          System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
 *      }
 *      System.out.println("");
 *  }
 *
 *
 * </pre> </font>
 *
 * java1.8+
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 *
 */
public final class Isomorphism extends BaseMapping implements Serializable {

    private final static boolean DEBUG = false;
    private final static ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Isomorphism.class);
    static final long serialVersionUID = 0x24845e5c5ae877L;
    private final Algorithm algorithmType;
    private double bondSensitiveMcGregorOut = -1;//mins
    private double bondInSensitiveMcGregor = -1;//mins

    /**
     *
     * IMP: Initialize query and target molecules
     * (MoleculeInitializer.initializeMolecule).Note: Here its assumed that
     * hydrogens are implicit and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector before
     * initializing calling this method.Please call MoleculeInitializer before
     * calling substructure search
     *
     * if (super.isMatchRings()) { try {
     * MoleculeInitializer.initializeMolecule(super.getQuery());
     * MoleculeInitializer.initializeMolecule(super.getTarget()); } catch
     * (CDKException ex) { } }
     *
     * Initialize query and target molecules.Note: Here its assumed that
     * hydrogens are implicit and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector before
     * initializing calling this method.
     *
     *
     *
     *
     * @param query query mol
     * @param target target mol This is the algorithm factory and entry port for
     * all the MCS algorithm in the Isomorphism supported algorithm
     * {@link org.openscience.smsd.interfaces.Algorithm} types: <OL> <lI>0:
     * Default,
     * <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3: CDKMCS </OL>
     * @param algorithmType {@link org.openscience.smsd.interfaces.Algorithm}
     * @param am
     * @param bm
     * @throws org.openscience.cdk.exception.CDKException
     */
    public Isomorphism(
            IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            AtomMatcher am, BondMatcher bm) throws CDKException {
        super(query, target, am, bm);
        this.algorithmType = algorithmType;
        mcsBuilder(super.getQuery(), super.getTarget());
        super.setSubgraph(isSubgraph());
    }

    /**
     *
     * IMP: Initialize query and target molecules
     * (MoleculeInitializer.initializeMolecule).Note: Here its assumed that
     * hydrogens are implicit and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector before
     * initializing calling this method.Please call MoleculeInitializer before
     * calling substructure search
     *
     * if (super.isMatchRings()) { try {
     * MoleculeInitializer.initializeMolecule(super.getQuery());
     * MoleculeInitializer.initializeMolecule(super.getTarget()); } catch
     * (CDKException ex) { } }
     *
     * Initialize query and target molecules.Note: Here its assumed that
     * hydrogens are implicit and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector before
     * initializing calling this method.
     *
     *
     *
     *
     * @param query query mol
     * @param target target mol This is the algorithm factory and entry port for
     * all the MCS algorithm in the Isomorphism supported algorithm
     * {@link org.openscience.smsd.interfaces.Algorithm} types: <OL> <lI>0:
     * Default,
     * <lI>1: MCSPlus, <lI>2: VFLibMCS, <lI>3: CDKMCS </OL>
     * @param algorithmType {@link org.openscience.smsd.interfaces.Algorithm}
     * @throws org.openscience.cdk.exception.CDKException
     */
    public Isomorphism(
            IQueryAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType) throws CDKException {
        super(query, target, AtomMatcher.forQuery(), BondMatcher.forQuery());
        this.algorithmType = algorithmType;
        mcsBuilder(super.getQuery(), super.getTarget());
        super.setSubgraph(isSubgraph());
    }

    private synchronized void mcsBuilder(IAtomContainer mol1, IAtomContainer mol2) throws CDKException {
        int rBondCount = mol1.getBondCount();
        int pBondCount = mol2.getBondCount();

        int rAtomCount = mol1.getAtomCount();
        int pAtomCount = mol2.getAtomCount();

        int expectedMaxGraphmatch = expectedMaxGraphmatch(mol1, mol2);

        if (DEBUG) {
            System.out.println("Expected match: " + expectedMaxGraphmatch);
            System.out.println("rAtomCount " + rAtomCount + ", rBondCount " + rBondCount);
            System.out.println("pAtomCount " + pAtomCount + ", pBondCount " + pBondCount);
        }
        if (expectedMaxGraphmatch == 1 || rAtomCount == 1 || pAtomCount == 1) {
            singleMapping();
        } else {
            chooseAlgorithm();
        }
    }

    private synchronized void chooseAlgorithm() throws CDKException {

        switch (algorithmType) {
            case CDKMCS:
                if (DEBUG) {
                    System.out.println("Calling CDKMCS ");
                }
                cdkMCSAlgorithm();
                if (DEBUG) {
                    System.out.println("Calling DONE CDKMCS ");
                }
                break;
            case DEFAULT:
                if (DEBUG) {
                    System.out.println("Calling DEFAULT ");
                }
                defaultMCSAlgorithm();
                if (DEBUG) {
                    System.out.println("Calling DONE DEFAULT ");
                }
                break;
            case MCSPlus:
                if (DEBUG) {
                    System.out.println("Calling MCSPlus ");
                }
                mcsPlusAlgorithm();
                if (DEBUG) {
                    System.out.println("Calling DONE MCSPlus ");
                }
                break;
            case VFLibMCS:
                if (DEBUG) {
                    System.out.println("Calling VFLibMCS ");
                }
                vfLibMCSAlgorithm();
                if (DEBUG) {
                    System.out.println("Calling DONE VFLibMCS ");
                }
                break;
        }

    }

    private synchronized boolean cdkMCSAlgorithm() {
        CDKMCSHandler mcs;
        if (getQuery() instanceof IQueryAtomContainer) {
            mcs = new CDKMCSHandler(getQuery(), getTarget(), atomMatcher, bondMatcher);
        } else {
            mcs = new CDKMCSHandler(getQuery(), getTarget(), atomMatcher, bondMatcher);
        }
        clearMaps();
        getMCSList().addAll(mcs.getAllAtomMapping());
        return mcs.isTimeout();
    }

    private synchronized boolean mcsPlusAlgorithm() throws CDKException {
        IResults mcs;
        int expectedMaxGraphmatch = expectedMaxGraphmatch(getQuery(), getTarget());
        if (getQuery() instanceof IQueryAtomContainer) {
            if (DEBUG) {
                System.out.println("org.openscience.smsd.algorithm.mcsplus2.MCSPlusMapper");
            }
            mcs = new org.openscience.smsd.algorithm.mcsplus2.MCSPlusMapper((IQueryAtomContainer) getQuery(), getTarget(), atomMatcher, bondMatcher);
        } else if (expectedMaxGraphmatch < 3) {
            if (DEBUG) {
                System.out.println("org.openscience.smsd.algorithm.mcsplus1.MCSPlusMapper");
            }
            mcs = new org.openscience.smsd.algorithm.mcsplus1.MCSPlusMapper(getQuery(), getTarget(), atomMatcher, bondMatcher);
        } else if (expectedMaxGraphmatch > 3) {
            if (DEBUG) {
                System.out.println("org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper");
            }
            mcs = new org.openscience.smsd.algorithm.mcsplus.MCSPlusMapper(getQuery(), getTarget(), atomMatcher, bondMatcher);
        } else {
            if (DEBUG) {
                System.out.println("org.openscience.smsd.algorithm.mcsplus2.MCSPlusMapper");
            }
            mcs = new org.openscience.smsd.algorithm.mcsplus2.MCSPlusMapper(getQuery(), getTarget(), atomMatcher, bondMatcher);
        }
        clearMaps();
        getMCSList().addAll(mcs.getAllAtomMapping());
        return false;
    }

    private synchronized boolean substructureAlgorithm() throws CDKException {
        if (DEBUG) {
            System.out.println("Check substructureAlgorithm");
        }
        Substructure mcs;
        if (getQuery() instanceof IQueryAtomContainer) {
            mcs = new Substructure((IQueryAtomContainer) getQuery(), getTarget(), atomMatcher, bondMatcher, true);
        } else {
            mcs = new Substructure(getQuery(), getTarget(), atomMatcher, bondMatcher, true);
        }
        clearMaps();
        if (mcs.isSubgraph()) {
            getMCSList().addAll(mcs.getAllAtomMapping());
        }
        return mcs.isSubgraph();
    }

    private synchronized void vfLibMCSAlgorithm() throws CDKException {
        VF2MCS mcs;
        if (getQuery() instanceof IQueryAtomContainer) {
            mcs = new VF2MCS((IQueryAtomContainer) getQuery(), getTarget(), atomMatcher, bondMatcher);
        } else {
            mcs = new VF2MCS(getQuery(), getTarget(), atomMatcher, bondMatcher);
        }
        clearMaps();
        getMCSList().addAll(mcs.getAllAtomMapping());
    }

    private synchronized void singleMapping() {
        SingleMappingHandler mcs;
        mcs = new SingleMappingHandler(getQuery(), getTarget(), atomMatcher);
        clearMaps();
        getMCSList().addAll(mcs.getAllAtomMapping());
    }

    private synchronized void defaultMCSAlgorithm() {
        try {
            boolean substructureAlgorithm = false;

            if (DEBUG) {
                System.out.println("defaultMCSAlgorithm - substructure check ");
            }
            if (getQuery().getAtomCount() > 1
                    && getTarget().getAtomCount() > 1) {
                substructureAlgorithm = substructureAlgorithm();
            }
            if (DEBUG) {
                System.out.println("defaultMCSAlgorithm - no substructure ");
            }
            if (!substructureAlgorithm) {
                int expectedMaxGraphmatch = expectedMaxGraphmatch(getQuery(), getTarget());
                boolean moleculeConnected = isMoleculeConnected(getQuery(), getTarget());
                if (DEBUG) {
                    System.out.println("Expected Match Size: " + expectedMaxGraphmatch);

                }
                if (moleculeConnected
                        && getQuery().getBondCount() > 1
                        && getTarget().getBondCount() > 1) {
                    if (DEBUG) {
                        System.out.println("defaultMCSAlgorithm - Calling CDKMCS ");
                    }
                    cdkMCSAlgorithm();
                    if (DEBUG) {
                        System.out.println("getFirstAtomMapping().getCount() " + getFirstAtomMapping().getCount());
                        System.out.println("defaultMCSAlgorithm - Done CDKMCS ");
                    }
                } else {
                    if (DEBUG) {
                        System.out.println("defaultMCSAlgorithm - Calling MCSPlus ");
                    }
                    mcsPlusAlgorithm();
                    if (DEBUG) {
                        System.out.println("getFirstAtomMapping().getCount() " + getFirstAtomMapping().getCount());
                        System.out.println("defaultMCSAlgorithm - - Calling MCSPlus ");
                    }
                }
            }
        } catch (CDKException e) {
            LOGGER.error(Level.SEVERE, null, e);
        }
    }

    /**
     *
     * @return true if query is a subgraph of the target
     */
    @Override
    public synchronized boolean isSubgraph() {

        float mappingSize;
        if (getMappingCount() > 0) {
            mappingSize = getAllAtomMapping().iterator().next().getCount();
        } else {
            return false;
        }
        int sourceAtomCount = getQuery().getAtomCount();
        int targetAtomCount = getTarget().getAtomCount();

        if (mappingSize == sourceAtomCount && mappingSize <= targetAtomCount) {
            if (mappingSize == 1) {
                return true;
            } else if (!getAllBondMaps().isEmpty()
                    && getAllBondMaps().iterator().next().size() == getQuery().getBondCount()) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the bondSensitiveMcGregorOut
     */
    public double getBondSensitiveMcGregorOut() {
        return bondSensitiveMcGregorOut;
    }

    /**
     * @param bondSensitiveMcGregorOut the bondSensitiveMcGregorOut to set
     */
    public void setBondSenSitiveMcGregorOut(double bondSensitiveMcGregorOut) {
        this.bondSensitiveMcGregorOut = bondSensitiveMcGregorOut;
    }

    /**
     * @return the bondInSensitiveMcGregor
     */
    public double getBondInSensitiveMcGregor() {
        return bondInSensitiveMcGregor;
    }

    /**
     * @param bondInSensitiveMcGregor the bondInSensitiveMcGregor to set
     */
    public void setBondInSenSitiveMcGregor(double bondInSensitiveMcGregor) {
        this.bondInSensitiveMcGregor = bondInSensitiveMcGregor;
    }
}
