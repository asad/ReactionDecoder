/* 
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.mcgregor.McGregor;

/**
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib
 * {@link org.openscience.smsd.algorithm.vflib.VF2MCS} and reports
 * MCS between run source and target graphs. Then these solutions are extended
 * using McGregor {@link org.openscience.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class BaseMCS {

    protected int countR;
    protected int countP;
    protected final IAtomContainer source;
    protected final IAtomContainer target;
    private final boolean shouldMatchRings;
    private final boolean matchBonds;
    private final boolean matchAtomType;
    protected final List<Map<IAtom, IAtom>> vfLibSolutions;
    final List<Map<Integer, Integer>> allLocalMCS;
    final List<AtomAtomMapping> allLocalAtomAtomMapping;
    private final static ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(BaseMCS.class);
    private final boolean DEBUG = false;

    BaseMCS(IAtomContainer source, IAtomContainer target, boolean matchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.allLocalAtomAtomMapping = new ArrayList<>();
        this.allLocalMCS = new ArrayList<>();
        this.shouldMatchRings = shouldMatchRings;
        this.matchBonds = matchBonds;
        this.matchAtomType = matchAtomType;
        this.vfLibSolutions = new ArrayList<>();
        this.source = source;
        this.target = target;
    }

    BaseMCS(IQueryAtomContainer source, IAtomContainer target) {
        this.allLocalAtomAtomMapping = new ArrayList<>();
        this.allLocalMCS = new ArrayList<>();
        this.shouldMatchRings = true;
        this.matchBonds = true;
        this.matchAtomType = true;
        this.vfLibSolutions = new ArrayList<>();
        this.source = source;
        this.target = target;
    }

    /**
     *
     * @param cliqueMap
     * @param mapGlobal
     * @return true if condition meet else false
     */
    protected synchronized boolean hasClique(
            Map<Integer, Integer> cliqueMap, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> storedMap : mapGlobal) {
            if (cliqueMap.size() < storedMap.size()) {
                return true;
            } else if (cliqueMap.equals(storedMap)) {
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param cliqueMap
     * @param mapGlobal
     * @return true if condition meet else false
     */
    protected synchronized boolean isCliquePresent(
            Map<Integer, Integer> cliqueMap, List<Map<Integer, Integer>> mapGlobal) {
        return mapGlobal.stream().anyMatch((storedMap) -> (cliqueMap.equals(storedMap)));
    }

    /**
     *
     * @param refinedMCSSeeds
     * @throws CDKException
     * @throws IOException
     */
    protected synchronized void extendCliquesWithMcGregor(
            List<Map<Integer, Integer>> refinedMCSSeeds) throws CDKException, IOException {
        List<List<Integer>> mappings = new ArrayList<>();
        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : refinedMCSSeeds) {
            Map<Integer, Integer> extendMapping = new TreeMap<>(firstPassMappings);
            McGregor mgit;
            if (source instanceof IQueryAtomContainer) {
                mgit = new McGregor((IQueryAtomContainer) source, target, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());
                //Start McGregor search
                mgit.startMcGregorIteration((IQueryAtomContainer) source, mgit.getMCSSize(), extendMapping);
            } else if (countR > countP) {
                mgit = new McGregor(source, target, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());

                //Start McGregor search
                mgit.startMcGregorIteration(source, mgit.getMCSSize(), extendMapping);
            } else {
                extendMapping.clear();
                mgit = new McGregor(target, source, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());
                ROPFlag = false;
                firstPassMappings.entrySet().stream().forEach((map) -> {
                    extendMapping.put(map.getValue(), map.getKey());
                });
                //Start McGregor search
                mgit.startMcGregorIteration(target, mgit.getMCSSize(), extendMapping);
            }
            mappings = mgit.getMappings();
        }
//        System.out.println("\nSol count after MG " + mappings.size());
        setMcGregorMappings(ROPFlag, mappings);
//        System.out.println("After set Sol count MG" + allMCS.size());
//        System.out.println("MCSSize " + vfMCSSize + "\n");
    }

    /**
     *
     * @param RONP
     */
    protected synchronized void setVFMappings(boolean RONP) {
//        System.out.println(" setVFMappings ");
        /*
         * Sort biggest clique to smallest
         */
        Collections.sort(vfLibSolutions, new Map2ValueComparator(SortOrder.DESCENDING));
        for (Map<IAtom, IAtom> solution : vfLibSolutions) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<>();
            solution.entrySet().stream().forEach((mapping) -> {
                IAtom qAtom;
                IAtom tAtom;
                Integer qIndex;
                Integer tIndex;

                if (RONP) {
                    qAtom = mapping.getKey();
                    tAtom = mapping.getValue();
                    qIndex = source.indexOf(qAtom);
                    tIndex = target.indexOf(tAtom);
                } else {
                    tAtom = mapping.getKey();
                    qAtom = mapping.getValue();
                    qIndex = source.indexOf(qAtom);
                    tIndex = target.indexOf(tAtom);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            });

            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())) {
                getLocalAtomMCSSolution().add(atomatomMapping);
                getLocalMCSSolution().add(indexindexMapping);
            }
        }

        if (DEBUG) {
            System.out.println("VF seed mappings stored count: " + getLocalMCSSolution().size());
        }
    }

    private synchronized void setMcGregorMappings(boolean RONP,
            List<List<Integer>> mappings) throws CDKException {
        int counter = 0;
        int solSize = 0;
        getLocalAtomMCSSolution().clear();
        getLocalMCSSolution().clear();
        for (List<Integer> mapping : mappings) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<>();
            for (int index = 0; index < mapping.size(); index += 2) {
                IAtom qAtom;
                IAtom tAtom;
                int qIndex;
                int tIndex;

                if (RONP) {
                    qAtom = source.getAtom(mapping.get(index));
                    tAtom = target.getAtom(mapping.get(index + 1));
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qAtom = source.getAtom(mapping.get(index + 1));
                    tAtom = target.getAtom(mapping.get(index));
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    throw new CDKException("Atom index pointing to NULL");
                }
            }
            if (indexindexMapping.size() > solSize) {
                solSize = indexindexMapping.size();
                getLocalAtomMCSSolution().clear();
                getLocalMCSSolution().clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())
                    && indexindexMapping.size() == solSize) {
                getLocalAtomMCSSolution().add(counter, atomatomMapping);
                getLocalMCSSolution().add(counter, indexindexMapping);
                counter++;
            }
        }

    }

    protected synchronized IAtomContainer getReactantMol() {
        return source;
    }

    protected synchronized IAtomContainer getProductMol() {
        return target;
    }

    /**
     * @return the shouldMatchRings
     */
    protected boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the shouldMatchBonds
     */
    protected synchronized boolean isBondMatchFlag() {
        return matchBonds;
    }

    /**
     * @return the allLocalMCS
     */
    private synchronized List<Map<Integer, Integer>> getLocalMCSSolution() {
        return Collections.synchronizedList(allLocalMCS);
    }

    /**
     * @return the allLocalAtomAtomMapping
     */
    private synchronized List<AtomAtomMapping> getLocalAtomMCSSolution() {
        return Collections.synchronizedList(allLocalAtomAtomMapping);
    }

    protected synchronized boolean isExtensionRequired(List<Map<Integer, Integer>> mcsSeeds) {
        int maxSize = 0;
        for (Map<Integer, Integer> map : mcsSeeds) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return this.source.getAtomCount() > maxSize && this.target.getAtomCount() > maxSize;
    }

    protected synchronized boolean isExtensionRequired() {
        int commonAtomCount = checkCommonAtomCount(getReactantMol(), getProductMol());
        int maxSize = 0;
        for (Map<Integer, Integer> map : allLocalMCS) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return commonAtomCount > maxSize;
    }

    private synchronized int checkCommonAtomCount(
            IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> atoms = new ArrayList<>();
        for (int i = 0; i < reactantMolecule.getAtomCount(); i++) {
            atoms.add(reactantMolecule.getAtom(i).getSymbol());
        }
        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {
            String symbol = productMolecule.getAtom(i).getSymbol();
            if (atoms.contains(symbol)) {
                atoms.remove(symbol);
                common++;
            }
        }
        return common;
    }

    /**
     * @return the matchAtomType
     */
    public boolean isMatchAtomType() {
        return matchAtomType;
    }
}
