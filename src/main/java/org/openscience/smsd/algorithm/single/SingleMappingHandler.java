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
 * MERCHANTABILITY or FITNESS FOR sourceAtom PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.single;

import java.util.*;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.interfaces.IResults;

/**
 * This is a handler class for single atom mapping
 * ({@link org.openscience.smsd.algorithm.single.SingleMapping}).
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class SingleMappingHandler implements IResults {

    private final ILoggingTool Logger
            = LoggingToolFactory.createLoggingTool(SingleMappingHandler.class);
    private List<AtomAtomMapping> allAtomMCS = null;
    private final IAtomContainer source;
    private final IAtomContainer target;
    private final AtomMatcher atomMatcher;

    /**
     *
     * @param source
     * @param target
     * @param am
     */
    public SingleMappingHandler(
            IAtomContainer source,
            IAtomContainer target,
            AtomMatcher am) {
        allAtomMCS = new ArrayList<>();
        this.source = source;
        this.target = target;
        this.atomMatcher = am;
        searchMCS();
    }

    /**
     *
     * @param source
     * @param target
     * @param am
     */
    public SingleMappingHandler(
            IQueryAtomContainer source,
            IAtomContainer target,
            AtomMatcher am) {
        allAtomMCS = new ArrayList<>();
        this.source = source;
        this.target = target;
        this.atomMatcher = am;
        searchMCS();
    }

    /**
     * Function is called by the main program and serves as a starting point for
     * the comparison procedure. {@inheritDoc}
     *
     */
    private synchronized void searchMCS() {
        SingleMapping singleMapping = new SingleMapping();
        List<Map<IAtom, IAtom>> mappings = null;
        try {
            if (target instanceof IQueryAtomContainer) {
                throw new CDKException("Target can't be IQueryAtomContainer");
            } else if (!(source instanceof IQueryAtomContainer)) {
                mappings = singleMapping.getOverLaps(source, target, atomMatcher);
            } else {
                mappings = singleMapping.getOverLaps((IQueryAtomContainer) source, target, atomMatcher);
            }
        } catch (CDKException ex) {
            Logger.error(Level.SEVERE, null, ex);
        }
        setAllAtomMapping(mappings);
        //setStereoScore();
    }

    private synchronized void setAllAtomMapping(List<Map<IAtom, IAtom>> mappings) {

        try {
            int counter = 0;
            for (Map<IAtom, IAtom> solution : mappings) {
                AtomAtomMapping atomMappings = new AtomAtomMapping(source, target);
                solution.entrySet().stream().forEach((map) -> {
                    IAtom sourceAtom = map.getKey();
                    IAtom targetAtom = map.getValue();
                    atomMappings.put(sourceAtom, targetAtom);
                });
                allAtomMCS.add(counter, atomMappings);
                counter++;
            }
        } catch (Exception I) {
            I.getCause();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }
}
