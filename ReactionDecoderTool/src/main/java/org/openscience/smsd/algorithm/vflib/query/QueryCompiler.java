/* Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
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
 *
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.openscience.smsd.algorithm.vflib.query;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomTypeMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultBondMatcher;
import org.openscience.smsd.algorithm.vflib.builder.VFQueryBuilder;
import org.openscience.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.smsd.algorithm.vflib.interfaces.IQueryCompiler;

/**
 * This class creates an template for MCS/substructure query. 
 * 
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class QueryCompiler implements IQueryCompiler {
    private static final Logger LOG = getLogger(QueryCompiler.class.getName());

    private IAtomContainer molecule = null;
    private IQueryAtomContainer queryMolecule = null;
    private final boolean shouldMatchBonds;
    private final boolean shouldMatchRings;
    private final boolean matchAtomType;

    /**
     * Construct query object from the molecule
     *
     * @param molecule
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public QueryCompiler(IAtomContainer molecule, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.setMolecule(molecule);
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomType = matchAtomType;
        this.shouldMatchBonds = shouldMatchBonds;
    }

    /**
     * Construct query object from the molecule
     *
     * @param molecule
     */
    public QueryCompiler(IQueryAtomContainer molecule) {
        this.setQueryMolecule(molecule);
        this.shouldMatchRings = true;
        this.matchAtomType = true;
        this.shouldMatchBonds = true;
    }

    /**
     * Set GraphMolecule
     *
     * @param molecule
     */
    private synchronized void setMolecule(IAtomContainer molecule) {
        this.molecule = molecule;
    }

    /**
     * Set GraphMolecule
     *
     * @param molecule
     */
    private synchronized void setQueryMolecule(IQueryAtomContainer molecule) {
        this.queryMolecule = molecule;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized IQuery compile() {
        return this.queryMolecule == null ? build(molecule) : build(queryMolecule);
    }

    private synchronized IQuery build(IAtomContainer queryMolecule) {
        VFQueryBuilder result = new VFQueryBuilder();
        for (IAtom atom : queryMolecule.atoms()) {
            AtomMatcher matcher = createAtomMatcher(atom);
            if (matcher != null) {
                result.addNode(matcher, atom);
            }
        }
        for (int i = 0; i < queryMolecule.getBondCount(); i++) {
            IBond bond = queryMolecule.getBond(i);
            IAtom atomI = bond.getAtom(0);
            IAtom atomJ = bond.getAtom(1);
            result.connect(result.getNode(atomI), result.getNode(atomJ), createBondMatcher(bond));
        }
        return result;
    }

    private synchronized AtomMatcher createAtomMatcher(IAtom atom) {
        if (isMatchAtomType()) {
            return new DefaultAtomTypeMatcher(atom, isShouldMatchRings());
        } else {
            return new DefaultAtomMatcher(atom, isShouldMatchRings());
        }
    }

    private synchronized BondMatcher createBondMatcher(IBond bond) {
        return new DefaultBondMatcher(bond, isBondMatchFlag(), isShouldMatchRings(), isMatchAtomType());
    }

    /**
     * @return the shouldMatchBonds
     */
    private synchronized boolean isBondMatchFlag() {
        return shouldMatchBonds;
    }

    /**
     * @return the shouldMatchRings
     */
    private synchronized boolean isShouldMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the matchAtomType
     */
    private boolean isMatchAtomType() {
        return matchAtomType;
    }
}
