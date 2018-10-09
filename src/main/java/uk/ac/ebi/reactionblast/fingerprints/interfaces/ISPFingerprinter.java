/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2011-2018       Syed Asad Rahman <asad @ ebi.ac.uk>
 *
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

package uk.ac.ebi.reactionblast.fingerprints.interfaces;

import java.util.BitSet;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * 
 */
public interface ISPFingerprinter extends org.openscience.cdk.fingerprint.IFingerprinter {

    /**
     * Generates a fingerprint of the default fingerprintLength for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @exception CDKException if there is a timeout in ring or aromaticity perception
     * @return A {@link BitSet} representing the fingerprint
     */
    @Override
    IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException;

    /**
     *
     * @param atomContainer
     * @return
     * @throws CDKException
     */
    @Override
    ICountFingerprint getCountFingerprint(IAtomContainer atomContainer) throws CDKException;

    /**
     *
     * @param atomContainer
     * @return
     * @throws CDKException
     */
    @Override
    Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException;

    /**
     *
     * @return
     */
    int getSearchDepth();

    /**
     *
     * @return
     */
    int getSize();

    /**
     * @return the respectFormalCharges
     */
    boolean isRespectFormalCharges();

    /**
     * Should match rings to rings and non-rings to non-rings
     *
     * @return the respect ring matches
     */
    boolean isRespectRingMatches();

    /**
     * @param respectFormalCharges the flag to set if formal charge is checked
     */
    void setRespectFormalCharges(boolean respectFormalCharges);

    /**
     * Ring matches are allowed and non-ring to ring matches are discarded
     *
     * @param respectRingMatches respect the ring-to-ring matches and discard non-ring to ring matches
     */
    void setRespectRingMatches(boolean respectRingMatches);
}
