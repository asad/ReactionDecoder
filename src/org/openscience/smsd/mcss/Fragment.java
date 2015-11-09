/*
 * Copyright (C) 2014 Syed Asad Rahman <asad at ebi.ac.uk>.
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
package org.openscience.smsd.mcss;

import java.io.Serializable;
import java.util.BitSet;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.ShortestPathFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;

/**
 * 
 * 
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public class Fragment implements Comparable<Fragment>, Serializable {

    private static final long serialVersionUID = 134634654886765L;
    private static final Logger LOG = getLogger(Fragment.class.getName());

    /**
     * Return SMILES
     *
     * @param ac
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static String toSmiles(IAtomContainer ac) throws CDKException {
        SmilesGenerator g = new SmilesGenerator().aromatic();
        return g.create(ac);
    }
    private final BitSet fingerprint;
    private final long fingerprintAsLong;
    private final IAtomContainer container;

    /**
     *
     * @param container
     * @throws CDKException
     */
    public Fragment(IAtomContainer container) throws CDKException {
        if (container == null) {
            throw new CDKException("NULL container not supported");
        }
        this.container = container;
        this.fingerprint = new ShortestPathFingerprinter().getBitFingerprint(container).asBitSet();
        this.fingerprintAsLong = convert(this.fingerprint);
    }

    /**
     *
     * @return
     */
    public synchronized IAtomContainer getContainer() {
        return container;
    }

    @Override
    public synchronized boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Fragment other = (Fragment) obj;

        if (this.getContainer() != other.getContainer() && (this.getContainer() == null
                || (this.getContainer().getAtomCount() != other.getContainer().getAtomCount()))) {
            return false;
        }

        if (this.getFingerprint() != other.getFingerprint() && (this.getFingerprint() == null
                || !this.getFingerprint().equals(other.getFingerprint()))) {
            return false;
        }
        return this.fingerprintAsLong == other.fingerprintAsLong;
    }

    @Override
    public synchronized int hashCode() {
        int hash = 3;
        hash = 47 * hash + (this.getFingerprint() != null ? this.getFingerprint().hashCode() : 0);
        hash = 47 * hash + (int) (this.fingerprintAsLong ^ (this.fingerprintAsLong >>> 32));
        return hash;
    }

    @Override
    public synchronized int compareTo(Fragment t) {

        if (this.fingerprintAsLong == t.fingerprintAsLong) {
            return 0;
        } else if (this.fingerprintAsLong > t.fingerprintAsLong) {
            return 1;
        } else {
            return -1;
        }
    }

    private synchronized long convert(BitSet bits) {
        long value = 0L;
        if (bits == null || bits.isEmpty()) {
            return value;
        }
        for (int i = 0; i < bits.length(); ++i) {
            value += bits.get(i) ? (1L << i) : 0L;
        }
        return value;
    }

    /**
     * @return the fingerprint
     */
    public BitSet getFingerprint() {
        return fingerprint;
    }
}
