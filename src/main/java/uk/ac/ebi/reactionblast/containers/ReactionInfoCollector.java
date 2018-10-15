/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.containers;

import java.io.Serializable;
import java.util.BitSet;

import org.openscience.cdk.interfaces.IReaction.Direction;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionInfoCollector implements Serializable {

    private static final long serialVersionUID = 878898788771L;
    private final String sourceDbID;
    private BitSet rorFp = null;
    private BitSet porFp = null;
    private IPatternFingerprinter cfFp = null;
    private IPatternFingerprinter ocFp = null;
    private IPatternFingerprinter stFP = null;
    private final Direction direction;

    /**
     *
     * @param sourceDbID
     * @param direction
     */
    public ReactionInfoCollector(String sourceDbID, Direction direction) {
        this.sourceDbID = sourceDbID;
        this.direction = direction;
    }

    /**
     *
     * @param rorFp
     * @param porFp
     */
    public synchronized void setStructuralFingerprints(BitSet rorFp, BitSet porFp) {
        this.rorFp = rorFp;
        this.porFp = porFp;
    }

    /**
     *
     * @param cfFp
     * @param ocFp
     * @param stFP
     */
    public synchronized void setBondChangeFingerprints(IPatternFingerprinter cfFp, IPatternFingerprinter ocFp, IPatternFingerprinter stFP) {
        this.cfFp = cfFp;
        this.ocFp = ocFp;
        this.stFP = stFP;
    }

    /**
     *
     * @return
     */
    public synchronized BitSet getRorFp() {
        return rorFp;
    }

    /**
     *
     * @return
     */
    public synchronized BitSet getPorFp() {
        return porFp;
    }

    /**
     *
     * @return
     */
    public synchronized IPatternFingerprinter getCfFp() {
        return cfFp;
    }

    /**
     *
     * @return
     */
    public synchronized IPatternFingerprinter getOcFp() {
        return ocFp;
    }

    /**
     *
     * @return
     */
    public synchronized IPatternFingerprinter getStFp() {
        return stFP;
    }

    /**
     *
     * @return
     */
    public synchronized Direction getDirection() {
        return direction;
    }

    /**
     *
     * @return
     */
    public synchronized String getSourceDbID() {
        return sourceDbID;
    }
}
