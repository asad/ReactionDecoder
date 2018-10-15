/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.container.helper;

import java.io.Serializable;
import java.util.Objects;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MolMapping extends Object implements Serializable {

    private static final long serialVersionUID = 1738327023703717L;
    private final String mol1;
    private final String mol2;
    private final Integer indexI;
    private final Integer indexJ;
    private Integer indexStep = 0;
    private boolean keggMapping;
    private boolean rBLASTMapping;
    private IAtomContainer matchedMol = null;
    private String matchedSMILES = null;

    /**
     *
     * @param mol1
     * @param mol2
     * @param indexI
     * @param indexJ
     */
    public MolMapping(String mol1, String mol2, Integer indexI, Integer indexJ) {
        this.mol1 = mol1;
        this.mol2 = mol2;
        this.indexI = indexI;
        this.indexJ = indexJ;
        this.keggMapping = false;
        this.rBLASTMapping = false;
    }

    /**
     *
     * @return
     */
    public String getQuery() {
        return mol1;
    }

    /**
     *
     * @return
     */
    public String getTarget() {
        return mol2;
    }

    /**
     * @return the keggMapping
     */
    public synchronized boolean isKeggMapping() {
        return keggMapping;
    }

    /**
     * @param keggMapping the keggMapping to set
     */
    public synchronized void setKeggMapping(boolean keggMapping) {
        this.keggMapping = keggMapping;
    }

    /**
     * @return the rBLASTMapping
     */
    public synchronized boolean isrBLASTMapping() {
        return rBLASTMapping;
    }

    /**
     * @param rBLASTMapping the rBLASTMapping to set
     */
    public synchronized void setReactionMapping(boolean rBLASTMapping) {
        this.rBLASTMapping = rBLASTMapping;
    }

    /**
     * @return the indexI
     */
    public synchronized Integer getIndexI() {
        return indexI;
    }

    /**
     * @return the indexJ
     */
    public synchronized Integer getIndexJ() {
        return indexJ;
    }

    /**
     * @return the matchedMol
     */
    public synchronized IAtomContainer getMatchedMol() {
        return matchedMol;
    }

    /**
     * @param matchedMol the matchedMol to set
     */
    public synchronized void setMatchedMol(IAtomContainer matchedMol) {
        this.matchedMol = matchedMol;
    }

    /**
     * @return the matchedSMILES
     */
    public synchronized String getMatchedSMILES() {
        return matchedSMILES;
    }

    /**
     * @param matchedSMILES the matchedSMILES to set
     * @param step index for the mapping process
     */
    public synchronized void setMatchedSMILES(String matchedSMILES, Integer step) {
        this.matchedSMILES = matchedSMILES;
        setIndexStep(step);
    }

    /**
     * @return the indexStep
     */
    public synchronized Integer getIndexStep() {
        return indexStep;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MolMapping other = (MolMapping) obj;
        if ((this.mol1 == null) ? (other.mol1 != null) : !this.mol1.equals(other.mol1)) {
            return false;
        }
        if ((this.mol2 == null) ? (other.mol2 != null) : !this.mol2.equals(other.mol2)) {
            return false;
        }
        if (!Objects.equals(this.indexI, other.indexI) && (this.indexI == null || !this.indexI.equals(other.indexI))) {
            return false;
        }
        if (!Objects.equals(this.indexJ, other.indexJ) && (this.indexJ == null || !this.indexJ.equals(other.indexJ))) {
            return false;
        }
        if (!Objects.equals(this.indexStep, other.indexStep) && (this.indexStep == null || !this.indexStep.equals(other.indexStep))) {
            return false;
        }
        if (this.keggMapping != other.keggMapping) {
            return false;
        }
        if (this.rBLASTMapping != other.rBLASTMapping) {
            return false;
        }
        if (this.matchedMol != other.matchedMol && (this.matchedMol == null || !this.matchedMol.equals(other.matchedMol))) {
            return false;
        }
        return !((this.matchedSMILES == null) ? (other.matchedSMILES != null) : !this.matchedSMILES.equals(other.matchedSMILES));
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 83 * hash + (this.mol1 != null ? this.mol1.hashCode() : 0);
        hash = 83 * hash + (this.mol2 != null ? this.mol2.hashCode() : 0);
        hash = 83 * hash + (this.indexI != null ? this.indexI.hashCode() : 0);
        hash = 83 * hash + (this.indexJ != null ? this.indexJ.hashCode() : 0);
        hash = 83 * hash + (this.indexStep != null ? this.indexStep.hashCode() : 0);
        hash = 83 * hash + (this.keggMapping ? 1 : 0);
        hash = 83 * hash + (this.rBLASTMapping ? 1 : 0);
        hash = 83 * hash + (this.matchedMol != null ? this.matchedMol.hashCode() : 0);
        hash = 83 * hash + (this.matchedSMILES != null ? this.matchedSMILES.hashCode() : 0);
        return hash;
    }

    /**
     * @param indexStep the indexStep to set
     */
    private synchronized void setIndexStep(Integer indexStep) {
        this.indexStep = indexStep;
    }

    @Override
    public String toString() {
        return "MolMapping{" + "mol1=" + mol1 + ", mol2=" + mol2 + ", indexI=" + indexI + ", indexJ=" + indexJ + ", indexStep=" + indexStep + ", keggMapping=" + keggMapping + ", rBLASTMapping=" + rBLASTMapping + ", matchedMol=" + matchedMol + ", matchedSMILES=" + matchedSMILES + '}';
    }
}
