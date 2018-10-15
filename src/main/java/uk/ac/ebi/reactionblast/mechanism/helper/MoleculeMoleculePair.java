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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Objects;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MoleculeMoleculePair implements Serializable, Comparable<MoleculeMoleculePair>, Comparator<MoleculeMoleculePair> {

    private static final long serialVersionUID = 107097779868968L;
    private final ReactantProductPair name;
    private final ReactantProductPair smarts;
    private final ReactantProductPair signature;
    private final String smirks;
    private ReactantProductPair smarts1;
    private ReactantProductPair signature1;
    private String smirks1;
    private ReactantProductPair smarts2;
    private ReactantProductPair signature2;
    private String smirks2;
    private ReactantProductPair smarts3;
    private ReactantProductPair signature3;
    private String smirks3;

    /**
     *
     * @param name
     * @param smarts
     * @param signature
     * @param smirks
     */
    public MoleculeMoleculePair(
            ReactantProductPair name,
            ReactantProductPair smarts,
            ReactantProductPair signature,
            String smirks) {
        this.name = name;
        this.smarts = smarts;
        this.signature = signature;
        this.smirks = smirks;
    }

    /**
     * @return the name
     */
    public ReactantProductPair getName() {
        return name;
    }

    /**
     * @return the smarts
     */
    public ReactantProductPair getSmarts() {
        return smarts;
    }

    /**
     * @return the signature
     */
    public ReactantProductPair getSignature() {
        return signature;
    }

    /**
     * @return the smirks
     */
    public String getSmirks() {
        return smirks;
    }

    /**
     * @return the smarts at level 1
     */
    public ReactantProductPair getSmarts1() {
        return smarts1;
    }

    /**
     * @param smarts1 the smarts at level 1 to set
     */
    public void setSmarts1(ReactantProductPair smarts1) {
        this.smarts1 = smarts1;
    }

    /**
     * @return the signature at level 1
     */
    public ReactantProductPair getSignature1() {
        return signature1;
    }

    /**
     * @param signature1 the signature1 to set
     */
    public void setSignature1(ReactantProductPair signature1) {
        this.signature1 = signature1;
    }

    /**
     * @return the smirks at level 1
     */
    public String getSmirks1() {
        return smirks1;
    }

    /**
     * @param smirks1 the smirks at level 1 to set
     */
    public void setSmirks1(String smirks1) {
        this.smirks1 = smirks1;
    }

    /**
     * @return the smarts at level 2
     */
    public ReactantProductPair getSmarts2() {
        return smarts2;
    }

    /**
     * @param smarts2 the smarts at level 2 to set
     */
    public void setSmarts2(ReactantProductPair smarts2) {
        this.smarts2 = smarts2;
    }

    /**
     * @return the signature at level 2
     */
    public ReactantProductPair getSignature2() {
        return signature2;
    }

    /**
     * @param signature2 the signature at level 2 to set
     */
    public void setSignature2(ReactantProductPair signature2) {
        this.signature2 = signature2;
    }

    /**
     * @return the smirks at level 2
     */
    public String getSmirks2() {
        return smirks2;
    }

    /**
     * @param smirks2 the smirks at level 2 to set
     */
    public void setSmirks2(String smirks2) {
        this.smirks2 = smirks2;
    }

    /**
     * @return the smarts at level 3
     */
    public ReactantProductPair getSmarts3() {
        return smarts3;
    }

    /**
     * @param smarts3 the smarts3 to set
     */
    public void setSmarts3(ReactantProductPair smarts3) {
        this.smarts3 = smarts3;
    }

    /**
     * @return the signature at level 3
     */
    public ReactantProductPair getSignature3() {
        return signature3;
    }

    /**
     * @param signature3 the signature at level 3 to set
     */
    public void setSignature3(ReactantProductPair signature3) {
        this.signature3 = signature3;
    }

    /**
     * @return the smirks at level 3
     */
    public String getSmirks3() {
        return smirks3;
    }

    /**
     * @param smirks3 the smirks at level 3 to set
     */
    public void setSmirks3(String smirks3) {
        this.smirks3 = smirks3;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Name:").append(this.name);
        sb.append("signature:").append(this.signature);
        sb.append("smarts:").append(this.smarts);
        sb.append("smirks:").append(this.smirks);
        return sb.toString();
    }

    @Override
    public int compareTo(MoleculeMoleculePair o) {
        String local = this.name + this.smirks;
        String other = o.getName() + o.getSmirks();
        return local.compareTo(other);
    }

    @Override
    public int compare(MoleculeMoleculePair o1, MoleculeMoleculePair o2) {
        return o1.compareTo(o2);
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 37 * hash + Objects.hashCode(this.name);
        hash = 37 * hash + Objects.hashCode(this.smarts);
        hash = 37 * hash + Objects.hashCode(this.signature);
        hash = 37 * hash + Objects.hashCode(this.smirks);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MoleculeMoleculePair other = (MoleculeMoleculePair) obj;
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        if (!Objects.equals(this.smarts, other.smarts)) {
            return false;
        }
        if (!Objects.equals(this.signature, other.signature)) {
            return false;
        }
        return Objects.equals(this.smirks, other.smirks);
    }
}
