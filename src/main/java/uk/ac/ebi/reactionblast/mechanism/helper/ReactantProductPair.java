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
public class ReactantProductPair implements Serializable,
        Comparable<ReactantProductPair>, Comparator<ReactantProductPair> {

    private static final long serialVersionUID = 19876565735478L;
    private final String query;
    private final String target;

    /**
     *
     * @param query
     * @param target
     */
    public ReactantProductPair(String query, String target) {
        this.query = query;
        this.target = target;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 97 * hash + Objects.hashCode(this.query);
        hash = 97 * hash + Objects.hashCode(this.target);
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
        final ReactantProductPair other = (ReactantProductPair) obj;
        if (!Objects.equals(this.query, other.query)) {
            return false;
        }
        return Objects.equals(this.target, other.target);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("R:").append(this.getQuery());
        sb.append(", P:").append(this.getTarget());
        return sb.toString();
    }

    /**
     * @return the query
     */
    public String getQuery() {
        return query;
    }

    /**
     * @return the target
     */
    public String getTarget() {
        return target;
    }

    @Override
    public int compareTo(ReactantProductPair o) {
        String local = this.query + this.target;
        String object = o.getQuery() + o.getTarget();
        return local.compareTo(object);
    }

    @Override
    public int compare(ReactantProductPair o1, ReactantProductPair o2) {
        return o1.compareTo(o2);
    }
}
