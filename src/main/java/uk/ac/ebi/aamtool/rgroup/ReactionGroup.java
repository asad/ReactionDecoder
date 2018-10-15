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
package uk.ac.ebi.aamtool.rgroup;

import java.util.Comparator;
import java.util.Objects;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class ReactionGroup implements Comparable<ReactionGroup>, Comparator<ReactionGroup> {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(ReactionGroup.class);

    private final String name;
    private boolean rGroupPresent;

    ReactionGroup(String name) {
        this.name = name;
        this.rGroupPresent = false;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 43 * hash + Objects.hashCode(this.name);
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
        final ReactionGroup other = (ReactionGroup) obj;
        return Objects.equals(this.name, other.name);
    }

    @Override
    public int compareTo(ReactionGroup o) {
        return this.name.compareTo(o.name);
    }

    @Override
    public int compare(ReactionGroup o1, ReactionGroup o2) {
        return o1.name.compareTo(o2.name);
    }

    /**
     * @return the isRGroupPresent
     */
    public boolean isRGroupPresent() {
        return rGroupPresent;
    }

    /**
     * @param isRGroupPresent the isRGroupPresent to set
     */
    public void setRGroupPresent(boolean isRGroupPresent) {
        this.rGroupPresent = isRGroupPresent;
    }
}
