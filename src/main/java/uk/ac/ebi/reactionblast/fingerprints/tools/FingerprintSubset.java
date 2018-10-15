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
package uk.ac.ebi.reactionblast.fingerprints.tools;

import java.io.Serializable;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FingerprintSubset implements Serializable {

    private static final long serialVersionUID = 4342623464361L;

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(FingerprintSubset.class);

    /**
     * Determine if this set is an (improper) subset of another.
     *
     * @param source the set we are testing for.
     * @param destination the set we are testing against.
     * @return source is a subset of destination, yes then return true else
     * false
     * @throws CDKException
     */
    public static boolean isSubset(BitSet source, BitSet destination) throws CDKException {
        boolean flag = false;
        if (source.cardinality() <= destination.cardinality()) {
            not_null(source);

            /* make a copy of the source set */
            BitSet copy_other = (BitSet) source.clone();

            /* and or in */
            copy_other.and(destination);

            /* if it hasn't changed, we were a subset */
            flag = copy_other.equals(source);
        }

        return flag;
    }

    /*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    /**
     * Determine if this set is an (improper) superset of another.
     *
     * @param source the set we are testing for.
     * @param destination the set we are testing against.
     * @return source is a superset of destination, yes then return true else
     * false
     * @throws CDKException
     */
    public static boolean isSuperSet(BitSet source, BitSet destination) throws CDKException {
        boolean flag = false;
        if (source.cardinality() >= destination.cardinality()) {

            not_null(source);

            /* make a copy of the source set */
            BitSet copy_other = (BitSet) source.clone();

            /* and or in */
            copy_other.and(destination);

            /* if it hasn't changed, we were a subset */
            flag = copy_other.equals(destination);
//            flag = copy_other.equals(destination);
        }
        return flag;
    }

    /**
     * Helper function to test for a null object and throw an exception if one
     * is found.
     *
     * @param obj the object we are testing.
     * @throws CDKException
     */
    protected static void not_null(Object obj) throws CDKException {
        if (obj == null) {
            throw new CDKException("Null object used in set operation");
        }
    }

    private FingerprintSubset() {
    }
}
