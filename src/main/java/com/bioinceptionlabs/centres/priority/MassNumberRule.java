/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package com.bioinceptionlabs.centres.priority;

import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.centres.Ligand;
import static com.bioinceptionlabs.centres.PriorityRule.Type.CONSTITUTIONAL;
import com.bioinceptionlabs.centres.priority.access.MassNumberAccessor;

/**
 * An abstract class for constitutional priority based on mass number. A mass
 * number accessor
 * ({@link com.bioinceptionlabs.centres.priority.access.MassNumberAccessor}) can be
 * provided to allow the comparator to work on a custom atom type.
 *
 * @author John May
 * @param <A>
 */
public class MassNumberRule<A>
        extends AbstractPriorityRule<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MassNumberRule.class);

    /**
     * Accessor used to get the atomic number from an atom.
     */
    private final MassNumberAccessor<A> accessor;

    /**
     * Constructs an mass number comparator that uses the provided accessor to
     * fetch the mass number for a given atom.
     *
     * @param accessor an accessor for the atom's mass number
     */
    public MassNumberRule(MassNumberAccessor<A> accessor) {
        super(CONSTITUTIONAL);
        this.accessor = accessor;
    }

    /**
     * Compares the ligands by their atoms mass numbers.
     *
     * @inheritDoc
     */
    @Override
    public int compare(Ligand<A> o1, Ligand<A> o2) {
        if (accessor == null || o1.getAtom() == null) {
            LOGGER.debug(accessor + " 1 NULL");
        }
        if (accessor == null || o2.getAtom() == null) {
            LOGGER.debug(accessor + " 2 NULL");
        }
        IAtom a = (IAtom) o1.getAtom();
        IAtom b = (IAtom) o2.getAtom();

        return getMassNumber(a) - getMassNumber(b);
    }

    public int getMassNumber(IAtom o) {

        if (o != null && !(o instanceof IPseudoAtom) && o.getMassNumber() == null) {
            try {
                int massNumber = Isotopes.getInstance().getMajorIsotope(o.getAtomicNumber()).getMassNumber();
                return massNumber;
            } catch (Exception e) {
                return 11;
            }
        } else if (o instanceof IPseudoAtom) {
            return 11;//less than carbon for 'R'
        }
        return 0;
    }
}
