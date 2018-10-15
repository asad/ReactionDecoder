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
package uk.ac.ebi.centres.priority;

import static java.lang.System.out;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Ligand;
import static uk.ac.ebi.centres.PriorityRule.Type.CONSTITUTIONAL;
import uk.ac.ebi.centres.priority.access.AtomicNumberAccessor;

/**
 * An abstract class for constitutional priority based on atomic number. An
 * atomic number accessor ({@link AtomicNumberAccessor}) can be provided to
 * allow the comparator to work on a custom atom type.
 *
 * @author John May
 * @param <A>
 */
public class AtomicNumberRule<A>
        extends AbstractPriorityRule<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(AtomicNumberRule.class);

    /**
     * Accessor used to get the atomic number from an atom.
     */
    private final AtomicNumberAccessor<A> accessor;

    /**
     * Constructs an atomic number comparator that uses the provided accessor to
     * fetch the atomic number for a given atom.
     *
     * @param accessor an accessor for the atom's atomic number
     */
    public AtomicNumberRule(AtomicNumberAccessor<A> accessor) {
        super(CONSTITUTIONAL);
        this.accessor = accessor;
    }

    /**
     * Compares the ligands by their atoms atomic numbers.
     *
     * @inheritDoc
     */
    @Override
    public int compare(Ligand<A> o1, Ligand<A> o2) {
        if (accessor == null || o1.getAtom() == null) {
            out.println(accessor + " 1 NULL");
        }
        if (accessor == null || o2.getAtom() == null) {
            out.println(accessor + " 2 NULL");
        }
        return accessor.getAtomicNumber(o1.getAtom()) - accessor.getAtomicNumber(o2.getAtom());
    }
}
