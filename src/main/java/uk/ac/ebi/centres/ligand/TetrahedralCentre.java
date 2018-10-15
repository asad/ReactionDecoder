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
package uk.ac.ebi.centres.ligand;

import java.util.Collection;
import static java.util.Collections.singleton;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Centre;
import uk.ac.ebi.centres.Descriptor;
import static uk.ac.ebi.centres.Descriptor.Type.PSEUDO_ASYMMETRIC;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;
import uk.ac.ebi.centres.Priority;
import uk.ac.ebi.centres.PriorityRule;
import uk.ac.ebi.centres.SignCalculator;
import static uk.ac.ebi.centres.descriptor.General.NONE;
import static uk.ac.ebi.centres.descriptor.General.UNKNOWN;
import static uk.ac.ebi.centres.descriptor.General.UNSPECIFIED;
import static uk.ac.ebi.centres.descriptor.Tetrahedral.R;
import static uk.ac.ebi.centres.descriptor.Tetrahedral.S;
import static uk.ac.ebi.centres.descriptor.Tetrahedral.r;
import static uk.ac.ebi.centres.descriptor.Tetrahedral.s;

/**
 * @author John May
 * @param <A>
 */
public class TetrahedralCentre<A>
        extends AbstractLigand<A>
        implements Centre<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(TetrahedralCentre.class);

    private final A atom;
    private A parent;

    /**
     *
     * @param descriptor
     * @param atom
     */
    public TetrahedralCentre(MutableDescriptor descriptor,
            A atom) {
        super(descriptor, 0);
        this.atom = atom;
        this.parent = atom;
    }

    /**
     *
     * @return
     */
    @Override
    public A getAtom() {
        return atom;
    }

    @Override
    public void setParent(A atom) {
        // don't have a parent here
        this.parent = atom;
    }

    /**
     *
     * @return
     */
    @Override
    public A getParent() {
        return this.parent;
    }

    @Override
    public Set<A> getAtoms() {
        return singleton(atom);
    }

    /**
     *
     * @param centres
     * @param rule
     * @param calculator
     * @return
     */
    @Override
    public int perceiveAuxiliary(Collection<Centre<A>> centres,
            PriorityRule<A> rule,
            SignCalculator<A> calculator) {

        Map<Ligand<A>, Descriptor> auxiliary = new HashMap<>(centres.size());

        getProvider().build();

        for (Centre<A> centre : centres) {

            // don't do aux perception on self
            if (centre == this) {
                continue;
            }

            // can only reroot on single atom centres
            if (centre.getAtoms().size() == 1) {

                for (Ligand<A> ligand : getProvider().getLigands(centre.getAtom())) {

                    getProvider().reroot(ligand);

                    Descriptor descriptor = centre.perceive(getProvider().getLigands(ligand),
                            rule,
                            calculator);

                    if (descriptor != UNKNOWN) {
                        auxiliary.put(ligand, descriptor);
                    }

                }
            }
        }

        // transfer auxiliary descriptors to their respective ligands
        for (Map.Entry<Ligand<A>, Descriptor> entry : auxiliary.entrySet()) {
            entry.getKey().setAuxiliary(entry.getValue());
        }

        // reroot on this
        getProvider().reroot(this);

        return auxiliary.size();

    }

    /**
     *
     * @param proximal
     * @param rule
     * @param calculator
     * @return
     */
    @Override
    public Descriptor perceive(List<Ligand<A>> proximal, PriorityRule<A> rule, SignCalculator<A> calculator) {

        if (proximal.size() < 3) {
            return NONE;
        }

        Priority priority = rule.prioritise(proximal);

        if (priority.isUnique()) {

            if (proximal.size() < 4) {
                proximal.add(this);
            }

            int sign = calculator.getSign(proximal.get(0),
                    proximal.get(1),
                    proximal.get(2),
                    proximal.get(3));

            boolean pseudo = priority.getType().equals(PSEUDO_ASYMMETRIC);

            return sign > 0 ? pseudo ? s
                    : S
                    : sign < 0
                            ? pseudo ? r
                                    : R
                            : UNSPECIFIED;

        }

        return UNKNOWN;
    }

    @Override
    public Descriptor perceive(PriorityRule<A> rule, SignCalculator<A> calculator) {

        return perceive(getLigands(), rule, calculator);

    }

    @Override
    public Boolean isParent(A atom) {
        return parent.equals(atom);
    }

    @Override
    public void dispose() {
        getProvider().dispose();
        setProvider(null);
    }
}
