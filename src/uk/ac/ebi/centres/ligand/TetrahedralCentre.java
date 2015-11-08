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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import uk.ac.ebi.centres.Centre;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;
import uk.ac.ebi.centres.Priority;
import uk.ac.ebi.centres.PriorityRule;
import uk.ac.ebi.centres.SignCalculator;
import uk.ac.ebi.centres.descriptor.General;
import uk.ac.ebi.centres.descriptor.Tetrahedral;

/**
 * @author John May
 */
public class TetrahedralCentre<A>
        extends AbstractLigand<A>
        implements Centre<A> {
    private static final Logger LOG = Logger.getLogger(TetrahedralCentre.class.getName());

    private final A atom;
    private A parent;

    public TetrahedralCentre(MutableDescriptor descriptor,
            A atom) {
        super(descriptor, 0);
        this.atom = atom;
        this.parent = atom;
    }

    @Override
    public A getAtom() {
        return atom;
    }

    @Override
    public void setParent(A atom) {
        // don't have a parent here
        this.parent = atom;
    }

    @Override
    public A getParent() {
        return this.parent;
    }

    @Override
    public Set<A> getAtoms() {
        return Collections.singleton(atom);
    }

    @Override
    public int perceiveAuxiliary(Collection<Centre<A>> centres,
            PriorityRule<A> rule,
            SignCalculator<A> calculator) {

        Map<Ligand<A>, Descriptor> auxiliary = new HashMap<Ligand<A>, Descriptor>(centres.size());

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

                    if (descriptor != General.UNKNOWN) {
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

    @Override
    public Descriptor perceive(List<Ligand<A>> proximal, PriorityRule<A> rule, SignCalculator<A> calculator) {

        if (proximal.size() < 3) {
            return General.NONE;
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

            boolean pseudo = priority.getType().equals(Descriptor.Type.PSEUDO_ASYMMETRIC);

            return sign > 0 ? pseudo ? Tetrahedral.s
                    : Tetrahedral.S
                    : sign < 0
                    ? pseudo ? Tetrahedral.r
                    : Tetrahedral.R
                    : General.UNSPECIFIED;


        }

        return General.UNKNOWN;
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
