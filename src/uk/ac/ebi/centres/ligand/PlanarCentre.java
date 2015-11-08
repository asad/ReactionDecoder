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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.ligand;

import com.google.common.collect.Sets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IElement;
import uk.ac.ebi.centres.Centre;
import uk.ac.ebi.centres.ConnectionProvider;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;
import uk.ac.ebi.centres.Priority;
import uk.ac.ebi.centres.PriorityRule;
import uk.ac.ebi.centres.SignCalculator;
import uk.ac.ebi.centres.descriptor.General;
import uk.ac.ebi.centres.descriptor.Planar;

/**
 * @author John May
 */
public class PlanarCentre<A> extends AbstractLigand<A> implements Centre<A> {

    private final AbstractLigand<A> first;
    private final AbstractLigand<A> second;
    private final Set<A> atoms;

    @SuppressWarnings("unchecked")
    public PlanarCentre(A first, A second,
            MutableDescriptor descriptor) {

        super(descriptor, 0);

        Ligand<A> self = this;

        // create two ligand delegates
        this.first = new NonterminalLigand<A>(descriptor, first, second, 0);
        this.second = new NonterminalLigand<A>(descriptor, second, first, 0);

        atoms = Sets.newHashSet(first, second);

    }

    @Override
    public void setProvider(ConnectionProvider<A> provider) {
        super.setProvider(provider);
        first.setProvider(provider);
        second.setProvider(provider);
    }

    @Override
    public List<Ligand<A>> getLigands() {
        List<Ligand<A>> ligands = new ArrayList<Ligand<A>>(16);
        ligands.addAll(first.getLigands());
        ligands.addAll(second.getLigands());
        return ligands;
    }

    @Override
    public A getAtom() {
        // might need a rethink...
        throw new NoSuchMethodError("Centre does not have a single atom");
    }

    @Override
    public String toString() {
        if (first.getAtom() instanceof IAtom) {
            return ((IElement) first.getAtom()).getSymbol() + "" + ((IChemObject) first.getAtom()).getProperty("number") + "="
                    + ((IElement) second.getAtom()).getSymbol() + "" + ((IChemObject) second.getAtom()).getProperty("number");
        }
        return first.getAtom().toString() + "=" + second.getAtom().toString();
    }

    @Override
    public Boolean isParent(Object atom) {
        return atoms.contains(atom);
    }

    @Override
    public Set<A> getAtoms() {
        return atoms;
    }

    @Override
    public A getParent() {
        throw new UnsupportedOperationException("Can't get parent on a planar centre");
    }

    @Override
    public void setParent(A atom) {
        throw new UnsupportedOperationException("Can't set parent on a planar centre");
    }

    @Override
    public int perceiveAuxiliary(Collection<Centre<A>> centres, PriorityRule<A> rule, SignCalculator<A> calculator) {
        // System.err.println("Auxiliary perception is not currently supported on planar centres");
        return 0;
    }

    @Override
    public Descriptor perceive(List<Ligand<A>> proximal, PriorityRule<A> rule, SignCalculator<A> calculator) {
        // can't do this type of perception for planar centres
        return General.UNKNOWN;
    }

    @Override
    public Descriptor perceive(PriorityRule<A> rule, SignCalculator<A> calculator) {

        List<Ligand<A>> firstLigands = first.getLigands();
        List<Ligand<A>> secondLigands = second.getLigands();

        if (firstLigands.isEmpty() || secondLigands.isEmpty()) {
            return General.NONE;
        }

        // check for pseudo
        Priority firstPriority = rule.prioritise(firstLigands);
        Priority secondPriority = rule.prioritise(secondLigands);

        if (!firstPriority.isUnique() || !secondPriority.isUnique()) {
            // we don't know whether it is none yet...
            return General.UNKNOWN;
        }

        int firstSign = calculator.getSign(firstLigands.iterator().next().getAtom(),
                first.getAtom(),
                second.getAtom());
        int secondSign = calculator.getSign(secondLigands.iterator().next().getAtom(),
                second.getAtom(),
                first.getAtom());

        if (firstSign == 0 || secondSign == 0) {
            return General.UNSPECIFIED;
        }

        boolean pseudo = firstPriority.getType().equals(Descriptor.Type.PSEUDO_ASYMMETRIC)
                && secondPriority.getType().equals(Descriptor.Type.PSEUDO_ASYMMETRIC);

        // also check for psuedo (from prioritise)
        return firstSign == secondSign ? pseudo ? Planar.e : Planar.E
                : pseudo ? Planar.z : Planar.Z;

    }

    @Override
    public void dispose() {
        getProvider().dispose();
        setProvider(null);
    }
    private static final Logger LOG = Logger.getLogger(PlanarCentre.class.getName());
}
