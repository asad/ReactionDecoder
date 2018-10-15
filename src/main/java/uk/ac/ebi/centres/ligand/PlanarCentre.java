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

import static com.google.common.collect.Sets.newHashSet;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Centre;
import uk.ac.ebi.centres.ConnectionProvider;
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
import static uk.ac.ebi.centres.descriptor.Planar.E;
import static uk.ac.ebi.centres.descriptor.Planar.Z;
import static uk.ac.ebi.centres.descriptor.Planar.e;
import static uk.ac.ebi.centres.descriptor.Planar.z;

/**
 * @author John May
 * @param <A>
 */
public class PlanarCentre<A> extends AbstractLigand<A> implements Centre<A> {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(PlanarCentre.class);

    private final AbstractLigand<A> first;
    private final AbstractLigand<A> second;
    private final Set<A> atoms;

    /**
     *
     * @param first
     * @param second
     * @param descriptor
     */
    @SuppressWarnings("unchecked")
    public PlanarCentre(A first, A second,
            MutableDescriptor descriptor) {

        super(descriptor, 0);

        Ligand<A> self = this;

        // create two ligand delegates
        this.first = new NonterminalLigand<>(descriptor, first, second, 0);
        this.second = new NonterminalLigand<>(descriptor, second, first, 0);

        atoms = newHashSet(first, second);

    }

    /**
     *
     * @param provider
     */
    @Override
    public void setProvider(ConnectionProvider<A> provider) {
        super.setProvider(provider);
        first.setProvider(provider);
        second.setProvider(provider);
    }

    @Override
    public List<Ligand<A>> getLigands() {
        List<Ligand<A>> ligands = new ArrayList<>(16);
        ligands.addAll(first.getLigands());
        ligands.addAll(second.getLigands());
        return ligands;
    }

    /**
     *
     * @return
     */
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

    /**
     *
     * @return
     */
    @Override
    public A getParent() {
        throw new UnsupportedOperationException("Can't get parent on a planar centre");
    }

    @Override
    public void setParent(A atom) {
        throw new UnsupportedOperationException("Can't set parent on a planar centre");
    }

    /**
     *
     * @param centres
     * @param rule
     * @param calculator
     * @return
     */
    @Override
    public int perceiveAuxiliary(Collection<Centre<A>> centres, PriorityRule<A> rule, SignCalculator<A> calculator) {
        // System.LOGGER.debug("Auxiliary perception is not currently supported on planar centres");
        return 0;
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
        // can't do this type of perception for planar centres
        return UNKNOWN;
    }

    @Override
    public Descriptor perceive(PriorityRule<A> rule, SignCalculator<A> calculator) {

        List<Ligand<A>> firstLigands = first.getLigands();
        List<Ligand<A>> secondLigands = second.getLigands();

        if (firstLigands.isEmpty() || secondLigands.isEmpty()) {
            return NONE;
        }

        // check for pseudo
        Priority firstPriority = rule.prioritise(firstLigands);
        Priority secondPriority = rule.prioritise(secondLigands);

        if (!firstPriority.isUnique() || !secondPriority.isUnique()) {
            // we don't know whether it is none yet...
            return UNKNOWN;
        }

        int firstSign = calculator.getSign(firstLigands.iterator().next().getAtom(),
                first.getAtom(),
                second.getAtom());
        int secondSign = calculator.getSign(secondLigands.iterator().next().getAtom(),
                second.getAtom(),
                first.getAtom());

        if (firstSign == 0 || secondSign == 0) {
            return UNSPECIFIED;
        }

        boolean pseudo = firstPriority.getType().equals(PSEUDO_ASYMMETRIC)
                && secondPriority.getType().equals(PSEUDO_ASYMMETRIC);

        // also check for psuedo (from prioritise)
        return firstSign == secondSign ? pseudo ? e : E
                : pseudo ? z : Z;

    }

    @Override
    public void dispose() {
        getProvider().dispose();
        setProvider(null);
    }
}
