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

import static com.google.common.collect.Sets.newHashSet;
import static java.lang.Boolean.FALSE;
import static java.util.Collections.EMPTY_SET;
import java.util.List;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IElement;
import uk.ac.ebi.centres.ConnectionProvider;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;
import static uk.ac.ebi.centres.descriptor.General.NONE;
import static uk.ac.ebi.centres.descriptor.General.UNKNOWN;
import uk.ac.ebi.centres.graph.Arc;

/**
 * @author John May
 * @param <A>
 */
public abstract class AbstractLigand<A> implements Ligand<A> {

    private Descriptor auxiliary = UNKNOWN;
    private ConnectionProvider<A> provider;
    private final Set<A> visited;
    private final MutableDescriptor descriptor;
    private final int distance;
    private boolean duplicate;
    private List<Ligand<A>> ligands;
    private Descriptor descriptorCache;

    /**
     *
     * @param provider
     * @param visited
     * @param descriptor
     * @param distance
     */
    public AbstractLigand(ConnectionProvider<A> provider,
            Set<A> visited,
            MutableDescriptor descriptor,
            int distance) {

        this.provider = provider;
        this.descriptor = descriptor;
        this.distance = distance;

        // optimise size for a load factor of 0.75
        this.visited = newHashSet(visited);

    }

    /**
     *
     * @param visited
     * @param descriptor
     * @param distance
     */
    public AbstractLigand(Set<A> visited,
            MutableDescriptor descriptor,
            int distance) {

        this.descriptor = descriptor;
        this.distance = distance;

        // optimise size for a load factor of 0.75
        this.visited = newHashSet(visited);

    }

    /**
     *
     * @param descriptor
     * @param distance
     */
    public AbstractLigand(MutableDescriptor descriptor,
            int distance) {

        this.descriptor = descriptor;
        this.distance = distance;

        this.visited = EMPTY_SET;

    }

    /**
     *
     * @return
     */
    public boolean isDuplicate() {
        return duplicate;
    }

    /**
     *
     * @param duplicate
     */
    public void setDuplicate(boolean duplicate) {
        this.duplicate = duplicate;
    }

    /**
     *
     * @return
     */
    public ConnectionProvider<A> getProvider() {
        return provider;
    }

    /**
     *
     * @param provider
     */
    public void setProvider(ConnectionProvider<A> provider) {
        this.provider = provider;
    }

    @Override
    public Boolean isVisited(A atom) {
        return visited.contains(atom);
    }

    @Override
    public Set<A> getVisited() {
        return visited;
    }

    @Override
    public void setDescriptor(Descriptor descriptor) {
        this.descriptor.set(descriptor);
    }

    @Override
    public Descriptor getDescriptor() {
        if (descriptorCache == null) {
            Descriptor descriptor = this.descriptor.get();
            if (descriptor == NONE) // cache access to NONE descriptors
            {
                descriptorCache = descriptor;
            }
            return descriptor;
        }
        return descriptorCache;
    }

    /**
     * @inheritDoc
     */
    @Override
    public List<Ligand<A>> getLigands() {
        if (ligands == null) {
            ligands = provider.getLigands(this);
        }
        return ligands;
    }

    public void reset() {
        ligands = null;
    }

    @Override
    public String toString() {
        A atom = getAtom();
        if (atom instanceof IAtom) {
            return ((IElement) atom).getSymbol() + "" + ((IChemObject) atom).getProperty("number");
        }
        return "Non CDK Atom";
    }

    /**
     *
     * @return
     */
    @Override
    public List<Arc<A>> getArcs() {
        return provider.getArcs(this);
    }

    /**
     *
     * @return
     */
    @Override
    public Arc<A> getParentArc() {
        return provider.getParentArc(this);
    }

    @Override
    public int getDistanceFromRoot() {
        return distance;
    }

    /**
     * @inheritDoc
     */
    @Override
    public Descriptor getAuxiliary() {
        return auxiliary;
    }

    /**
     * @inheritDoc
     */
    @Override
    public void setAuxiliary(Descriptor descriptor) {
        this.auxiliary = descriptor;
    }

    @Override
    public int getDepth() {
        Arc<A> arc = getParentArc();
        return arc == null ? 0 : arc.getDepth();
    }

    /**
     *
     * @return
     */
    @Override
    public boolean isBranching() {
        return FALSE;
    }

    /**
     *
     * @return
     */
    @Override
    public boolean isTerminal() {
        return FALSE;
    }
}
