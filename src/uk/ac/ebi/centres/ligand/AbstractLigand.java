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

import com.google.common.collect.Sets;
import org.openscience.cdk.interfaces.IAtom;
import uk.ac.ebi.centres.ConnectionProvider;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;
import uk.ac.ebi.centres.descriptor.General;
import uk.ac.ebi.centres.graph.Arc;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * @author John May
 */
public abstract class AbstractLigand<A> implements Ligand<A> {

    private Descriptor auxiliary = General.UNKNOWN;
    private ConnectionProvider<A> provider;
    private final Set<A> visited;
    private final MutableDescriptor descriptor;
    private final int distance;
    private boolean duplicate;
    private List<Ligand<A>> ligands;
    private Descriptor descriptorCache;

    public AbstractLigand(ConnectionProvider<A> provider,
            Set<A> visited,
            MutableDescriptor descriptor,
            int distance) {

        this.provider = provider;
        this.descriptor = descriptor;
        this.distance = distance;

        // optimise size for a load factor of 0.75
        this.visited = Sets.newHashSet(visited);

    }

    public AbstractLigand(Set<A> visited,
            MutableDescriptor descriptor,
            int distance) {

        this.descriptor = descriptor;
        this.distance = distance;

        // optimise size for a load factor of 0.75
        this.visited = Sets.newHashSet(visited);

    }

    public AbstractLigand(MutableDescriptor descriptor,
            int distance) {

        this.descriptor = descriptor;
        this.distance = distance;

        this.visited = Collections.EMPTY_SET;

    }

    public boolean isDuplicate() {
        return duplicate;
    }

    public void setDuplicate(boolean duplicate) {
        this.duplicate = duplicate;
    }

    public ConnectionProvider<A> getProvider() {
        return provider;
    }

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
            if (descriptor == General.NONE) // cache access to NONE descriptors
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
            return ((IAtom) atom).getSymbol() + "" + ((IAtom) atom).getProperty("number");
        }
        return "Non CDK Atom";
    }

    @Override
    public List<Arc<A>> getArcs() {
        return provider.getArcs(this);
    }

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

    @Override
    public boolean isBranching() {
        return Boolean.FALSE;
    }

    @Override
    public boolean isTerminal() {
        return Boolean.FALSE;
    }
}
