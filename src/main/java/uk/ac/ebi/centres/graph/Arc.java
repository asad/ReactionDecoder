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
package uk.ac.ebi.centres.graph;

import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.Ligand;
import uk.ac.ebi.centres.MutableDescriptor;

/**
 * @author John May
 * @param <A>
 */
public class Arc<A> {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Arc.class);

    private Ligand<A> tail;
    private Ligand<A> head;
    private MutableDescriptor descriptor;
    private int depth = 0;

    /**
     *
     * @param tail
     * @param head
     * @param descriptor
     */
    public Arc(Ligand<A> tail, Ligand<A> head,
            MutableDescriptor descriptor) {
        this.tail = tail;
        this.head = head;
        this.descriptor = descriptor;
    }

    /**
     * @param tail
     * @param head
     * @param descriptor
     * @param depth 1 = tail is closed then head, -1 = head is closer then head, 0 = same plane. -1 = wedge bond when
     * tail is root atom
     */
    public Arc(Ligand<A> tail,
            Ligand<A> head,
            MutableDescriptor descriptor,
            int depth) {
        this.tail = tail;
        this.head = head;
        this.descriptor = descriptor;
        this.depth = depth;
    }

    /**
     *
     * @return
     */
    public int getDepth() {
        return depth;
    }

    /**
     *
     * @return
     */
    public Descriptor getDescriptor() {
        return descriptor.get();
    }

    /**
     *
     * @return
     */
    public Ligand<A> getHead() {
        return this.head;
    }

    /**
     *
     * @return
     */
    public Ligand<A> getTail() {
        return this.tail;
    }

    /**
     *
     */
    public void transpose() {
        Ligand<A> tmp = tail;
        tail = head;
        head = tmp;
        depth *= -1; // invert the sign
        head.setParent(tail.getAtom());
        head.reset(); // need to reset any caches
        tail.reset();
    }

    @Override
    public String toString() {
        return tail + " -> " + head;
    }
}
