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
package uk.ac.ebi.centres.cdk;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Descriptor;
import uk.ac.ebi.centres.DescriptorManager;
import uk.ac.ebi.centres.MutableDescriptor;
import static uk.ac.ebi.centres.descriptor.General.UNKNOWN;

/**
 * @author John May
 */
public class CDKManager implements DescriptorManager<IAtom> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(DescriptorManager.class);

    private final IAtomContainer container;
    private final Map<IChemObject, MutableDescriptor> map = new HashMap<>();

    /**
     *
     * @param container
     */
    public CDKManager(IAtomContainer container) {
        this.container = container;
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public MutableDescriptor getDescriptor(IAtom atom) {
        return _getDescriptor(atom);
    }

    /**
     *
     * @param first
     * @param second
     * @return
     */
    @Override
    public MutableDescriptor getDescriptor(IAtom first, IAtom second) {
        return _getDescriptor(container.getBond(first, second));
    }

    private MutableDescriptor _getDescriptor(IChemObject chemObject) {
        MutableDescriptor mutator = map.get(chemObject);
        if (mutator == null) {
            mutator = new ProxyMutator(chemObject);
            map.put(chemObject, mutator);
        }
        return mutator;
    }

    /**
     *
     */
    @Override
    public void clear() {
        map.clear();
    }

    class ProxyMutator extends MutableDescriptor {

        private final IChemObject chemObject;

        ProxyMutator(IChemObject chemObject) {
            this.chemObject = chemObject;
            chemObject.setProperty("descriptor", UNKNOWN);
        }

        @Override
        public synchronized Descriptor get() {
            return (Descriptor) chemObject.getProperty("descriptor");
        }

        @Override
        public synchronized void set(Descriptor descriptor) {
            chemObject.setProperty("descriptor", descriptor);
        }
    }
}
