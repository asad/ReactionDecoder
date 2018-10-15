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
package uk.ac.ebi.centres;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import static java.util.concurrent.Executors.newSingleThreadExecutor;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.centres.descriptor.General.NONE;
import static uk.ac.ebi.centres.descriptor.General.UNKNOWN;

/**
 * @author John May
 * @param <A>
 */
public class DefaultPerceptor<A> implements Perceptor<A> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(DefaultPerceptor.class);

    private final CentrePerceptor<A> mainPerceptor;
    private final CentrePerceptor<A> auxPerceptor;
    private ExecutorService executor = newSingleThreadExecutor();
    private long timeout = 250;

    /**
     *
     * @param rule
     * @param auxRule
     * @param calculator
     */
    public DefaultPerceptor(final PriorityRule<A> rule,
            final PriorityRule<A> auxRule,
            final SignCalculator<A> calculator) {

        // create the main and aux perceptors
        this.mainPerceptor = new CentrePerceptor<A>(rule) {
            @Override
            public Descriptor perceive(Centre<A> centre, Collection<Centre<A>> centres) {
                return centre.perceive(rule, calculator);
            }
        };
        this.auxPerceptor = new CentrePerceptor<A>(auxRule) {
            @Override
            public Descriptor perceive(Centre<A> centre, Collection<Centre<A>> centres) {
                // only attempt re-perception if there were auxiliary labels defined
                return centre.perceiveAuxiliary(centres, rule, calculator) != 0
                        ? centre.perceive(auxRule, calculator)
                        : UNKNOWN;
            }
        };
    }

    private List<Centre<A>> _perceive(Collection<Centre<A>> unperceived,
            CentrePerceptor<A> perceptor) {

        List<Centre<A>> perceived = new ArrayList<>();
        Map<Centre<A>, Descriptor> map = new LinkedHashMap<>();

        do {

            map.clear();

            unperceived.forEach((centre) -> {
                Descriptor descriptor = perceptor.perceive(centre, unperceived);
                if (descriptor != UNKNOWN) {
                    map.put(centre, descriptor);
                }
            });

            // transfer descriptors
            map.entrySet().stream().map((entry) -> {
                unperceived.remove(entry.getKey());
                return entry;
            }).map((entry) -> {
                perceived.add(entry.getKey());
                return entry;
            }).map((entry) -> {
                entry.getKey().dispose();
                return entry;
            }).forEachOrdered((entry) -> {
                entry.getKey().setDescriptor(entry.getValue());
            });

        } while (!map.isEmpty());

        return perceived;

    }

    /**
     *
     * @param provider
     * @param manager
     */
    @Override
    public void perceive(final CentreProvider<A> provider, final DescriptorManager<A> manager) {

        // timeout fo the centre provider incase we have a huge molecule and the spanning tree can't
        // be constructed
        Collection<Centre<A>> unperceived = provider.getCentres(manager);

        if (unperceived.isEmpty()) {
            return;
        }

        // could switch to only use this on large molecule
        List<Centre<A>> perceived = _perceive(unperceived, mainPerceptor);

        // no centres perceived, perform auxiliary perception
        if (!unperceived.isEmpty() && perceived.isEmpty()) {
            perceived.addAll(_perceive(unperceived, auxPerceptor));
        }

        // set all unperceived centres to 'none'
        for (Centre<A> centre : unperceived) {
            centre.setDescriptor(NONE);
            centre.dispose();
        }

        unperceived.clear();
        unperceived = null;
        manager.clear();

    }

    /**
     * Shutdown the internal executor
     */
    @Override
    public void shutdown() {
        executor.shutdownNow();
    }

    abstract class CentrePerceptor<A> {

        private final PriorityRule<A> rule;

        protected CentrePerceptor(PriorityRule<A> rule) {
            this.rule = rule;
        }

        public abstract Descriptor perceive(Centre<A> centre, Collection<Centre<A>> centres);
    }
}
