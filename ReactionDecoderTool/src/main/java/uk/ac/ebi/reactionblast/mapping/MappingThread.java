/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.mapping;

import java.util.concurrent.Callable;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class MappingThread implements Callable<Reactor> {
    private static final ILoggingTool logger = createLoggingTool(MappingThread.class);
    private static final Logger LOG = getLogger(MappingThread.class.getName());

    private final IReaction cleanedReaction;
    private final IMappingAlgorithm algorithm;
    private final boolean removeHydrogen;

    /**
     *
     * @param reaction to be mapped (only balanced reactions are mapped)
     * @param removeHydrogen true (map without hydrogen, for faster mapping)
     * else false for complete with Hydrogen
     * @param algorithm
     *
     * @return Mapped Object
     */
    MappingThread(String message, IReaction cleanedReaction, IMappingAlgorithm algorithm, boolean removeHydrogen) {
        this.cleanedReaction = cleanedReaction;
        this.algorithm = algorithm;
        this.removeHydrogen = removeHydrogen;
        logger.info("|++++++++++++++++++++++++++++|");
        logger.info("|Atom Atom Mapping Tool Initialized for " + message);
    }

    @Override
    public Reactor call() throws Exception {
        try {
            Reactor reactor;
            reactor = new Reactor(cleanedReaction, removeHydrogen, algorithm);
            logger.info("|Done " + reactor.getAlgorithm() + " |");
            return reactor;
        } catch (Exception ex) {
            throw ex;
        }
    }
}
