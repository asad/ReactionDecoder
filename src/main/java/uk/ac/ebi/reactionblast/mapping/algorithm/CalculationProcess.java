/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.algorithm;
//~--- non-JDK imports --------------------------------------------------------

import java.io.Serializable;
import static java.lang.System.getProperty;
import static java.util.Collections.synchronizedSortedMap;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.mapping.algorithm.GameTheoryFactory.make;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;
import uk.ac.ebi.reactionblast.mapping.interfaces.IGameTheory;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MAX;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIN;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.MIXTURE;
import static uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm.RINGS;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CalculationProcess extends IsomeraseHandler implements Serializable {

    static final String NEW_LINE = getProperty("line.separator");
    private final static ILoggingTool LOGGER
            = createLoggingTool(CalculationProcess.class);
    private static final long serialVersionUID = 0x4a0bba049L;
    private final boolean removeHydrogen;
    private int delta = 0;
    private MoleculeMoleculeMapping reactionBlastMolMapping;
    private final IMappingAlgorithm algorithm;

    /**
     *
     * @param removeHydrogen
     * @param reaction
     * @param algorithm
     * @throws org.openscience.cdk.exception.Intractable
     */
    public CalculationProcess(
            boolean removeHydrogen,
            IReaction reaction,
            IMappingAlgorithm algorithm) throws Intractable {

        /*
         * This case handles rings cases where 6 membered ring reduces to 5 membered rings Example KEGG reaction R01432
         * of Isomerase class
         */
        super(reaction);

//        System.out.println("I am CalculationProcess");
        this.removeHydrogen = removeHydrogen;
        LOGGER.debug(NEW_LINE + "|++++++++++++++++++++++++++++|");
        LOGGER.debug("Performing Atom-Atom Mapping ....... " + reaction.getID() + " .......");
        LOGGER.debug(NEW_LINE + "|++++++++++++++++++++++++++++|");
        this.algorithm = algorithm;
        run();
    }

    private synchronized void run() {
        switch (algorithm) {
            case MIN:
                LOGGER.debug("Processing Reaction for Local Minimum: ");
                delta = (int) calRelation(reaction, MIN);
                break;
            case MAX:
                LOGGER.debug("Processing Reaction for Global Minimum: ");
                delta = (int) calRelation(reaction, MAX);
                break;
            case MIXTURE:
                LOGGER.debug("Processing Reaction for Max-Mixture Model: ");
                delta = (int) calRelation(reaction, MIXTURE);
                break;
            case RINGS:
                LOGGER.debug("Processing Reaction for Ring Model: ");
                delta = (int) calRelation(reaction, RINGS);
                break;
        }
    }

    /**
     *
     * @return
     */
    public synchronized IReaction getMappedReaction() {
        return reaction;
    }

    private synchronized double calRelation(IReaction reaction, IMappingAlgorithm theory) {
        try {
            Map<Integer, IAtomContainer> educts
                    = synchronizedSortedMap(new TreeMap<>());
            for (int i = 0; i < reaction.getReactantCount(); i++) {
                educts.put(i, reaction.getReactants().getAtomContainer(i));
            }

            Map<Integer, IAtomContainer> products
                    = synchronizedSortedMap(new TreeMap<>());
            for (int i = 0; i < reaction.getProductCount(); i++) {
                products.put(i, reaction.getProducts().getAtomContainer(i));
            }

            GameTheoryMatrix EDSH
                    = new GameTheoryMatrix(theory, reaction, removeHydrogen);

            IGameTheory gameTheory = make(theory,
                    reaction,
                    removeHydrogen,
                    educts,
                    products,
                    EDSH);

            this.reactionBlastMolMapping = gameTheory.getReactionMolMapping();
            EDSH.Clear();

            return gameTheory.getDelta();
        } catch (Exception e) {
            LOGGER.error(e);
            return -1;
        }
    }

    /**
     * @return the delta
     */
    public synchronized int getDelta() {
        return delta;
    }

    /**
     * @return the reactionBlastMolMapping
     */
    public synchronized MoleculeMoleculeMapping getReactionBlastMolMapping() {
        return reactionBlastMolMapping;
    }
}
