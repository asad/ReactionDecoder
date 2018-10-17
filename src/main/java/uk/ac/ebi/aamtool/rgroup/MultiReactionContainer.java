/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.aamtool.rgroup;

import static uk.ac.ebi.aamtool.rgroup.ECRgroupFrequency.DEBUG;
import static java.lang.System.out;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.aromaticity.Aromaticity;
import static org.openscience.cdk.aromaticity.ElectronDonation.daylight;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static org.openscience.cdk.tools.manipulator.ReactionManipulator.getAllAtomContainers;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCircularFragment;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class MultiReactionContainer {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MultiReactionContainer.class);

    private final Set<ReactionGroup> reaction;
    private final String enzyme;
    private boolean RGroup;
    private final FingerprintType fp;

    MultiReactionContainer(String enzyme) {
        this.enzyme = enzyme;
        this.reaction = new TreeSet<>();
        this.RGroup = false;
        fp = new FingerprintType(new TreeSet<>(), new TreeSet<>(), new TreeSet<>());

    }

    public void addReaction(IReaction r, String name) {
        ReactionGroup rg = new ReactionGroup(name);
        if (!reaction.contains(rg)) {
            calculateCommonFingerprint(r);

            boolean local_r_group_finder = false;
            List<IAtomContainer> allAtomContainers = getAllAtomContainers(r);
            for (IAtomContainer a : allAtomContainers) {
                if (isRGroupPresent(a)) {
                    local_r_group_finder = true;
                    break;
                }
            }
            rg.setRGroupPresent(local_r_group_finder);
            if (!isRGroup() && local_r_group_finder) {
                RGroup = true;
            }
            getReaction().add(rg);
        }
    }

    private boolean isRGroupPresent(IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            if (a instanceof IPseudoAtom) {
                if (a.getSymbol().contains("R")) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @return the enzyme
     */
    public String getEnzyme() {
        return enzyme;
    }

    /**
     * @return the enzyme
     */
    public String getEnzyme1Level() {
        return enzyme.split("\\.")[0];
    }

    /**
     * @return the enzyme
     */
    public String getEnzyme2Level() {
        return enzyme.split("\\.")[1];
    }

    /**
     * @return the enzyme
     */
    public String getEnzyme3Level() {
        return enzyme.split("\\.")[2];
    }

    /**
     * @return the enzyme
     */
    public String getEnzyme4Level() {
        return enzyme.split("\\.")[3];
    }

    /**
     * @return the RGroup
     */
    public boolean isRGroup() {
        return RGroup;
    }

    private void calculateCommonFingerprint(IReaction reaction) {
        Set<String> l = new HashSet<>();
        Set<String> r = new HashSet<>();
        SmilesGenerator sm = new SmilesGenerator(
                SmiFlavor.Unique
                | SmiFlavor.UseAromaticSymbols
                | SmiFlavor.Stereo);

        Aromaticity aromaticity = new Aromaticity(daylight(),
                Cycles.or(Cycles.all(),
                        Cycles.or(Cycles.relevant(),
                                Cycles.essential())));

        for (IAtomContainer a : reaction.getReactants().atomContainers()) {
            IAtomContainer ac = null;
            ac = removeHydrogensExceptSingleAndPreserveAtomID(a);
            for (int i = 0; i < ac.getAtomCount(); i++) {
                try {
                    aromaticity.apply(ac);
                    IAtomContainer circularFragment = getCircularFragment(ac, i, 1);
                    String smiles = sm.create(circularFragment);
                    l.add(smiles);
                    getAllFP().add(smiles);
                    circularFragment = getCircularFragment(ac, i, 2);
                    smiles = sm.create(circularFragment);
                    l.add(smiles);
                    getAllFP().add(smiles);

                    circularFragment = getCircularFragment(ac, i, 3);
                    smiles = sm.create(circularFragment);
                    l.add(smiles);
                    getAllFP().add(smiles);

                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }

        for (IAtomContainer a : reaction.getProducts().atomContainers()) {
            IAtomContainer ac = null;
            ac = removeHydrogensExceptSingleAndPreserveAtomID(a);
            for (int i = 0; i < ac.getAtomCount(); i++) {
                try {
                    aromaticity.apply(ac);
                    IAtomContainer circularFragment = getCircularFragment(ac, i, 1);
                    String smiles = sm.create(circularFragment);
                    r.add(smiles);
                    getAllFP().add(smiles);

                    circularFragment = getCircularFragment(ac, i, 2);
                    smiles = sm.create(circularFragment);
                    r.add(smiles);
                    getAllFP().add(smiles);

                    circularFragment = getCircularFragment(ac, i, 3);
                    smiles = sm.create(circularFragment);
                    r.add(smiles);
                    getAllFP().add(smiles);

                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
        }
        Set<String> common = new HashSet<>(l);
        boolean intersection = common.retainAll(r);

        if (DEBUG) {
            out.println("intersection " + common);
        }

        Set<String> difference = new TreeSet<>(l);
        difference.addAll(r);
        boolean removeAll = difference.removeAll(common);
        if (DEBUG) {
            out.println("difference " + difference);
        }
        /*
         if no reaction is present then add all the patterns
         */
        if (getReaction().isEmpty()) {
            getCommonCommonFP().addAll(common);
        } else {
            getCommonCommonFP().retainAll(common);
        }

        /*
         if no reaction is present then add all the patterns
         */
        if (getReaction().isEmpty()) {
            getCommonDifferenceFP().addAll(difference);
        } else {
            getCommonDifferenceFP().retainAll(difference);
        }
    }

    /**
     * @return the commonommonCommonFP
     */
    public Set<String> getCommonCommonFP() {
        return fp.getCommonCommonFP();
    }

    /**
     * @return the reactions
     */
    public Set<ReactionGroup> getReaction() {
        return reaction;
    }

    /**
     * @return the reaction count
     */
    public int getReactionCount() {
        return reaction.size();
    }

    /**
     * @return the commonommonDifferencommoneFP
     */
    public Set<String> getCommonDifferenceFP() {
        return fp.getCommonDifferenceFP();
    }

    /**
     * @return the commonommonDifferencommoneFP
     */
    public Set<String> getAllFP() {
        return fp.getAllPatternsFP();
    }
}
