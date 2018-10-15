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
package uk.ac.ebi.reactionblast.stereo.ebi;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.centres.cdk.CDKPerceptor;
import uk.ac.ebi.centres.descriptor.Planar;
import uk.ac.ebi.centres.descriptor.Tetrahedral;
import uk.ac.ebi.centres.descriptor.Trigonal;
import uk.ac.ebi.reactionblast.mechanism.helper.Utility;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 * Tool for comparing chiralities.
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 */
public abstract class StereoCenteralityTool extends Utility {

    private static final long serialVersionUID = 17867606807697859L;

    private static IAtom getAtomByID(String id, IAtomContainer ac) {
        for (IAtom a : ac.atoms()) {
            if (a.getID().equals(id)) {
                return a;
            }
        }
        return null;
    }

    /**
     * This Chirality is based on the 2D with stereo code written by John May in
     * our collaboration. Note: Explicit Hydrogens should be added before
     * calling.
     *
     * @param reaction
     * @return
     * @throws CDKException
     * @throws java.lang.CloneNotSupportedException
     */
    public static Map<IAtom, IStereoAndConformation> getChirality2D(IReaction reaction) throws CDKException, CloneNotSupportedException {
        Map<IAtom, IStereoAndConformation> chiralityMap = new HashMap<>();
        CDKPerceptor perceptor = new CDKPerceptor();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
//            System.LOGGER.debug("R 2D CDK based stereo perception for " + ac.getID());
            Map<IAtom, IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH, perceptor);
//            System.LOGGER.debug("R 2D CDK based stereo " + chirality2D.size());
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().stream().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            IAtomContainer containerWithoutH = removeHydrogensExceptSingleAndPreserveAtomID(ac);
//            System.LOGGER.debug("P 2D CDK based stereo perception for " + ac.getID());
            Map<IAtom, IStereoAndConformation> chirality2D = getChirality2D(containerWithoutH, perceptor);
//            System.LOGGER.debug("P 2D CDK based stereo " + chirality2D.size());
            if (!chirality2D.isEmpty()) {
                chirality2D.entrySet().stream().forEach((m) -> {
                    IAtom atomByID = getAtomByID(m.getKey().getID(), ac);
                    if (atomByID != null) {
                        atomByID.setProperty("Stereo", m.getValue());
                        chiralityMap.put(atomByID, m.getValue());
                    }
                });
            }
        }
        return chiralityMap;
    }

    /**
     *
     * @param ac
     * @param perceptor
     * @return
     */
    public static Map<IAtom, IStereoAndConformation> getChirality2D(IAtomContainer ac, CDKPerceptor perceptor) {
        Map<IAtom, IStereoAndConformation> chiralityMap = new HashMap<>();
        perceptor.perceive(ac);
        for (IAtom atom : ac.atoms()) {
            if (Tetrahedral.R.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.R);
            } else if (Tetrahedral.S.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.S);
            } else if (Planar.E.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.E);
            } else if (Planar.Z.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.Z);
            } else if (Trigonal.Re.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.P);
            } else if (Trigonal.Si.equals(atom.getProperty("descriptor"))) {
                chiralityMap.put(atom, IStereoAndConformation.M);
            } else {
                chiralityMap.put(atom, IStereoAndConformation.NONE);
            }
        }
        for (IBond bond : ac.bonds()) {
            if (Planar.E.equals(bond.getProperty("descriptor"))) {
                chiralityMap.put(bond.getAtom(0), IStereoAndConformation.E);
                chiralityMap.put(bond.getAtom(1), IStereoAndConformation.E);
            } else if (Planar.Z.equals(bond.getProperty("descriptor"))) {
                chiralityMap.put(bond.getAtom(0), IStereoAndConformation.Z);
                chiralityMap.put(bond.getAtom(1), IStereoAndConformation.Z);
            } else if (Trigonal.Re.equals(bond.getProperty("descriptor"))) {
                chiralityMap.put(bond.getAtom(0), IStereoAndConformation.P);
                chiralityMap.put(bond.getAtom(1), IStereoAndConformation.P);
            } else if (Trigonal.Si.equals(bond.getProperty("descriptor"))) {
                chiralityMap.put(bond.getAtom(0), IStereoAndConformation.M);
                chiralityMap.put(bond.getAtom(1), IStereoAndConformation.M);
            }
        }

        chiralityMap.keySet().stream().forEach((atom) -> {
            atom.setProperty("Stereo", chiralityMap.get(atom));
        });
        return chiralityMap;
    }
}
