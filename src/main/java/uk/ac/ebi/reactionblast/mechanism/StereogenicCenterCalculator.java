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
package uk.ac.ebi.reactionblast.mechanism;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.stereo.IStereoAndConformation;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.E;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.EITHER;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.NONE;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.R;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.S;
import static uk.ac.ebi.reactionblast.stereo.IStereoAndConformation.Z;

/**
 * This class creates a Stereo matrix for a set of molecules loosely based on a
 * series of paper by Ugi. (I.Ugi et al., J. Chem. Inf. Comput. Sci. 1994, 34,
 * 3-16)
 *
 * Prelog and Helmchen. Basic Principles of the CIP-System and Proposals for a
 * Revision. Angewandte Chemie International Edition 21 (1982) 567-683 Perdih
 * and Rmzinger. Stereo-chemistry and Sequence Rules A Proposal for Modification
 * of Cahn-Ingold-Prelog System. Tetrahedron: Asymmetry Vol 5 (1994) 835-861
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class StereogenicCenterCalculator implements Serializable {

    private static final long serialVersionUID = -1420740601548197863L;
    //    /**
//     * Returns type of Chirality changes
//     *
//     * @param a
//     * @param b
//     * @return
//     */
//    public boolean isChiralityChange(InchiStereoAndConformation a, InchiStereoAndConformation b) {
//        if (a.equals(InchiStereoAndConformation.PLUS) && b.equals(InchiStereoAndConformation.MINUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.PLUS) && b.equals(InchiStereoAndConformation.UNKNOWN)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.PLUS) && b.equals(InchiStereoAndConformation.NONE)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.MINUS) && b.equals(InchiStereoAndConformation.PLUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.MINUS) && b.equals(InchiStereoAndConformation.UNKNOWN)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.MINUS) && b.equals(InchiStereoAndConformation.NONE)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.NONE) && b.equals(InchiStereoAndConformation.UNKNOWN)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.UNKNOWN) && b.equals(InchiStereoAndConformation.NONE)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.UNKNOWN) && b.equals(InchiStereoAndConformation.PLUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.UNKNOWN) && b.equals(InchiStereoAndConformation.MINUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.NONE) && b.equals(InchiStereoAndConformation.PLUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.NONE) && b.equals(InchiStereoAndConformation.MINUS)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.OPPOSITE) && b.equals(InchiStereoAndConformation.TOGETHER)) {
//            return true;
//        } else if (a.equals(InchiStereoAndConformation.TOGETHER) && b.equals(InchiStereoAndConformation.OPPOSITE)) {
//            return true;
//        }
//        return false;
//    }
//

    /**
     *
     * @param reaction
     * @param chirality2DCDK
     * @param chirality2DChemaxon
     * @param chirality3DChemaxon
     * @return
     */
    public synchronized List<StereoChange> compare(
            IReaction reaction,
            Map<IAtom, IStereoAndConformation> chirality2DCDK,
            Map<IAtom, IStereoAndConformation> chirality2DChemaxon,
            Map<IAtom, IStereoAndConformation> chirality3DChemaxon) {
        List<StereoChange> stereoChangeList = new ArrayList<>();
        List<IAtom> queryAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                queryAtoms.add(a);
            }
        }
        List<IAtom> targetAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                targetAtoms.add(a);
            }
        }

        queryAtoms.stream().forEach((IAtom atomQ) -> {
            targetAtoms.stream().filter((atomT) -> (atomQ.getID().equals(atomT.getID()) && !atomQ.getSymbol().equalsIgnoreCase("H"))).forEach((atomT) -> {
                IStereoAndConformation rAtom2DCDKStereo = chirality2DCDK.get(atomQ);
                IStereoAndConformation pAtom2DCDKStereo = chirality2DCDK.get(atomT);
                IStereoAndConformation rAtom2DStereo = chirality2DChemaxon.get(atomQ);
                IStereoAndConformation pAtom2DStereo = chirality2DChemaxon.get(atomT);
                IStereoAndConformation rAtom3DStereo = chirality3DChemaxon.get(atomQ);
                IStereoAndConformation pAtom3DStereo = chirality3DChemaxon.get(atomT);
                if (isStereogenicChange(rAtom2DCDKStereo, pAtom2DCDKStereo)) {
                    StereoChange sc = new StereoChange(rAtom2DCDKStereo, pAtom2DCDKStereo, atomQ, atomT);
                    stereoChangeList.add(sc);
                }
                if (isStereogenicChange(rAtom3DStereo, pAtom3DStereo)) {
                    StereoChange sc = new StereoChange(rAtom3DStereo, pAtom3DStereo, atomQ, atomT);
                    stereoChangeList.add(sc);
                }
            });
        });
        return stereoChangeList;
    }

    /**
     *
     * @param reaction
     * @param chirality2DCDK
     * @param chirality2DChemaxon
     * @return
     */
    public synchronized List<StereoChange> compare(IReaction reaction, Map<IAtom, IStereoAndConformation> chirality2DCDK, Map<IAtom, IStereoAndConformation> chirality2DChemaxon) {
        List<StereoChange> stereoChangeList = new ArrayList<>();
        List<IAtom> queryAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                queryAtoms.add(a);
            }
        }
        List<IAtom> targetAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                targetAtoms.add(a);
            }
        }

        queryAtoms.stream().forEach((IAtom atomQ) -> {
            targetAtoms.stream().filter((atomT) -> (atomQ.getID().equals(atomT.getID()) && !atomQ.getSymbol().equalsIgnoreCase("H"))).forEach((IAtom atomT) -> {
                IStereoAndConformation rAtom2DCDKStereo = chirality2DCDK.get(atomQ);
                IStereoAndConformation pAtom2DCDKStereo = chirality2DCDK.get(atomT);
                IStereoAndConformation rAtom2DStereo = chirality2DChemaxon.get(atomQ);
                IStereoAndConformation pAtom2DStereo = chirality2DChemaxon.get(atomT);
                //                    System.out.println("atomQ " + atomQ.getID() + " S: " + atomQ.getSymbol());
//                    System.out.println("atomT " + atomT.getID() + " S: " + atomT.getSymbol());
//
//                    System.out.println("atomQ " + chirality2DCDK.containsKey(atomQ));
//                    System.out.println("atomT " + chirality2DCDK.containsKey(atomT));
                if (isStereogenicChange(rAtom2DStereo, pAtom2DStereo)) {
                    StereoChange sc = new StereoChange(rAtom2DStereo, pAtom2DStereo, atomQ, atomT);
                    stereoChangeList.add(sc);
                }
                if (isStereogenicChange(rAtom2DCDKStereo, pAtom2DCDKStereo)) {
                    StereoChange sc = new StereoChange(rAtom2DCDKStereo, pAtom2DCDKStereo, atomQ, atomT);
                    stereoChangeList.add(sc);
                }
            });
        });
        return stereoChangeList;
    }

    /**
     *
     * @param reaction
     * @param chirality2DCDK
     * @return
     */
    public synchronized List<StereoChange> compare(IReaction reaction, Map<IAtom, IStereoAndConformation> chirality2DCDK) {

        List<StereoChange> stereoChangeList = new ArrayList<>();
        List<IAtom> queryAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                queryAtoms.add(a);
            }
        }
        List<IAtom> targetAtoms = new ArrayList<>();
        for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
            for (IAtom a : ac.atoms()) {
                targetAtoms.add(a);
            }
        }
        queryAtoms.forEach((IAtom atomQ) -> {
            targetAtoms.stream().filter((atomT) -> (atomQ.getID().equals(atomT.getID()) && !atomQ.getSymbol().equalsIgnoreCase("H"))).forEachOrdered((atomT) -> {
                IStereoAndConformation rAtom2DCDKStereo = chirality2DCDK.get(atomQ);
                IStereoAndConformation pAtom2DCDKStereo = chirality2DCDK.get(atomT);
//                    System.out.println("atomQ " + atomQ.getID() + " S: " + atomQ.getSymbol());
//                    System.out.println("atomT " + atomT.getID() + " S: " + atomT.getSymbol());
//
//                    System.out.println("atomQ " + chirality2DCDK.containsKey(atomQ));
//                    System.out.println("atomT " + chirality2DCDK.containsKey(atomT));
                if (isStereogenicChange(rAtom2DCDKStereo, pAtom2DCDKStereo)) {
                    StereoChange sc = new StereoChange(rAtom2DCDKStereo, pAtom2DCDKStereo, atomQ, atomT);
                    stereoChangeList.add(sc);
                }
            });
        });
        return stereoChangeList;
    }

    /**
     * Returns type of stereo changes
     *
     * @param a
     * @param b
     * @return
     */
    public synchronized boolean isStereogenicChange(IStereoAndConformation a, IStereoAndConformation b) {
        if (a.equals(S) && b.equals(NONE)) {
            return true;
        } else if (a.equals(R) && b.equals(NONE)) {
            return true;
        } else if (b.equals(S) && a.equals(NONE)) {
            return true;
        } else if (b.equals(R) && a.equals(NONE)) {
            return true;
        } else if (a.equals(R) && b.equals(S)) {
            return true;
        } else if (a.equals(S) && b.equals(R)) {
            return true;
        } else if (a.equals(S) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(R) && b.equals(EITHER)) {
            return true;
        } else if (b.equals(S) && a.equals(EITHER)) {
            return true;
        } else if (b.equals(R) && a.equals(EITHER)) {
            return true;
        } else if (a.equals(EITHER) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(NONE) && b.equals(EITHER)) {
            return true;
        } else if (a.equals(EITHER) && b.equals(NONE)) {
            return true;
        } else if (a.equals(Z) && b.equals(E)) {
            return true;
        } else if (a.equals(E) && b.equals(Z)) {
            return true;
        }
        return false;
    }
}
