package org.openscience.smsd.tools;

/**
 *
 * Copyright (C) 2006-2016 Syed Asad Rahman {asad@ebi.ac.uk}
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received atom copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
import java.io.Serializable;
import java.util.HashMap;
import java.util.Objects;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.suppressHydrogens;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

/**
 * Class that handles some customised features for atom containers.
 * <p>
 * This is an extension of CDK GraphAtomContainer. Some part of this code was
 * taken from CDK source code and modified.</p>
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class ExtAtomContainerManipulator extends AtomContainerManipulator implements Serializable {

    static final Logger logger = Logger.getLogger(ExtAtomContainerManipulator.class.getName());
    static final long serialVersionUID = 1786786539472837495L;

    /**
     * Create an copy of the {@code org} structure with explicit hydrogens
     * removed. Stereochemistry is updated but up and down bonds in a depiction
     * may need to be recalculated (see. StructureDiagramGenerator).
     *
     * @param org The AtomContainer from which to remove the hydrogens
     * @return The molecule without hydrogens.
     * @cdk.keyword hydrogens, removal, suppress
     * @see #copyAndSuppressedHydrogens
     */
    public static IAtomContainer removeHydrogensExceptSingleAndPreserveAtomID(IAtomContainer org) {
        return copyAndSuppressedHydrogens(org);
    }

    /**
     * Copy the input container and suppress any explicit hydrogens. Only
     * hydrogens that can be represented as a hydrogen count value on the atom
     * are suppressed. If a copy is not needed please use {@link
     * #suppressHydrogens}.
     *
     * @param org the container from which to remove hydrogens
     * @return a copy of the input with suppressed hydrogens
     * @see #suppressHydrogens
     */
    public static IAtomContainer copyAndSuppressedHydrogens(IAtomContainer org) {

        /*
         * @ASAD: IMP STEP to avoid unset Hydrogen arror:
         * Set implicit Hydrogen count
         */
        try {
            IAtomContainer clone = org.clone();
            clone.setID(org.getID());
            for (int i = 0; i < org.getAtomCount(); i++) {
                int implicitHydrogenCount = getImplicitHydrogenCount(clone.getAtom(i));
                clone.getAtom(i).setImplicitHydrogenCount(implicitHydrogenCount);
                clone.getAtom(i).setProperties(org.getAtom(i).getProperties());
                clone.getAtom(i).setFlags(org.getAtom(i).getFlags());
            }
            return suppressHydrogens(clone);
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException("atom container could not be cloned");
        }
    }

    private static void printAtoms(IAtomContainer mol) {
        System.out.print("Atom: ");
        for (IAtom a : mol.atoms()) {

            System.out.print(a.getSymbol());
            System.out.print("[" + a.getFormalCharge() + "]");
            if (a.getID() != null) {
                System.out.print("[" + a.getID() + "]");
            }

        }
        System.out.println();
        System.out.println();
    }

    /**
     * Modules for cleaning a molecule
     *
     * @param molecule_orignal
     * @return cleaned GraphAtomContainer
     */
    public synchronized static IAtomContainer checkAndCleanMolecule(IAtomContainer molecule_orignal) {
        boolean isMarkush = false;
        IAtomContainer molecule = molecule_orignal;
        for (IAtom atom : molecule.atoms()) {
            if (atom.getSymbol().equals("R")) {
                isMarkush = true;
                break;
            }
        }

        if (isMarkush) {
            logger.log(Level.WARNING, "Skipping Markush structure for sanity check");
        }

        // Check for salts and such
        if (!ConnectivityChecker.isConnected(molecule)) {
            // lets see if we have just two parts if so, we assume its a salt and just work
            // on the larger part. Ideally we should have a check to ensure that the smaller
            //  part is a metal/halogen etc.
            IAtomContainerSet fragments = ConnectivityChecker.partitionIntoMolecules(molecule);
            if (fragments.getAtomContainerCount() > 2) {
                logger.log(Level.WARNING, "More than 2 components. Skipped");
            } else {
                IAtomContainer frag1 = fragments.getAtomContainer(0);
                IAtomContainer frag2 = fragments.getAtomContainer(1);
                if (frag1.getAtomCount() > frag2.getAtomCount()) {
                    molecule = frag1;
                } else {
                    molecule = frag2;
                }
            }
        }
        aromatizeMolecule(molecule);
        return molecule;
    }

    /**
     * This function finds rings and uses aromaticity detection code to
     * aromatize the molecule.
     *
     * @param mol input molecule
     */
    public static void aromatizeMolecule(IAtomContainer mol) {
        try {
            // need to find rings and aromaticity again since added H's
            IRingSet ringSet = null;
            try {
                AllRingsFinder arf = new AllRingsFinder();
                ringSet = arf.findAllRings(mol);
                RingSetManipulator.markAromaticRings(ringSet);
            } catch (CDKException e) {
                logger.log(Level.WARNING, "Error in find and assigning rings in the molecule. ", mol.getID());
            }

            try {
                try {
                    // figure out which atoms are in aromatic rings:
                    percieveAtomTypesAndConfigureAtoms(mol);
                    aromatizeCDK(mol);
                } catch (CDKException e) {
                    aromatizeDayLight(mol);
                }
            } catch (CDKException e) {
                logger.log(Level.WARNING, "Error in aromaticity dectection. ", mol.getID());
            }

            if (ringSet == null) {
                return;
            }

            // only atoms in 6 membered rings are aromatic
            // determine largest ring that each atom is atom part of
            for (int i = 0; i <= mol.getAtomCount() - 1; i++) {

                mol.getAtom(i).setFlag(CDKConstants.ISAROMATIC, false);

                jloop:
                for (int j = 0; j <= ringSet.getAtomContainerCount() - 1; j++) {
                    //logger.debug(i+"\t"+j);
                    IRing ring = (IRing) ringSet.getAtomContainer(j);
                    if (!ring.getFlag(CDKConstants.ISAROMATIC)) {
                        continue jloop;
                    }

                    boolean haveatom = ring.contains(mol.getAtom(i));

                    //logger.debug("haveatom="+haveatom);
                    if (haveatom && ring.getAtomCount() == 6) {
                        mol.getAtom(i).setFlag(CDKConstants.ISAROMATIC, true);
                    }
                }
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, "Aromaticity detection failed for molecule. ", mol.getID());
        }
    }

    /**
     * Returns deep copy of the molecule
     *
     * @param container
     * @return deep copy of the mol
     * @throws CloneNotSupportedException
     */
    public static IAtomContainer cloneWithIDs(IAtomContainer container) throws CloneNotSupportedException {
        //setNullHCountToZero(container);
        IAtomContainer ac = new AtomContainer(container).clone();/*Set IDs as CDK clone doesn't*/
        for (int i = 0; i < ac.getAtomCount(); i++) {
            ac.getAtom(i).setID(container.getAtom(i).getID());
            if (ac.getAtom(i).getProperties() == null) {
                ac.getAtom(i).setProperties(new HashMap<>());
            }
        }

        for (int i = 0; i < ac.getBondCount(); i++) {
            if (ac.getBond(i).getProperties() == null) {
                ac.getBond(i).setProperties(new HashMap<>());
            }
        }

        ac.setProperties(container.getProperties());
        ac.setFlags(container.getFlags());
        ac.setID(container.getID());
        ac.notifyChanged();
        return ac;
    }

    /**
     * Returns deep copy of the molecule
     *
     * @param container
     * @return deep copy of the mol
     * @throws CloneNotSupportedException
     */
    public static IAtomContainer newInstanceWithIDs(IAtomContainer container) throws CloneNotSupportedException {
        //setNullHCountToZero(container);
        IAtomContainer ac = container.getBuilder().newInstance(IAtomContainer.class, container);
        /*Set IDs as CDK clone doesn't*/
        for (int i = 0; i < ac.getAtomCount(); i++) {
            ac.getAtom(i).setID(container.getAtom(i).getID());
            if (ac.getAtom(i).getProperties() == null) {
                ac.getAtom(i).setProperties(new HashMap<>());
            }
        }

        for (int i = 0; i < ac.getBondCount(); i++) {
            if (ac.getBond(i).getProperties() == null) {
                ac.getBond(i).setProperties(new HashMap<>());
            }
        }
        ac.setProperties(container.getProperties());
        ac.setFlags(container.getFlags());
        ac.setID(container.getID());
        ac.notifyChanged();
        return ac;
    }

    /**
     * Returns The number of explicit hydrogens for a given IAtom.
     *
     * @param atomContainer
     * @param atom
     * @return The number of explicit hydrogens on the given IAtom.
     */
    public static int getExplicitHydrogenCount(IAtomContainer atomContainer, IAtom atom) {
        int hCount = 0;
        hCount = atomContainer.getConnectedAtomsList(atom).stream().map((iAtom) -> iAtom).filter((connectedAtom)
                -> (connectedAtom.getSymbol().equals("H"))).map((_item) -> 1).reduce(hCount, Integer::sum);
        return hCount;
    }

    /**
     * Returns The number of Implicit Hydrogen Count for a given IAtom.
     *
     * @param atom
     * @return Implicit Hydrogen Count
     */
    public static int getImplicitHydrogenCount(IAtom atom) {
        return Objects.equals(atom.getImplicitHydrogenCount(), CDKConstants.UNSET) ? 0 : atom.getImplicitHydrogenCount();
    }

    /**
     * The summed implicit + explicit hydrogens of the given IAtom.
     *
     * @param atomContainer
     * @param atom
     * @return The summed implicit + explicit hydrogens of the given IAtom.
     */
    public static int getHydrogenCount(IAtomContainer atomContainer, IAtom atom) {
        return getExplicitHydrogenCount(atomContainer, atom) + getImplicitHydrogenCount(atom);
    }

    /**
     * Returns IAtomContainer without Hydrogen. If an GraphAtomContainer has
     * atom single atom which is atom Hydrogen then its not removed.
     *
     * @param atomContainer
     * @return IAtomContainer without Hydrogen. If an GraphAtomContainer has
     * atom single atom which is atom Hydrogen then its not removed.
     */
    public static IAtomContainer convertExplicitToImplicitHydrogens(IAtomContainer atomContainer) {
        IAtomContainer mol = atomContainer.getBuilder().newInstance(IAtomContainer.class, atomContainer);
        if (mol.getAtomCount() > 1) {
            mol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
        } else if (mol.getAtomCount() == 1 && !atomContainer.atoms().iterator().next().getSymbol().equalsIgnoreCase("H")) {
            /*
             * Pseudo-atoms and unconfigured atoms safetynet
             */
            convertImplicitToExplicitHydrogens(mol);
            mol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
        } else if (mol.getAtomCount() == 1 && atomContainer.atoms().iterator().next().getSymbol().equalsIgnoreCase("H")) {
            System.err.println("WARNING: single hydrogen atom removal not supported!");
        }
        mol.setProperties(atomContainer.getProperties());
        mol.setFlags(atomContainer.getFlags());
        if (atomContainer.getID() != null) {
            mol.setID(atomContainer.getID());
        }
        return mol;
    }

    /**
     * Convenience method to perceive atom types for all <code>IAtom</code>s in
     * the <code>IAtomContainer</code>, using the
     * <code>CDKAtomTypeMatcher</code>. If the matcher finds atom matching atom
     * type, the <code>IAtom</code> will be configured to have the same
     * properties as the <code>IAtomType</code>. If no matching atom type is
     * found, no configuration is performed.
     *
     * @param container
     * @throws CDKException
     */
    public static void percieveAtomTypesAndConfigureAtoms(IAtomContainer container) throws CDKException {
        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        for (IAtom atom : container.atoms()) {
            if (!(atom instanceof IPseudoAtom)) {
                try {
                    IAtomType matched = matcher.findMatchingAtomType(container, atom);
                    if (matched != null) {
                        AtomTypeManipulator.configure(atom, matched);
                    }
                } catch (CDKException e) {
                    logger.log(Level.WARNING,
                            "Failed to find Matching AtomType! {0}{1}", new Object[]{atom.getSymbol(), e});
                }
            }
        }
    }

    /**
     *
     * @param gMol
     * @return a new mol with explicit Hydrogens
     * @throws CloneNotSupportedException
     */
    public static IAtomContainer addExplicitH(IAtomContainer gMol) throws CloneNotSupportedException {
        IAtomContainer mol = cloneWithIDs(gMol);
//        fixDativeBonds(mol);
        CDKHydrogenAdder hydAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
        try {
            percieveAtomTypesAndConfigureAtoms(mol);
            for (IAtom a : mol.atoms()) {
                if (!(a instanceof IPseudoAtom)) {
                    try {
                        hydAdder.addImplicitHydrogens(mol, a);
                    } catch (Exception e) {
                        a.setImplicitHydrogenCount(0);
                        System.err.println("WARNING: Error in adding Hydrogen" + ":" + a.getSymbol());
                        logger.log(Level.WARNING, "This might effect the final calculations!");
                    }
                } else {
                    a.setImplicitHydrogenCount(0);
                }
            }
        } catch (CDKException ex) {
            logger.log(Level.SEVERE, null, ex);
        }
        convertImplicitToExplicitHydrogens(mol);
        return mol;

    }

    /**
     *
     * @param molecule
     * @throws CDKException
     */
    public static void aromatizeDayLight(IAtomContainer molecule) throws CDKException {
        ElectronDonation model = ElectronDonation.daylight();
        CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.relevant());
        Aromaticity aromaticity = new Aromaticity(model, cycles);
        try {
            aromaticity.apply(molecule);
        } catch (CDKException e) {
            logger.log(Level.WARNING,
                    "Aromaticity detection failed due to presence of unset "
                    + "atom hybridisation", molecule.getID());
        }
    }

    /**
     *
     * @param molecule
     * @throws CDKException
     */
    public static void aromatizeCDK(IAtomContainer molecule) throws CDKException {
        ElectronDonation model = ElectronDonation.cdk();
        CycleFinder cycles = Cycles.cdkAromaticSet();
        Aromaticity aromaticity = new Aromaticity(model, cycles);
        ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        try {
            aromaticity.apply(molecule);
        } catch (CDKException e) {
            logger.log(Level.WARNING,
                    "Aromaticity detection failed due to presence of unset "
                    + "atom hybridisation", molecule.getID());
        }
    }

    /**
     * This method is a workaround by assigning dative bonds to single
     *
     * @param mol
     */
    public static void fixDativeBonds(IAtomContainer mol) {
        if (!(mol instanceof IQueryAtomContainer)) {
            for (IBond bond : mol.bonds()) {
                if (bond.getOrder() == IBond.Order.UNSET) {
                    bond.setOrder(IBond.Order.SINGLE);
                }
            }
        }
    }

    /**
     * Set all null hydrogen counts to 0. Generally hydrogen counts are present
     * and if not we add them. However the molecule being tested can't include
     * hydrogen counts as then fingerprints don't line up (substructure
     * filtering). The previous behaviour of the SMARTS matching was to treat
     * null hydrogens as 0 - the new behaviour is to complain about it.
     *
     * @param mol molecule to zero out hydrogen counts
     */
    static void setNullHCountToZero(IAtomContainer mol) {
        for (IAtom a : mol.atoms()) {
            if (a.getImplicitHydrogenCount() == null) {
                a.setImplicitHydrogenCount(0);
            }
        }
    }
}
