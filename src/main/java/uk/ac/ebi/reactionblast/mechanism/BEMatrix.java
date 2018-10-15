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
import static java.lang.System.getProperty;
import static java.util.Collections.unmodifiableList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.E_OR_Z;
import static org.openscience.cdk.interfaces.IBond.Stereo.NONE;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_OR_DOWN;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import uk.ac.ebi.reactionblast.tools.ValencyCalculator;

/**
 * This class create the BEMatrix of a set of molecule according to the
 * DU-Theory. (I.Ugi et al., J. Chem. Inf. Comput. Sci. 1994, 34, 3-16)
 *
 * @author Syed Asad Rahman<asad@ebi.ac.uk>
 * @author Lorenzo Baldacci {lorenzo@ebi.ac.uk|lbaldacc@csr.unibo.it}
 */
public class BEMatrix extends EBIMatrix implements Serializable {

    private static final long serialVersionUID = -1420740601548197863L;

    private IAtomContainerSet myMoleculeSet = null;
    private List<IBond> bonds = null;
    private List<IAtom> atomArray = null;
    private final boolean withoutH;
    private final Map<IAtom, IAtom> mappings;

    /**
     * Class constructor. Creates an empty BEMatrix. Generates the BEMatrix for
     * a given IAtomContainerSet. H atoms are not considered when creating the
     * matrix.
     *
     * @param molSet The IAtomContainerSet for which the BEMatrix is required
     * @param skipHydrogen
     * @param bonds
     * @param mappings AAM
     */
    public BEMatrix(boolean skipHydrogen,
            IAtomContainerSet molSet,
            List<IBond> bonds,
            Map<IAtom, IAtom> mappings) {
        super(0, 0);
        this.withoutH = skipHydrogen;
        this.atomArray = new LinkedList<>();
        this.myMoleculeSet = molSet;
        this.bonds = bonds;
        this.mappings = mappings;
    }

    /**
     *
     * @throws CDKException
     */
    void setMatrixAtoms() throws CDKException {
        //System.out.println("H " + withoutH);
        initMatrix(0.);
        atomArray.clear();
        for (IAtomContainer container : myMoleculeSet.atomContainers()) {
            for (IAtom atom : container.atoms()) {
                if (withoutH && atom.getSymbol().matches("H")) {
                    continue;
                }
                if (!mappings.containsKey(atom) && !mappings.containsValue(atom)) {
                    continue;
                }
                atomArray.add(atom);
            }
        }
        setMatrix();
    }

    /**
     * The methods returns the order of the bond that connects two given atoms.
     * It returns 0 if the atoms are not linked by a bond
     *
     * @param a1 The first IAtom atom
     * @param a2 The second IAtom atom
     * @return The bond order
     */
    public int getOrder(IAtom a1, IAtom a2) {
        return (int) getValue(getIndexOfAtomID(a1.getID()), getIndexOfAtomID(a2.getID()));
    }

    private void setMatrix() throws CDKException {
//        reSizeMatrix(atomArray.size(), atomArray.size());
        reSizeMatrix(atomArray.size() + 1, atomArray.size() + 1);
        //free valence electrons on the diagonal
        for (int i = 0; i < atomArray.size(); i++) {
            setValue(i, i, getFreeValenceElectrons(atomArray.get(i)));
        }
        for (int i = 0; i < atomArray.size(); i++) {
            for (int j = 0; j < atomArray.size(); j++) {
                if (i != j) {
                    setValue(i, j, getBondOrder(atomArray.get(i), atomArray.get(j)));
                    setValue(j, i, getBondOrder(atomArray.get(i), atomArray.get(j)));
                }
            }
        }
        //Setting lone pairs
        for (int i = 0; i < atomArray.size(); i++) {
            setValue(atomArray.size(), i, 100);
            setValue(i, atomArray.size(), 100);
        }
        setValue(atomArray.size(), atomArray.size(), 200);
    }

    /**
     * This method performs the sorting of the atoms of the BEMatrix according
     * to the ArrayList passed as parameter. The method pivots also the matrix
     * according to the new sorting.
     *
     * @param orderedAtomArray ArrayList containing the new order of atoms
     * @return Canonical atom index
     * @throws CDKException
     */
    public int[] orderAtomArray(List<IAtom> orderedAtomArray) throws CDKException {
        int[] canonicalIndex = new int[orderedAtomArray.size()];
        /*
         This condition was changed from exception to error by Asad 
         to accomodate unbalanced reactions
         */
        if (orderedAtomArray.size() != atomArray.size()) {
//            System.LOGGER.debug("The matrix has not been ordered, " + atomArray.size() + " !=" + orderedAtomArray.size());
            throw new CDKException("The matrix has not been ordered: " + atomArray.size() + " !=" + orderedAtomArray.size());
        }
        for (IAtom orderedAtom : orderedAtomArray) {
            if (getIndexOfAtomID(orderedAtom.getID()) == -1) {
                throw new CDKException("The matrix has not been ordered");
            }
        }
        for (int i = 0; i < orderedAtomArray.size(); i++) {
            int di = getIndexOfAtomID(orderedAtomArray.get(i).getID());
            if (di != i) {
                pivot(di, i);
            }
            canonicalIndex[i] = di;
        }
        return canonicalIndex;
    }

    private int getIndexOfAtomID(String atomID) {
        int ind = -1;
        for (int i = 0; i < atomArray.size(); i++) {
            if (atomArray.get(i).getID().equals(atomID)) {
                ind = i;
            }
        }
        return ind;
    }

    /**
     * Perform the pivoting of the BEMatrix switching also the atoms in position
     * i1th and i2th.
     *
     * @param i1 Row index of the pivoting
     * @param i2 Column index of the pivoting
     */
    @Override
    public void pivot(int i1, int i2) {
        //label pivot
        IAtom appA = atomArray.get(i1);
        atomArray.set(i1, atomArray.get(i2));
        atomArray.set(i2, appA);
        double appD = 0.0d;
        //column exchange
        for (int i = 0; i < getRowDimension(); i++) {
            appD = getValue(i, i1);
            setValue(i, i1, getValue(i, i2));
            setValue(i, i2, appD);
        }
        //row exchange
        for (int i = 0; i < getColumnDimension(); i++) {
            appD = getValue(i1, i);
            setValue(i1, i, getValue(i2, i));
            setValue(i2, i, appD);
        }
    }

    private double getFreeValenceElectrons(IAtom a) throws CDKException {
        double freeValEle = 0;
        for (int i = 0; i < myMoleculeSet.getAtomContainerCount(); i++) {
            IAtomContainer mol = myMoleculeSet.getAtomContainer(i);
            if (mol.contains(a)) {
                freeValEle = ValencyCalculator.getFreeValenceElectrons(mol, a, withoutH);
            }
        }
        return freeValEle;
    }

    /**
     *
     * @param a
     * @param b
     * @return
     * @throws CDKException
     */
    public double getBondOrder(IAtom a, IAtom b) throws CDKException {
        double bondOrder = 0;
        IBond bond;
        for (int i = 0; i < myMoleculeSet.getAtomContainerCount(); i++) {
            IAtomContainer m = myMoleculeSet.getAtomContainer(i);
            bond = m.getBond(a, b);
            if (bond != null) {
                return convertBondOrder(bond);
            }
        }
        return bondOrder;
    }

    /**
     *
     * @param a
     * @param b
     * @return
     * @throws CDKException
     */
    public int getBondStereo(IAtom a, IAtom b) throws CDKException {
        int bondOrder = 0;
        IBond bond = getBond(a, b);
        if (bond != null) {
            return convertBondStereo(bond);
        }

        return bondOrder;
    }

    /**
     *
     * @throws CDKException
     */
    public void setAromaticBond() throws CDKException, CDKException {
        for (int i = 0; i < myMoleculeSet.getAtomContainerCount(); i++) {
            IAtomContainer m = myMoleculeSet.getAtomContainer(i);
            aromatizeMolecule(m);
        }
    }

    /**
     * Returns the ArrayList containing the atoms of the BEMatrix
     *
     * @return An ArrayList containing the atoms of the BEMatrix
     */
    public List<IAtom> getAtoms() {
        return unmodifiableList(atomArray);
    }

    /**
     * Returns the atom at the position pos in [0,..].
     *
     * @param pos The position of the atom to be retrieved.
     * @return The atom at the position pos.
     * @throws CDKException
     */
    public IAtom getAtom(int pos) throws CDKException {
        if (pos >= atomArray.size()) {
            throw new CDKException("Passed index out of range");
        }
        return atomArray.get(pos);
    }

    /**
     *
     * @param a first atom
     * @param b second atom
     * @return bond between a and b
     */
    public IBond getBond(IAtom a, IAtom b) {
        IBond bond = null;
        for (IBond localBond : bonds) {
            if (localBond.contains(a) && localBond.contains(b)) {
                bond = localBond;
            }
        }
        return bond;
    }

    /**
     *
     * @param at
     * @return
     */
    public IAtomContainer getAtomContainer(IAtom at) {
        IAtomContainer retMol = null;
        for (int i = 0; i < myMoleculeSet.getAtomContainerCount(); i++) {
            IAtomContainer mol = myMoleculeSet.getAtomContainer(i);
            for (int j = 0; j < mol.getAtomCount(); j++) {
                IAtom ma = mol.getAtom(j);
                if (ma.getID().equals(at.getID())) {
                    retMol = mol;
                    break;
                }
            }
            if (retMol != null) {
                break;
            }
        }
        return retMol;
    }

    /**
     *
     * @param bond
     * @return
     */
    public double convertBondOrder(IBond bond) {
        double value;
        switch (bond.getOrder()) {
            case QUADRUPLE:
                value = 4.0;
                break;
            case TRIPLE:
                value = 3.0;
                break;
            case DOUBLE:
                value = 2.0;
                break;
            case SINGLE:
                value = 1.0;
                break;
            default:
                value = 1.0;
        }
        return value;
    }

    /**
     *
     * @param bond
     * @return
     */
    public int convertBondStereo(IBond bond) {
        int value;
        switch (bond.getStereo()) {
            case UP:
                value = 1;
                break;
            case UP_INVERTED:
                value = 1;
                break;
            case DOWN:
                value = 6;
                break;
            case DOWN_INVERTED:
                value = 6;
                break;
            case UP_OR_DOWN:
                value = 4;
                break;
            case UP_OR_DOWN_INVERTED:
                value = 4;
                break;
            case E_OR_Z:
                value = 3;
                break;
            default:
                value = 0;
        }
        return value;
    }

    /**
     *
     * @param stereoValue
     * @return
     */
    public IBond.Stereo convertStereo(int stereoValue) {
        IBond.Stereo stereo = NONE;

        switch (stereoValue) {
            case 1:
                // up bond
                stereo = UP;
                break;
            case 6:
                // down bond
                stereo = DOWN;
                break;
            case 0:
                // bond has no stereochemistry
                stereo = NONE;
                break;
            case 4:
                //up or down bond
                stereo = UP_OR_DOWN;
                break;
            case 3:
                //e or z undefined
                stereo = E_OR_Z;
                break;
            default:
                break;
        }

        return stereo;
    }

    /**
     * @return the bonds
     */
    public List<IBond> getBonds() {
        return unmodifiableList(bonds);
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");
        result.append(atomArray.size()).append(NEW_LINE);
        atomArray.stream().forEach((atom) -> {
            result.append(atom.getSymbol()).append(atom.getID()).append("\t");
        });
        result.append(NEW_LINE);
        for (int i = 0; i < this.getRowDimension(); i++) {
            for (int j = 0; j < this.getColumnDimension(); j++) {
                result.append(this.getValue(i, j)).append("\t");
            }
            result.append(NEW_LINE);
        }
        return result.toString();
    }

    @Override
    public synchronized Object clone() throws CloneNotSupportedException {
        return super.clone(); //To change body of generated methods, choose Tools | Templates.
    }
}
