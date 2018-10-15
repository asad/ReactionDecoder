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
import static java.lang.Math.abs;
import static java.lang.System.getProperty;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.mechanism.helper.AtomAtomMappingContainer;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;

/**
 * This class create the RMatrix of a reaction according to the DU-Theory.
 * (I.Ugi et al., J. Chem. Inf. Comput. Sci. 1994, 34, 3-16).
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author Lorenzo Baldacci {lorenzo@ebi.ac.uk|lbaldacc@csr.unibo.it}
 */
public final class RMatrix extends EBIMatrix implements Serializable {

    private static final long serialVersionUID = 7057060562283378684L;
    private static final ILoggingTool LOGGER = createLoggingTool(RMatrix.class);

    private BEMatrix reactantBEMatrix = null;
    private BEMatrix productBEMatrix = null;
    private AtomAtomMappingContainer myMapping = null;
    private boolean DEBUG = false;

    /**
     * Class constructor. Generates the RMatrix of a reaction given the
     * BEMatrices of Reactants and Products.
     *
     * @param reactantBE Reactants BEMatrix
     * @param productBE Products BEMatrix
     * @param mapping Atom-Atom mappings between reactant and products atoms
     * @throws CDKException
     */
    public RMatrix(BEMatrix reactantBE, BEMatrix productBE, AtomAtomMappingContainer mapping) throws CDKException {
        super(reactantBE.getRowDimension(), reactantBE.getRowDimension());
        /*
         * Asad Commented this part to check mapping with Hydrogen and partial mapping
         */
//        if (reactantBE.getRowDimension() != productBE.getRowDimension()) {
////            throw new CDKException("This reaction is unbalanced. The bond change(s) "
////                    + "information in the R-Matrix may be baised.");
//            throw new CDKException("Either this reaction is unbalanced or mapping(s) are incorrect. "
//                    + "\nThe bond change(s)information may be biased.");
//        }
        /*
         Comment this to allow mapping of reactions even if the atoms are unmapped
         */
        int expectedOverlap = countAtomOverlap(reactantBE.getAtoms(), productBE.getAtoms());
        if (DEBUG) {
            out.println("expectedOverlap " + expectedOverlap + ", " + mapping.getSize());
        }
        if (expectedOverlap != mapping.getSize()) {
            LOGGER.debug("Core Reactant Atoms: " + (reactantBE.getRowDimension() - 1));
            LOGGER.debug("Core Product Atoms: " + (productBE.getRowDimension() - 1));
            LOGGER.debug("Mapping Atoms: " + mapping.getSize());
            throw new CDKException("Unable to construct a reaction matrix; "
                    + ((productBE.getRowDimension() - (mapping.getSize() + 1)))
                    + " atom(s) remain unmapped!.");
        }
        initMatrix(0.);

        this.reactantBEMatrix = reactantBE;
        this.productBEMatrix = productBE;
        this.myMapping = mapping;
        try {
            reactantBEMatrix.setAromaticBond();
            productBEMatrix.setAromaticBond();
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        ArrayList<IAtom> orderedBEMatrixAtomArray = new ArrayList<>();

        for (int i = 0; i < reactantBEMatrix.getRowDimension() - 1; i++) {
            IAtom atomR = reactantBEMatrix.getAtom(i);
            IAtom mappedProductAtom = myMapping.getMappedProductAtom(atomR);
            if (mappedProductAtom != null) {
                orderedBEMatrixAtomArray.add(mappedProductAtom);
            }
        }

        int[] canonicalOrderedAtomArray = productBEMatrix.orderAtomArray(orderedBEMatrixAtomArray);
        for (int i = 0; i < getMappedAtomCount(); i++) {
            String p_id_I = productBEMatrix.getAtom(i).getID();
            String r_id_I = reactantBEMatrix.getAtom(i).getID();
            for (int j = 0; j < getMappedAtomCount(); j++) {
                String p_id_J = productBEMatrix.getAtom(j).getID();
                String r_id_J = reactantBEMatrix.getAtom(j).getID();
                /*
                 Match ids for unbalanced reactions
                 */
                if (r_id_I.equals(p_id_I) && r_id_J.equals(p_id_J)) {
                    double value = productBEMatrix.getValue(i, j) - reactantBEMatrix.getValue(i, j);
                    boolean aromaticFlag = isAromaticChange(i, j);
                    if (aromaticFlag && value != 0.0) {
                        super.setValue(i, j, 0.0);
                    } else {
                        super.setValue(i, j, value);
                    }
                }
            }
        }
        if (DEBUG) {
            System.out.println("BE-React " + reactantBE.toString());
            System.out.println("BE-Prod " + productBE.toString());
            System.out.println("R " + toString());
        }
    }

    private synchronized boolean isAromaticChange(int IndexI, int IndexJ) throws CDKException {

        IAtom ra1 = getReactantBEMatrix().getAtom(IndexI);
        IAtom pa1 = getProductBEMatrix().getAtom(IndexI);
        IAtom ra2 = getReactantBEMatrix().getAtom(IndexJ);
        IAtom pa2 = getProductBEMatrix().getAtom(IndexJ);

        IBond rb = getReactantBEMatrix().getAtomContainer(ra1).getBond(ra1, ra2);
        IBond pb = getProductBEMatrix().getAtomContainer(pa1).getBond(pa1, pa2);
        if (rb != null && pb != null) {
            if ((rb.getFlag(ISINRING) && pb.getFlag(ISINRING))
                    && (rb.getFlag(ISAROMATIC) && pb.getFlag(ISAROMATIC))) {
                return true;
            }
        }
        return false;

    }

    /**
     * Returns the product atom in the idx-th position
     *
     * @param idx the position for which the product atom is needed
     * @return the product atom in the idx-th position, returns null if the idx
     * index is out of bounds.
     * @throws CDKException
     */
    public synchronized IAtom getProductAtom(int idx) throws CDKException {
        IAtom ret = null;
        if ((idx < getProductBEMatrix().getRowDimension()) && (idx > -1)) {
            ret = getProductBEMatrix().getAtom(idx);
        }
        return ret;
    }

    /**
     * Returns the reactant atom in the idx-th position
     *
     * @param idx the position for which the reactant atom is needed
     * @return the reactant atom in the idx-th position, returns null if the idx
     * index is out of bounds.
     * @throws CDKException
     */
    public synchronized IAtom getReactantAtom(int idx) throws CDKException {
        IAtom ret = null;
        if ((idx < getReactantBEMatrix().getRowDimension()) && (idx > -1)) {
            ret = getReactantBEMatrix().getAtom(idx);
        }
        return ret;
    }

    /**
     * Returns the ArrayList containing the reactant atoms of the RMatrix
     *
     * @return The ArrayList containing the reactant atoms of the RMatrix
     */
    public synchronized List<IAtom> getReactantsAtomArray() {
        return getReactantBEMatrix().getAtoms();
    }

    /**
     * Returns the ArrayList containing the product atoms of the RMatrix
     *
     * @return The ArrayList containing the product atoms of the RMatrix
     */
    public synchronized List<IAtom> getProductsAtomArray() {
        return getProductBEMatrix().getAtoms();
    }

    /**
     *
     * @param atomID1
     * @param atomID2
     * @return
     * @throws CDKException
     */
    public synchronized int getValueByReactantAtoms(String atomID1, String atomID2) throws CDKException {
        int res = 0;
        for (int i = 0; i < getRowDimension() - 1; i++) {
            for (int j = 0; j < getColumnDimension() - 1; j++) {
                if (getReactantBEMatrix().getAtom(i).getID().equals(atomID1)
                        && getReactantBEMatrix().getAtom(j).getID().equals(atomID2)) {
                    res = (int) getValue(i, j);
                }
            }
        }
        return res;
    }

    /**
     *
     * @param atomID1
     * @param atomID2
     * @return
     * @throws CDKException
     */
    public synchronized int getValueByProductAtoms(String atomID1, String atomID2) throws CDKException {
        int res = 0;
        for (int i = 0; i < getRowDimension() - 1; i++) {
            for (int j = 0; j < getColumnDimension() - 1; j++) {
                if (getProductBEMatrix().getAtom(i).getID().equals(atomID1)
                        && getProductBEMatrix().getAtom(j).getID().equals(atomID2)) {
                    res = (int) getValue(i, j);
                }
            }
        }
        return res;
    }

    /**
     *
     * @return
     */
    public synchronized int getAbsChanges() {
        int acc = 0;
        for (int i = 0; i < getRowDimension(); i++) {
            for (int j = 0; j < getColumnDimension(); j++) {
                acc += abs((int) getValue(i, j));
            }
        }
        return acc;
    }

    /**
     *
     * @return
     */
    protected synchronized int getMappedAtomCount() {
        return getMyMapping().getSize();
    }

    /**
     *
     * @return Atom count without Hydrogens
     */
    protected synchronized int getAtomCountWithoutHydrogens() {
        return getMyMapping().getSizeNoHydrogens();
    }

    /**
     * @return the reactantBEMatrix
     */
    public synchronized BEMatrix getReactantBEMatrix() {
        return reactantBEMatrix;
    }

    /**
     * @param reactantBEMatrix the reactantBEMatrix to set
     */
    public synchronized void setReactantBEMatrix(BEMatrix reactantBEMatrix) {
        this.reactantBEMatrix = reactantBEMatrix;
    }

    /**
     * @return the productBEMatrix
     */
    public synchronized BEMatrix getProductBEMatrix() {
        return productBEMatrix;
    }

    /**
     * @param productBEMatrix the productBEMatrix to set
     */
    public synchronized void setProductBEMatrix(BEMatrix productBEMatrix) {
        this.productBEMatrix = productBEMatrix;
    }

    /**
     * @return the myMapping
     */
    public synchronized AtomAtomMappingContainer getMyMapping() {
        return myMapping;
    }

    /**
     * @param myMapping the myMapping to set
     */
    public synchronized void setMyMapping(AtomAtomMappingContainer myMapping) {
        this.myMapping = myMapping;
    }

    @Override
    public synchronized String toString() {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = getProperty("line.separator");
        result.append("\t");
        for (int i = 0; i < this.getRowDimension() - 1; i++) {
            result.append("\t").append(i);
        }
        result.append(NEW_LINE);
        result.append("\t");
        for (int i = 0; i < this.getRowDimension() - 1; i++) {
            try {
                result.append("\t").append(this.getReactantBEMatrix().getAtom(i).getSymbol())
                        .append(this.getReactantBEMatrix().getAtom(i).getID());
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        result.append(NEW_LINE);
        result.append("\t");
        for (int i = 0; i < this.getRowDimension() - 1; i++) {
            try {
                result.append("\t").append(this.getProductBEMatrix().getAtom(i).getSymbol())
                        .append(this.getProductBEMatrix().getAtom(i).getID());
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        result.append(NEW_LINE);
        for (int i = 0; i < this.getRowDimension() - 1; i++) {
            if (i == this.getRowDimension() - 1) {
                result.append("\t");
            } else {
                try {
                    result.append(this.getReactantBEMatrix().getAtom(i).getSymbol());
                    result.append(this.getReactantBEMatrix().getAtom(i).getID());
                    result.append("\t");
                    result.append(this.getProductBEMatrix().getAtom(i).getSymbol());
                    result.append(this.getProductBEMatrix().getAtom(i).getID());
                    result.append("\t");
                } catch (CDKException ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
            }
            for (int j = 0; j < this.getColumnDimension() - 1; j++) {
                result.append(this.getValue(i, j)).append("\t");
            }
            result.append(NEW_LINE);
        }
        return result.toString();
    }

    private int countAtomOverlap(List<IAtom> atomsE, List<IAtom> atomsP) {

        Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
        Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();
        Map<String, Integer> atomOverlap = new TreeMap<>();

        int leftHandAtomCount = 0;

        leftHandAtomCount = atomsE.stream().filter((a) -> !(a.getSymbol().equals("H"))).map((a) -> {
            if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                atomUniqueCounter1.put(a.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                atomUniqueCounter1.put(a.getSymbol(), counter);
            }
            return a;
        }).map((_item) -> 1).reduce(leftHandAtomCount, Integer::sum);

        int rightHandAtomCount = 0;

        rightHandAtomCount = atomsP.stream().filter((b) -> !(b.getSymbol().equals("H"))).map((b) -> {
            if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                atomUniqueCounter2.put(b.getSymbol(), 1);
            } else {
                int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                atomUniqueCounter2.put(b.getSymbol(), counter);
            }
            return b;
        }).map((_item) -> 1).reduce(rightHandAtomCount, Integer::sum);

        atomUniqueCounter1.keySet().stream().filter((s) -> (atomUniqueCounter2.containsKey(s))).forEach((String s) -> {
            Integer overlap = atomUniqueCounter1.get(s) <= atomUniqueCounter2.get(s)
                    ? atomUniqueCounter1.get(s) : atomUniqueCounter2.get(s);
            atomOverlap.put(s, overlap);
        });
        int total = 0;
        total = atomOverlap.values().stream().map((i) -> i).reduce(total, Integer::sum);
        if (DEBUG) {
            out.println("atomUniqueCounter1 " + leftHandAtomCount);
            out.println("atomUniqueCounter2 " + rightHandAtomCount);
            out.println("overlap " + total);
        }

        return total;
    }

    @Override
    public synchronized Object clone() throws CloneNotSupportedException {
        return super.clone(); //To change body of generated methods, choose Tools | Templates.
    }
}
