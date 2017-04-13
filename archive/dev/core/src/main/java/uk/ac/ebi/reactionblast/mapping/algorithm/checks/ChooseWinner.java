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

 /*
 * Maximize.java
 *
 * Created on February 3, 2006, 12:06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package uk.ac.ebi.reactionblast.mapping.algorithm.checks;

import java.io.IOException;
import java.io.Serializable;
import static java.lang.Double.MIN_VALUE;
import java.util.ArrayList;
import static java.util.Collections.synchronizedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.tools.BondEnergies;
import uk.ac.ebi.reactionblast.mapping.algorithm.Holder;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class ChooseWinner extends Selector implements Serializable {

    private static final long serialVersionUID = 0x296558709L;
    private static final Logger LOG = getLogger(ChooseWinner.class.getName());
    private EBIMatrix stereoMatrix;

    private EBIMatrix similarityMatrix = null;
    private List<Cells> crossMappingTracer = null;
    private Map<Integer, IAtomContainer> educts = null;
    private Map<Integer, IAtomContainer> products = null;
    private double[][] localEnergyMatrix = null;
    private double[][] localCarbonMatrix = null;

    //~--- constructors -------------------------------------------------------
    /**
     * Creates a new instance of Maximize
     *
     * @param eductNameList
     * @param productNameList
     */
    public ChooseWinner(List<String> eductNameList,
            List<String> productNameList) {

        this.rowSize = eductNameList.size();
        this.colSize = productNameList.size();
        this.flagMatrix = new boolean[rowSize][colSize];
        this.localEnergyMatrix = new double[rowSize][colSize];
        this.localCarbonMatrix = new double[rowSize][colSize];
    }

    /**
     * @return the stereoMatrix
     */
    public EBIMatrix getStereoMatrix() {
        return stereoMatrix;
    }

    /**
     * @param stereoMatrix the stereoMatrix to set
     */
    public void setStereoMatrix(EBIMatrix stereoMatrix) {
        this.stereoMatrix = stereoMatrix;
    }

    /**
     * @return the energyMatrix
     */
    public double[][] getEnergyMatrix() {
        return localEnergyMatrix;
    }

    /**
     *
     * @param eductMap
     * @param productMap
     * @param mHolder
     */
    public synchronized void searchWinners(Map<Integer, IAtomContainer> eductMap, Map<Integer, IAtomContainer> productMap, Holder mHolder) {
        initFlagMatrix();
        this.educts = eductMap;
        this.products = productMap;
        this.similarityMatrix = mHolder.getGraphSimilarityMatrix();
        this.setStereoMatrix(mHolder.getStereoMatrix());
        this.crossMappingTracer = synchronizedList(new ArrayList<Cells>());

        boolean isMappingFesiable = checkStatusFlag();
        List<Double> scores = new ArrayList<>();

//        System.out.println("isMappingFesiable " + isMappingFesiable);
        if (isMappingFesiable) {
            calculateEnergy(mHolder.getBestMatchContainer());
            calculateCarbonMatches(mHolder.getBestMatchContainer());
            double maximumSimilarity = 0.0;

            boolean maxValueI;
            boolean maxValueJ;

            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {

                    double similarity = similarityMatrix.getValue(i, j);
                    //System.out.println("similarity "+similarity);
                    //matrix.
                    if (similarity > MIN_VALUE) {
                        maxValueI = isMajorSubgraphRow(similarityMatrix, i, j);
                        maxValueJ = isMajorSubgraphColumn(similarityMatrix, i, j);

//                        System.out.println("similarity " + similarity);
//                        System.out.println("i= " + i + ", j=" + j);
//                        System.out.println(" maxValue I " + maxValueI);
//                        System.out.println(" maxValue J " + maxValueJ);
                        if (maxValueI && maxValueJ) {

                            if (similarity > maximumSimilarity) {
                                maximumSimilarity = similarity;
                                initFlagMatrix();
                                scores.clear();
                            }
                            if (similarity == maximumSimilarity) {
                                this.flagMatrix[i][j] = true;
                                scores.add(similarity);
                            }
                        }
                    }
                }
            }
        }
//        System.out.println("resolveDeadLocks");
        resolveDeadLocks(scores);
//        System.out.println("setWinOverFlags");
        setWinOverFlags();
//        System.out.println("Done");
    }

    /*
    * @return true if a cell in the this.flagMatrix
    * was set true else false
     */
    /**
     *
     * @return
     */
    public synchronized boolean getFlag() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                if (this.flagMatrix[i][j]) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     *
     * @return this.flagMatrix containing the selected cells
     *
     */
    public synchronized boolean[][] getFlagMatrix() {
        return this.flagMatrix;
    }

    private synchronized void initFlagMatrix() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                this.flagMatrix[i][j] = false;
            }
        }
    }

    private synchronized void initFlagMatrix(boolean[][] boolMatrix, int rowSize, int colSize) {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                boolMatrix[i][j] = false;
            }
        }
    }

    private synchronized boolean checkStatusFlag() {
        for (int i = 0; i < rowSize; i++) {
            for (int j = 0; j < colSize; j++) {
                if (similarityMatrix.getValue(i, j) > MIN_VALUE) {
                    return true;
                }
            }
        }
        return false;
    }

    private synchronized void setWinOverFlags() {
        educts.keySet().stream().forEach((Integer indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {
                Cells cell = new Cells();
                cell.indexI = indexI;
                cell.indexJ = indexJ;
                cell.eductName = educts.get(indexI).getID();
                cell.productName = products.get(indexJ).getID();
                if (ChooseWinner.this.flagMatrix[indexI][indexJ]) {
                    if (cell.eductName.equalsIgnoreCase(cell.productName)) {
                    } else {
                        ChooseWinner.this.flagMatrix[indexI][indexJ] = !checkTwinMapping(cell);
                    }
                }
            });
        });
    }

    private synchronized boolean checkTwinMapping(Cells refCell) {
        boolean _statusFlag = false;
        for (Cells cell : crossMappingTracer) {

            //System.out.println("Visitor AtomContainer " + moleculeName + " i:" + _value.get(1) + " j:" + j);
            if (cell.eductName.equals(refCell.eductName)
                    && cell.productName.equals(refCell.productName)) {
                if (cell.indexI == refCell.indexI
                        || cell.indexJ == refCell.indexJ) {
                    _statusFlag = true;
                    break;
                }
            }
        }

        if (!_statusFlag) {
            crossMappingTracer.add(refCell);
        }
        return _statusFlag;
    }

    private synchronized void resolveDeadLocks(List<Double> scores) {
        boolean[][] deadlockFreeFlagMatrix = new boolean[rowSize][colSize];
        initFlagMatrix(deadlockFreeFlagMatrix, rowSize, colSize);
        scores.stream().map((score) -> new DeadLockResolver().resolver(score)).forEach((choosenCell) -> {
            deadlockFreeFlagMatrix[choosenCell.indexI][choosenCell.indexJ] = true;
        });
        this.flagMatrix = deadlockFreeFlagMatrix;
    }

    private void calculateEnergy(BestMatch bestMatchContainer) {
        //System.out.println("Calculating energy ");
        educts.keySet().stream().forEach((Integer indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {
                localEnergyMatrix[indexI][indexJ] = -1;
            });
        });
        educts.keySet().stream().forEach((Integer indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {
                try {
                    AtomAtomMapping atomMatch = bestMatchContainer.getAtomMatch(indexI, indexJ);
                    if (!atomMatch.isEmpty()) {
                        localEnergyMatrix[indexI][indexJ] = getMappedMoleculeEnergies(atomMatch);
//                        System.out.println(" localEnergyMatrix[indexI][indexJ]  " + localEnergyMatrix[indexI][indexJ]);
                    }
                } catch (CDKException | IOException ex) {
                    Logger.getLogger(ChooseWinner.class.getName()).log(Level.SEVERE, null, ex);
                }
            });
        });
    }

    private void calculateCarbonMatches(BestMatch bestMatchContainer) {
        //System.out.println("Calculating carbon ");
        educts.keySet().stream().forEach((Integer indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {
                localCarbonMatrix[indexI][indexJ] = -1;
            });
        });
        educts.keySet().stream().forEach((Integer indexI) -> {
            products.keySet().stream().forEach((indexJ) -> {

                AtomAtomMapping atomMatch;
                try {
                    atomMatch = bestMatchContainer.getAtomMatch(indexI, indexJ);
                    if (!atomMatch.isEmpty()) {
                        localCarbonMatrix[indexI][indexJ] = getMappedCarbonsInMolecules(atomMatch);
//                        System.out.println(" localCarbonMatrix[indexI][indexJ]  " + localCarbonMatrix[indexI][indexJ]);
                    }
                } catch (IOException ex) {
                    Logger.getLogger(ChooseWinner.class.getName()).log(Level.SEVERE, null, ex);
                }
            });
        });
    }

    private double getMappedCarbonsInMolecules(AtomAtomMapping atomMatch) {
        double carbonScore = 0.0;
        for (IAtom a : atomMatch.getQuery().atoms()) {
            for (IAtom b : atomMatch.getQuery().atoms()) {
                if (a != b && atomMatch.getMappingsByAtoms().containsKey(a)
                        && atomMatch.getMappingsByAtoms().containsKey(b)) {
                    IBond bond = atomMatch.getQuery().getBond(a, b);
                    if (bond != null && a.getSymbol().equals("C") && b.getSymbol().equals("C")) {
                        carbonScore += 1;
                    }
                }
            }
        }
        return carbonScore;
    }

    /**
     * Chosen cell of the matrix is stored here
     */
    class Cells {

        String eductName;
        String productName;
        int indexI;
        int indexJ;

    }

    /**
     * Resolves deadlocks if more than one cell clashes with same scores. The
     * decision is then made on the max. stereo score and min. energy score
     */
    class DeadLockResolver {

        private synchronized double getMaxCarbon(List<ChooseWinner.Cells> choosenCells) {
            double max = -999;
            for (ChooseWinner.Cells cell : choosenCells) {
                double val = localCarbonMatrix[cell.indexI][cell.indexJ];
                if (val > max) {
                    max = val;
                }
            }
            return max;
        }

        private synchronized double getMaxStereo(List<ChooseWinner.Cells> choosenCells) {
            double max = -999;
            for (ChooseWinner.Cells cell : choosenCells) {
                double val = getStereoMatrix().getValue(cell.indexI, cell.indexJ);
                if (val > max) {
                    max = val;
                }
            }
            return max;
        }

        private synchronized double getMinEnergy(List<ChooseWinner.Cells> choosenCells) {
            double min = 999999;
            for (ChooseWinner.Cells cell : choosenCells) {
                double val = getEnergyMatrix()[cell.indexI][cell.indexJ];
                if (val > -1 && val < min) {
                    min = val;
                }
            }
            return min;
        }

        /**
         * Returns best match for the given score
         *
         * @param choosenScore
         * @return
         */
        public synchronized ChooseWinner.Cells resolver(double choosenScore) {
            List<ChooseWinner.Cells> choosenCells = new ArrayList<>();
            for (int i = 0; i < rowSize; i++) {
                for (int j = 0; j < colSize; j++) {
                    if (flagMatrix[i][j] && choosenScore > MIN_VALUE
                            && similarityMatrix.getValue(i, j) == choosenScore) {
                        ChooseWinner.Cells cells = new ChooseWinner.Cells();
                        cells.indexI = i;
                        cells.indexJ = j;
                        choosenCells.add(cells);
                    }
                }
            }
            double maxStereo = getMaxStereo(choosenCells);
            double minEnergy = getMinEnergy(choosenCells);
            double maxCarbon = getMaxCarbon(choosenCells);
//            System.out.println("------------------------------");
//            System.out.println("maxStereo " + maxStereo);
//            System.out.println("minEnergy " + minEnergy);
//            System.out.println("maxCarbon " + maxCarbon);
//            System.out.println("------------------------------");

            ChooseWinner.Cells chosenCell = choosenCells.listIterator().next();

            /*
            * Filter by Condition 1
             */
            for (ChooseWinner.Cells cell : choosenCells) {
                double stereoVal = stereoMatrix.getValue(cell.indexI, cell.indexJ);
                double energyVal = getEnergyMatrix()[cell.indexI][cell.indexJ];
                double carbonVal = localCarbonMatrix[cell.indexI][cell.indexJ];
//                System.out.println("stereoVal " + stereoVal);
//                System.out.println("energyVal " + energyVal);
//                System.out.println("carbonVal " + carbonVal);
//                System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%");

                if (stereoVal >= maxStereo && energyVal <= minEnergy && carbonVal >= maxCarbon) {
//                    System.out.println("Cond 1");
                    chosenCell = cell;
                    return chosenCell;
                }
            }
            /*
             * Filter by Condition 2
             */
            for (ChooseWinner.Cells cell : choosenCells) {
                double stereoVal = stereoMatrix.getValue(cell.indexI, cell.indexJ);
                double energyVal = getEnergyMatrix()[cell.indexI][cell.indexJ];
//                System.out.println("stereoVal " + stereoVal);
//                System.out.println("energyVal " + energyVal);

                if (stereoVal >= maxStereo && energyVal <= minEnergy) {
//                    System.out.println("Cond 2");
                    chosenCell = cell;
                    return chosenCell;
                }
            }
            /*
             * Filter by Condition 3
             */
            for (ChooseWinner.Cells cell : choosenCells) {
                double energyVal = getEnergyMatrix()[cell.indexI][cell.indexJ];
//                System.out.println("energyVal " + energyVal);

                if (energyVal <= minEnergy) {
//                    System.out.println("Cond 3");
                    chosenCell = cell;
                    return chosenCell;
                }
            }

            return chosenCell;
        }
    }

    private synchronized Double getMappedMoleculeEnergies(AtomAtomMapping mcsAtomSolution) throws CDKException {

//        System.out.println("\nSort By Energies");
        double totalBondEnergy = -9999.0;

        IAtomContainer educt = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, mcsAtomSolution.getQuery());
        IAtomContainer product = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, mcsAtomSolution.getTarget());

        for (int i = 0; i < educt.getAtomCount(); i++) {
            educt.getAtom(i).setFlag(999, false);
        }

        for (int i = 0; i < product.getAtomCount(); i++) {
            product.getAtom(i).setFlag(999, false);
        }

        Map<IAtom, IAtom> mappingsByAtoms = mcsAtomSolution.getMappingsByAtoms();
        mappingsByAtoms.entrySet().stream().map((mapping) -> {
            mapping.getKey().setFlag(999, true);
            return mapping;
        }).forEach((mapping) -> {
            mapping.getValue().setFlag(999, true);
        });
        totalBondEnergy = getEnergy(educt, product);


        /*
         * Reset the flag
         */
        for (int i = 0; i < educt.getAtomCount(); i++) {
            educt.getAtom(i).setFlag(999, false);
        }

        for (int i = 0; i < product.getAtomCount(); i++) {
            product.getAtom(i).setFlag(999, false);
        }

        return totalBondEnergy;
    }

    private synchronized static double getEnergy(IAtomContainer educt, IAtomContainer product) throws CDKException {
        Double eEnergy = 0.0;
        BondEnergies bondEnergy = BondEnergies.getInstance();
        for (int i = 0; i < educt.getBondCount(); i++) {
            IBond bond = educt.getBond(i);
            eEnergy += getBondEnergy(bond, bondEnergy);
        }
        Double pEnergy = 0.0;
        for (int j = 0; j < product.getBondCount(); j++) {
            IBond bond = product.getBond(j);
            pEnergy += getBondEnergy(bond, bondEnergy);
        }
        return (eEnergy + pEnergy);
    }

    private synchronized static double getBondEnergy(IBond bond, BondEnergies bondEnergy) {
        double energy = 0.0;
        if ((bond.getAtom(0).getFlag(999) == true && bond.getAtom(1).getFlag(999) == false)
                || (bond.getAtom(0).getFlag(999) == false && bond.getAtom(1).getFlag(999) == true)) {
            int val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
            energy = val;
        }
        return energy;
    }
}
