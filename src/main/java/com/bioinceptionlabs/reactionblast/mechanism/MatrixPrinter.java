/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.reactionblast.mechanism.BEMatrix;
import com.bioinceptionlabs.reactionblast.mechanism.RMatrix;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class MatrixPrinter extends Object {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MatrixPrinter.class);


    /**
     * This method prints the matrix to the standard output
     *
     * @param rMatrix R-Matrix to be Printed
     */
    public static void printReactionMatrix(RMatrix rMatrix) {
        try {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                sb.append("\t\t").append(i);
            }
            sb.append(System.lineSeparator());
            for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                sb.append("\t\t").append(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol())
                        .append(rMatrix.getReactantBEMatrix().getAtom(i).getID());
            }
            sb.append(System.lineSeparator());
            for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                sb.append("\t\t").append(rMatrix.getProductBEMatrix().getAtom(i).getSymbol())
                        .append(rMatrix.getProductBEMatrix().getAtom(i).getID());
            }
            sb.append(System.lineSeparator());
            for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                if (i == rMatrix.getRowDimension() - 1) {
                    sb.append("\t");
                } else {
                    sb.append(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol())
                            .append(rMatrix.getReactantBEMatrix().getAtom(i).getID())
                            .append("\t").append(rMatrix.getProductBEMatrix().getAtom(i).getSymbol())
                            .append(rMatrix.getProductBEMatrix().getAtom(i).getID());
                }
                for (int j = 0; j < rMatrix.getColumnDimension(); j++) {
                    sb.append("\t").append(rMatrix.getValue(i, j));
                }
                sb.append(System.lineSeparator());
            }
            LOGGER.debug(sb.toString());
        } catch (CDKException ex) {
            LOGGER.debug("A CDKException has been arisen while printing the RMatrix");
        }
    }

    /**
     *
     * @param outputFile
     * @param rMatrix R-Matrix
     * @throws IOException
     */
    public static void writeReactionMatrix(File outputFile, RMatrix rMatrix) throws IOException {
        try (BufferedWriter matrixFileWriter = new BufferedWriter(new FileWriter(outputFile))) {
            matrixFileWriter.newLine();
            try {
                for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                    matrixFileWriter.write("\t" + i);
                }
                matrixFileWriter.newLine();
                for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                    matrixFileWriter.write("\t" + rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                            + rMatrix.getReactantBEMatrix().getAtom(i).getID());
                }
                matrixFileWriter.newLine();
                for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                    matrixFileWriter.write("\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                            + rMatrix.getProductBEMatrix().getAtom(i).getID());
                }
                matrixFileWriter.newLine();
                for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                    if (i == rMatrix.getRowDimension() - 1) {
                        matrixFileWriter.write("\t");
                    } else {
                        matrixFileWriter.write(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                                + rMatrix.getReactantBEMatrix().getAtom(i).getID()
                                + "\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                                + rMatrix.getProductBEMatrix().getAtom(i).getID());
                    }
                    for (int j = 0; j < rMatrix.getColumnDimension(); j++) {
                        matrixFileWriter.write("\t" + rMatrix.getValue(i, j));
                    }
                    matrixFileWriter.newLine();
                }

            } catch (CDKException ex) {
                LOGGER.debug("A CDKException has been arisen while printing the RMatrix");
            }
        }
    }

    /**
     * This method prints the matrix to the standard output
     *
     * @param beMatrix
     *
     */
    public static void printBEMatrix(BEMatrix beMatrix) {
        List<IAtom> atomArray = beMatrix.getAtoms();
        StringBuilder sb = new StringBuilder();
        sb.append(atomArray.size()).append(System.lineSeparator());
        for (int i = 0; i < atomArray.size(); i++) {
            sb.append(atomArray.get(i).getSymbol()).append(atomArray.get(i).getID()).append("\t");
        }
        sb.append(System.lineSeparator());
        for (int i = 0; i < beMatrix.getRowDimension(); i++) {
            for (int j = 0; j < beMatrix.getColumnDimension(); j++) {
                sb.append(beMatrix.getValue(i, j)).append("\t");
            }
            sb.append(System.lineSeparator());
        }
        LOGGER.debug(sb.toString());
    }

    /**
     *
     * @param outputFile
     * @param beMatrix BE-Matrix
     * @throws IOException
     */
    public static void writeBEMatrix(File outputFile, BEMatrix beMatrix) throws IOException {
        List<IAtom> atomArray = beMatrix.getAtoms();
        try (BufferedWriter matrixFileWriter = new BufferedWriter(new FileWriter(outputFile))) {
            matrixFileWriter.newLine();
            matrixFileWriter.write(atomArray.size());
            matrixFileWriter.newLine();
            for (int i = 0; i < atomArray.size(); i++) {
                matrixFileWriter.write(atomArray.get(i).getSymbol() + atomArray.get(i).getID() + "\t");
            }
            matrixFileWriter.newLine();
            for (int i = 0; i < beMatrix.getRowDimension(); i++) {
                for (int j = 0; j < beMatrix.getColumnDimension(); j++) {
                    matrixFileWriter.write(beMatrix.getValue(i, j) + "\t");
                }
                matrixFileWriter.newLine();
            }
        }
    }

    MatrixPrinter() {
    }
}
