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
package uk.ac.ebi.reactionblast.mechanism.helper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.System.out;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.reactionblast.mechanism.BEMatrix;
import uk.ac.ebi.reactionblast.mechanism.RMatrix;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
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
            for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                out.print("\t\t" + i);
            }
            out.println();
            for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                out.print("\t\t" + rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                        + rMatrix.getReactantBEMatrix().getAtom(i).getID());
            }
            out.println();
            for (int i = 0; i < rMatrix.getRowDimension() - 1; i++) {
                out.print("\t\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                        + rMatrix.getProductBEMatrix().getAtom(i).getID());
            }
            out.println();
            for (int i = 0; i < rMatrix.getRowDimension(); i++) {
                if (i == rMatrix.getRowDimension() - 1) {
                    out.print("\t");
                } else {
                    out.print(rMatrix.getReactantBEMatrix().getAtom(i).getSymbol()
                            + rMatrix.getReactantBEMatrix().getAtom(i).getID()
                            + "\t" + rMatrix.getProductBEMatrix().getAtom(i).getSymbol()
                            + rMatrix.getProductBEMatrix().getAtom(i).getID());
                }
                for (int j = 0; j < rMatrix.getColumnDimension(); j++) {
                    out.print("\t" + rMatrix.getValue(i, j));
                }
                out.println();
            }
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
        FileWriter matrixFile = new FileWriter(outputFile);
        try (BufferedWriter matrixFileWriter = new BufferedWriter(matrixFile)) {
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
                matrixFileWriter.close();

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
        out.println(atomArray.size());
        for (int i = 0; i < atomArray.size(); i++) {
            out.print(atomArray.get(i).getSymbol() + atomArray.get(i).getID() + "\t");
        }
        out.println();
        for (int i = 0; i < beMatrix.getRowDimension(); i++) {
            for (int j = 0; j < beMatrix.getColumnDimension(); j++) {
                out.print(beMatrix.getValue(i, j) + "\t");
            }
            out.println();
        }
    }

    /**
     *
     * @param outputFile
     * @param beMatrix BE-Matrix
     * @throws IOException
     */
    public static void writeBEMatrix(File outputFile, BEMatrix beMatrix) throws IOException {
        List<IAtom> atomArray = beMatrix.getAtoms();
        FileWriter matrixFile = new FileWriter(outputFile);
        try (BufferedWriter matrixFileWriter = new BufferedWriter(matrixFile)) {
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
            matrixFileWriter.close();
        }
    }

    MatrixPrinter() {
    }
}
