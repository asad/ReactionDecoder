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
package com.bioinceptionlabs.reactionblast.tools;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import static java.io.StreamTokenizer.TT_EOF;
import static java.io.StreamTokenizer.TT_EOL;
import static java.io.StreamTokenizer.TT_WORD;
import static java.lang.Double.valueOf;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.System.arraycopy;
import static java.lang.System.getProperty;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import static java.util.Locale.US;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static java.lang.Math.hypot;
import static java.lang.Math.pow;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Jama = Java EBIMatrix class.
 * <P>
 * The Java EBIMatrix datalass provides the fundamental operations of numerical
 * linear algebra. Various constructors create Matrices from two dimensional
 * arrays of double precision doubleing point numbers. Various "gets" and "sets"
 * provide access to submatrices and matrix elements. Several methods implement
 * basic matrix arithmetic, including matrix addition and multiplication, matrix
 * norms, and element-by-element array operations. Methods for reading and
 * printing matrices are also included. All the operations in this version of
 * the EBIMatrix datalass involve real matrices. dataomplex matrices may be
 * handled in a future version.
 * <P>
 * Five fundamental matrix decompositions, which consist of pairs or triples of
 * matrices, permutation vectors, and the like, produce results in five
 * decomposition classes. These decompositions are accessed by the EBIMatrix
 * class to compute solutions of simultaneous linear equations, determinants,
 * inverses and other matrix functions. The five decompositions are:
 * <P>
 * <UL>
 * <LI>dataholesky Decomposition of symmetric, positive definite matrices.
 * <LI>LU Decomposition of rectangular matrices.
 * <LI>QR Decomposition of rectangular matrices. <LI>Singular Value
 * Decomposition of rectangular matrices.
 * <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square
 * matrices. </UL> <DL> <DT><B>Example of use:</B></DT>
 * <P>
 * <DD>Solve a linear system matrix x = b and compute the residual norm, ||b -
 * matrix x||.
 * <P>
 * <
 * PRE>
 * double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}}; EBIMatrix matrix = new
 * EBIMatrix(vals); EBIMatrix b = EBIMatrix.random(3,1); EBIMatrix x =
 * matrix.solve(b); EBIMatrix r = matrix.times(x).minus(b); double rnorm =
 * r.normInf();
 * </PRE></DD> </DL>
 *
 * @author The MathWorks, Inc. and the National Institute of Standards and
 * Technology.
 * @version 5 August 1998
 */
public class EBIMatrix extends Object implements Cloneable, java.io.Serializable {

    static final String NEW_LINE = getProperty("line.separator");

    private static final long serialVersionUID = 19787786981017786L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(EBIMatrix.class);

    /**
     * Solves a linear equation system with Gauss elimination.
     *
     * @param matrix
     * @param vector
     * @return
     * @keyword Gauss elimination
     */
    public static List<Double> elimination(EBIMatrix matrix, List<Double> vector) {
        int i, j, k, ipvt;
        int n = vector.size();
        int[] pivot = new int[n];
        double c, temp;
        //double[] x = new double[n];
        EBIMatrix a = matrix.duplicate();
        List<Double> b = new ArrayList<>(vector);
        for (j = 0; j < (n - 1); j++) {
            c = abs(a.matrix[j][j]);
            pivot[j] = j;
            ipvt = j;
            for (i = j + 1; i < n; i++) {
                if (abs(a.matrix[i][j]) > c) {
                    c = abs(a.matrix[i][j]);
                    ipvt = i;
                }
            }

            // Exchanges rows when necessary
            if (pivot[j] != ipvt) {
                pivot[j] = ipvt;
                pivot[ipvt] = j;
                for (k = 0; k < n; k++) {
                    temp = a.matrix[j][k];
                    a.matrix[j][k] = a.matrix[pivot[j]][k];
                    a.matrix[pivot[j]][k] = temp;
                }

                temp = b.get(j);
                b.set(j, b.get(pivot[j]));
                b.set(pivot[j], temp);
            }

            // Store multipliers
            for (i = j + 1; i < n; i++) {
                a.matrix[i][j] /= a.matrix[j][j];
            }

            // Give elements below the diagonal a zero value
            for (i = j + 1; i < n; i++) {
                for (k = j + 1; k < n; k++) {
                    a.matrix[i][k] -= a.matrix[i][j] * a.matrix[j][k];
                }
                b.set(i, b.get(i) - a.matrix[i][j] * b.get(j));

                a.matrix[i][j] = 0.0; // Not necessary
            }
        }
        // Rueckwaertseinsetzen (which is?)
        List<Double> result = new ArrayList<>(n);
        result.set(n - 1, b.get(n - 1) / a.matrix[n - 1][n - 1]);
        for (j = n - 2; j >= 0; j--) {
            result.set(j, b.get(j));
            for (k = n - 1; k > j; k--) {
                result.set(j, result.get(j) - result.get(k) * a.matrix[j][k]);
            }

            result.set(j, result.get(j) / a.matrix[j][k]);
        }
        return result;
    }

    /**
     * Generate matrix with random elements
     *
     * @param m
     * @param n
     * @return An rows-by-columns matrix with uniformly distributed random
     * elements.
     */
    public static EBIMatrix random(int m, int n) {
        EBIMatrix A = new EBIMatrix(m, n);
        double[][] X = A.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = Math.random();
            }
        }
        return A;
    }

    /**
     * Generate identity matrix
     *
     * @param m
     * @param n
     * @return An rows-by-columns matrix with ones on the diagonal and zeros
     * elsewhere.
     */
    public static EBIMatrix identity(int m, int n) {
        EBIMatrix A = new EBIMatrix(m, n);
        double[][] X = A.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = (i == j ? 1.0 : 0.0);
            }
        }
        return A;
    }

    /**
     * Read a matrix from a stream. The format is the same the print method, so
     * printed matrices can be read back in (provided they were printed using US
     * Locale). Elements are separated by whitespace, all the elements for each
     * row appear on a single line, the last row is followed by a blank line.
     *
     * @param input the input stream.
     * @return
     * @throws java.io.IOException
     */
    @SuppressWarnings(value = {"empty-statement", "unchecked"})
    public static EBIMatrix read(BufferedReader input) throws java.io.IOException {
        StreamTokenizer tokenizer = new StreamTokenizer(input);

        // Although StreamTokenizer will parse numbers, it doesn't recognize
        // scientific notation (E or D); however, double.valueOf does.
        // The strategy here is to disable StreamTokenizer's number parsing.
        // We'll only get whitespace delimited words, EOL's and EOF's.
        // These words should all be numbers, for double.valueOf to parse.
        tokenizer.resetSyntax();
        tokenizer.wordChars(0, 255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);
        ArrayList v = new ArrayList();

        // Ignore initial empty lines
        while (tokenizer.nextToken() == TT_EOL) {
            ;
        }
        if (tokenizer.ttype == TT_EOF) {
            throw new java.io.IOException("Unexpected EOF on matrix read.");
        }
        do {
            v.add(valueOf(tokenizer.sval)); // Read & store 1st row.
        } while (tokenizer.nextToken() == TT_WORD);

        int n = v.size();  // Now we've got the number of columns!
        double row[] = new double[n];
        for (int j = 0; j < n; j++) {
            // extract the elements of the 1st row.
            row[j] = ((Number) v.get(j)).doubleValue();
        }
        v.removeAll(v);
        v.add(row);  // Start storing rows instead of columns.
        while (tokenizer.nextToken() == TT_WORD) {
            // While non-empty lines
            v.add(row = new double[n]);
            int j = 0;
            do {
                if (j >= n) {
                    throw new java.io.IOException("Row " + v.size() + " is too long.");
                }
                valueOf(tokenizer.sval);
            } while (tokenizer.nextToken() == TT_WORD);
            if (j < n) {
                throw new java.io.IOException("Row " + v.size() + " is too short.");
            }
        }
        int m = v.size();  // Now we've got the number of rows.
        double[][] A = new double[m][];
        for (int i = 0; i < v.size(); i++) {
            A[i] = (double[]) v.get(i);// duplicate the rows out of the vector
        }
        return new EBIMatrix(A);
    }

    /*
     * ------------------------ Public Methods ------------------------
     */
    /**
     *
     * @param A
     * @return
     */
    public static EBIMatrix constructWithCopy(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        EBIMatrix X = new EBIMatrix(m, n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
            arraycopy(A[i], 0, C[i], 0, n);
        }
        return X;
    }
    /*
     * ------------------------ Class variables ------------------------
     */
    /**
     * Array for internal storage of elements.
     *
     * @serial internal array storage.
     */
    private double[][] matrix;
    /**
     * Row and column dimensions.
     *
     * @serial row dimension.
     * @serial column dimension.
     */
    private int rows, columns;

    /*
     * ------------------------ Constructors ------------------------
     */
    /**
     * creates a new EBIMatrix.
     *
     * @param rows
     * @param columns
     */
    public EBIMatrix(int rows, int columns) {
        this.rows = rows;
        this.columns = columns;
        matrix = new double[rows][columns];
    }

    /**
     * data construct a matrix from a 2-D array.
     *
     * @param A
     * @see #constructWithCopy
     */
    public EBIMatrix(double[][] A) {
        rows = A.length;
        columns = A[0].length;
        for (int i = 0; i < rows; i++) {
            if (A[i].length != columns) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
        }
        matrix = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            arraycopy(A[i], 0, matrix[i], 0, columns);
        }
    }

    /**
     * data construct a matrix quickly without checking arguments.
     *
     * @param A
     * @param m
     * @param n
     */
    public EBIMatrix(double[][] A, int m, int n) {
        this.matrix = A;
        this.rows = m;
        this.columns = n;
    }

    /*
     * ------------------------ Constructors ------------------------
     */
    /**
     * dataonstruct an rows-by-columns constant matrix.
     *
     * @param m
     * @param n
     * @param s Fill the matrix with this scalar value.
     */
    public EBIMatrix(int m, int n, double s) {
        this.rows = m;
        this.columns = n;
        matrix = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = s;
            }
        }
    }

    /**
     * data construct a matrix from a one-dimensional packed array
     *
     * @param vals One-dimensional array of doubles, packed by columns (ala
     * Fortran).
     * @param m
     */
    public EBIMatrix(double vals[], int m) {
        this.rows = m;
        this.columns = (m != 0 ? vals.length / m : 0);
        if (m * columns != vals.length) {
            throw new IllegalArgumentException("Array length must be a multiple of m.");
        }
        matrix = new double[m][columns];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = vals[i + j * m];
            }
        }
    }

    /**
     *
     * @param v default value for the Matrix cells
     */
    public void initMatrix(double v) {
        for (int i = 0; i < rows; i++) {
            java.util.Arrays.fill(matrix[i], v);
        }
    }

    /**
     * Get a single element.
     *
     * @param i Row index.
     * @param j dataolumn index.
     * @return matrix(i,j)
     * @exception ArrayIndexOutOfBoundsException
     */
    public double getValue(int i, int j) {

        double val = -1.0d;
        if (i >= 0 && i < rows && j >= 0 && j < columns) {
            val = matrix[i][j];
        } else {
            LOGGER.debug("Error: Array out of bounds [" + i + "," + j + "] for [" + rows + "," + columns + "]");
        }
        return val;
    }

    /**
     * Get a single element without bounds checking.
     * Use this on hot paths where indices are already known to be valid,
     * to avoid the overhead of bounds checking and logging in {@link #getValue}.
     *
     * @param i Row index.
     * @param j Column index.
     * @return matrix(i,j)
     */
    public double getValueUnsafe(int i, int j) {
        return matrix[i][j];
    }

    /**
     * Make a deep duplicate of a matrix
     *
     * @return
     */
    public EBIMatrix duplicate() {
        EBIMatrix result = new EBIMatrix(rows, columns);
        double[][] data = result.getArray();
        for (int i = 0; i < rows; i++) {
            arraycopy(matrix[i], 0, data[i], 0, columns);
        }
        return result;
    }

    /**
     * data clone the EBIMatrix object.
     *
     * @return
     * @throws java.lang.CloneNotSupportedException
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        return this.duplicate();
    }

    /**
     * data duplicate the internal two-dimensional array.
     *
     * @return Two-dimensional array duplicate of matrix elements.
     */
    public double[][] getArrayCopy() {
        double[][] C = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            arraycopy(matrix[i], 0, C[i], 0, columns);
        }
        return C;
    }

    /**
     * Make a one-dimensional column packed duplicate of the internal array.
     *
     * @return EBIMatrix elements packed in a one-dimensional array by columns.
     */
    public double[] getColumnPackedCopy() {
        double[] vals = new double[rows * columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                vals[i + j * rows] = matrix[i][j];
            }
        }
        return vals;
    }

    /**
     * Make a one-dimensional row packed duplicate of the internal array.
     *
     * @return EBIMatrix elements packed in a one-dimensional array by rows.
     */
    public double[] getRowPackedCopy() {
        double[] vals = new double[rows * columns];
        for (int i = 0; i < rows; i++) {
            arraycopy(matrix[i], 0, vals, i * columns, columns);
        }
        return vals;
    }

    /**
     * Access the internal two-dimensional array.
     *
     * @return Pointer to the two-dimensional array of matrix elements.
     */
    public double[][] getArray() {
        return matrix;
    }

    /**
     * Get a sub-matrix.
     *
     * @param rowStart Initial row index
     * @param rowEnd Final row index
     * @param colStart Initial column index
     * @param colEnd Final column index
     * @return matrix(rowStart:row,colStart:colEnd)
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public EBIMatrix getMatrix(int rowStart, int rowEnd, int colStart, int colEnd) {
        EBIMatrix X = new EBIMatrix(rowEnd - rowStart + 1, colEnd - colStart + 1);
        double[][] B = X.getArray();
        try {
            for (int i = rowStart; i <= rowEnd; i++) {
                for (int j = colStart; j <= colEnd; j++) {
                    B[i - rowStart][j - colStart] = matrix[i][j];
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /**
     * EBIMatrix transpose.
     *
     * @return matrix'
     */
    public EBIMatrix transpose() {
        EBIMatrix X = new EBIMatrix(columns, rows);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[j][i] = matrix[i][j];
            }
        }
        return X;
    }

    /**
     *
     * @param row
     * @param col
     * @param value
     * @return
     */
    public boolean setValue(int row, int col, double value) {

        double val = value;
        boolean flag = false;

        if (row >= 0 && row < rows && col >= 0 && col < columns) {
            flag = true;
            matrix[row][col] = val;
        } else {
            try {

                throw new CDKException("Array out of Bound");
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        return flag;
    }

    /**
     * Get row dimension.
     *
     * @return rows, the number of rows.
     */
    public int getRowDimension() {
        return this.rows;
    }

    /**
     * Get column dimension.
     *
     * @return columns, the number of columns.
     */
    public int getColumnDimension() {
        return this.columns;
    }

    /**
     *
     * @return
     */
    public List<Double> getDiagonalElements() {

        List<Double> val = new ArrayList<>();

        if (rows == columns) {

            for (int i = 0; i < rows; i++) {

                for (int j = 0; j < columns; j++) {

                    if (i == j) {

                        val.add(matrix[i][j]);

                    }

                }

            }
        } else {

            LOGGER.debug("Row =/= Columns");
        }

        return val;

    }

    /**
     * Get a submatrix.
     *
     * @param r Array of row indices.
     * @param c Array of column indices.
     * @return matrix(r(:),c(:))
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public EBIMatrix getMatrix(int[] r, int[] c) {
        EBIMatrix X = new EBIMatrix(r.length, c.length);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i][j] = matrix[r[i]][c[j]];
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /**
     * Get a submatrix.
     *
     * @param rowStart Initial row index
     * @param rowEnd Final row index
     * @param c Array of column indices.
     * @return matrix(rowStart:row,c(:))
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public EBIMatrix getMatrix(int rowStart, int rowEnd, int[] c) {
        EBIMatrix X = new EBIMatrix(rowEnd - rowStart + 1, c.length);
        double[][] B = X.getArray();
        try {
            for (int i = rowStart; i <= rowEnd; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i - rowStart][j] = matrix[i][c[j]];
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /**
     * Get a submatrix.
     *
     * @param r Array of row indices.
     * @param colStart
     * @param colEnd
     * @return matrix(r(:),colStart:colEnd)
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public EBIMatrix getMatrix(int[] r, int colStart, int colEnd) {
        EBIMatrix X = new EBIMatrix(r.length, colEnd - colStart + 1);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = colStart; j <= colEnd; j++) {
                    B[i][j - colStart] = matrix[r[i]][j];
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /**
     * Set a single element.
     *
     * @param i Row index.
     * @param j dataolumn index.
     * @param s matrix(i,j).
     * @exception ArrayIndexOutOfBoundsException
     */
    public void set(int i, int j, double s) {
        matrix[i][j] = s;
    }

    /**
     * Set a submatrix.
     *
     * @param rowStart Initial row index
     * @param colStart Initial column index
     * @param colEnd Final column index
     * @param rowEnd Final row index
     * @param X EBIMatrix
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     *
     */
    public void setMatrix(int rowStart, int rowEnd, int colStart, int colEnd, EBIMatrix X) {
        try {
            for (int i = rowStart; i <= rowEnd; i++) {
                for (int j = colStart; j <= colEnd; j++) {
                    matrix[i][j] = X.getValue(i - rowStart, j - colStart);
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }

    /**
     * Set a submatrix.
     *
     * @param r Array of row indices.
     * @param c Array of column indices.
     * @param X matrix(r(:),c(:))
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public void setMatrix(int[] r, int[] c, EBIMatrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    matrix[r[i]][c[j]] = X.getValue(i, j);
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices" + e);
        }
    }

    /**
     * Set a submatrix.
     *
     * @param r Array of row indices.
     * @param colStart Initial column index
     * @param colEnd Final column index
     * @param X matrix(r(:),colStart:colEnd)
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public void setMatrix(int[] r, int colStart, int colEnd, EBIMatrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = colStart; j <= colEnd; j++) {
                    matrix[r[i]][j] = X.getValue(i, j - colStart);
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }

    /**
     * Set a submatrix.
     *
     * @param rowStart Initial row index
     * @param rowEnd Final row index
     * @param c Array of column indices.
     * @param X matrix(rowStart:row,c(:))
     * @exception ArrayIndexOutOfBoundsException Submatrix indices
     */
    public void setMatrix(int rowStart, int rowEnd, int[] c, EBIMatrix X) {
        try {
            for (int i = rowStart; i <= rowEnd; i++) {
                for (int j = 0; j < c.length; j++) {
                    matrix[i][c[j]] = X.getValue(i - rowStart, j);
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }

    /**
     *
     * @param iPos index
     * @param jPos index
     * @return
     */
    public boolean is_element_max_in_column(int iPos, int jPos) {

        boolean flag = true;
        double refValue = matrix[iPos][jPos];

        for (int j = 0; j < columns; j++) {
            if (j != jPos) {
                double searchValue = matrix[iPos][j];

                if (searchValue > refValue) {
                    flag = false;

                }
            }
        }

        return flag;
    }

    /**
     *
     * @param iPos index
     * @param jPos index
     * @return
     */
    public boolean is_element_min_in_column(int iPos, int jPos) {

        boolean flag = true;
        double refValue = matrix[iPos][jPos];

        for (int j = 0; j < columns; j++) {
            if (j != jPos) {
                double searchValue = matrix[iPos][j];

                if (searchValue < refValue) {
                    flag = false;

                }
            }
        }

        return flag;
    }

    /**
     *
     * @param RowSize Size of the new Matrix Row
     * @param ColSize Size of the new Matrix dataoloumn
     */
    public void reSizeMatrix(int RowSize, int ColSize) {

        this.rows = RowSize;
        this.columns = ColSize;
        matrix = new double[rows][columns];
    }

    /**
     *
     * @param iPos index
     * @param jPos index
     * @return
     */
    public boolean is_element_max_in_row(int iPos, int jPos) {

        boolean flag = true;

        double refValue = matrix[iPos][jPos];

        for (int i = 0; i < rows; i++) {
            if (i != iPos) {
                double searchValue = matrix[i][jPos];

                if (searchValue > refValue) {
                    flag = false;

                }
            }
        }

        return flag;
    }

    /**
     *
     * @param iPos index
     * @param jPos index
     * @return
     */
    public boolean is_element_min_in_row(int iPos, int jPos) {

        boolean flag = true;

        double refValue = matrix[iPos][jPos];

        for (int i = 0; i < rows; i++) {
            if (i != iPos) {
                double searchValue = matrix[i][jPos];

                if (searchValue < refValue) {
                    flag = false;

                }
            }
        }

        return flag;
    }

    /**
     *
     * @param coloumn1
     * @param coloumn2
     */
    public void swapColumns(int coloumn1, int coloumn2) {

        double tempValue;
        //column exchange

        if (coloumn1 <= getColumnDimension() && coloumn2 <= getColumnDimension()) {
            for (int i = 0; i < rows; i++) {
                tempValue = matrix[i][coloumn1];
                matrix[i][coloumn1] = matrix[i][coloumn2];
                matrix[i][coloumn2] = tempValue;
                tempValue = 0.0;
            }
        } else {
            try {
                throw new CDKException("Index out of range" + coloumn1 + ", " + coloumn2);
            } catch (CDKException ex) {
                LOGGER.error(ex);
            }
        }
    }

    /**
     *
     * @param row1
     * @param row2
     * @throws org.openscience.cdk.exception.CDKException
     */
    public void swapRows(int row1, int row2) throws CDKException {

        double tempValue;
        //row exchange

        if (row1 <= getRowDimension() && row2 <= getRowDimension()) {
            for (int i = 0; i < columns; i++) {
                tempValue = matrix[row1][i];
                matrix[row1][i] = matrix[row2][i];
                matrix[row2][i] = tempValue;
            }
        } else {
            throw new CDKException("Index out of range" + row1 + ", " + row2);
        }

    }

    /**
     *
     * @param row chosen row
     * @param col chosen col
     */
    public void pivot(int row, int col) {

        //label pivot
        double tempValue;
        //column exchange
        for (int i = 0; i < rows; i++) {
            tempValue = matrix[i][row];
            matrix[i][row] = matrix[i][col];
            matrix[i][col] = tempValue;
        }
        //row exchange
        for (int i = 0; i < columns; i++) {
            tempValue = matrix[row][i];
            matrix[row][i] = matrix[col][i];
            matrix[col][i] = tempValue;
        }
    }

    /**
     * Normalizes the vectors of this matrix.
     *
     * @param S
     * @return
     */
    public EBIMatrix normalize(EBIMatrix S) {
        int p, q, i, j;
        double length;
        EBIMatrix result = duplicate();
        for (p = 0; p < columns; p++) // Loops over all vectors
        {
            // Calculates the normalization factor
            length = 0;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < rows; j++) {
                    length += result.matrix[i][p] * result.matrix[j][p] * S.matrix[i][j];
                }
            }

            length = sqrt(length);

            // Normalizes the vector
            if (length != 0d) {
                for (q = 0; q < rows; q++) {
                    result.matrix[q][p] /= length;
                }
            } else {
                LOGGER.debug("Warning(orthonormalize):" + (p + 1) + ". Vector has length null");
            }
        }
        return result;
    }

    /**
     * Multiplies a scalar with this EBIMatrix.
     *
     * @param a
     * @return
     */
    public EBIMatrix mul(double a) {
        EBIMatrix result = new EBIMatrix(rows, columns);
        int i, j;
        for (i = 0; i < rows; i++) {
            for (j = 0; j < columns; j++) {
                result.matrix[i][j] = matrix[i][j] * a;
            }
        }

        return result;
    }

    /**
     * Multiplies a Vector with this EBIMatrix.
     *
     * @param a
     * @return
     */
    public List<Double> mul(List<Double> a) {
        if ((a == null)
                || (columns != a.size())) {
            return null;
        }

        List<Double> result = new ArrayList<>(rows);
        int i, j;
        double sum;
        for (i = 0; i < rows; i++) {
            sum = 0;
            for (j = 0; j < columns; j++) {
                sum += (matrix[i][j]) * (a.get(j));
            }
            result.add(i, sum);
        }
        return result;
    }

    /**
     * Multiplies this EBIMatrix with another one.
     *
     * @param b
     * @return
     */
    public EBIMatrix mul(EBIMatrix b) {
        if ((b == null)
                || (columns != b.getRowDimension())) {
            return null;
        }

        EBIMatrix result = new EBIMatrix(rows, b.getColumnDimension());
        int i, j, k;
        double sum;
        for (i = 0; i < rows; i++) {
            for (k = 0; k < b.getColumnDimension(); k++) {
                sum = 0;
                for (j = 0; j < columns; j++) {
                    sum += matrix[i][j] * b.matrix[j][k];
                }
                result.matrix[i][k] = sum;
            }
        }

        return result;
    }

    /*
    * ------------------------ Private Methods ------------------------
     */
    /**
     * dataheck if size(matrix) == size(B) *
     */
    private void checkMatrixDimensions(EBIMatrix B) {
        if (B.getRowDimension() != rows || B.getColumnDimension() != columns) {
            throw new IllegalArgumentException("EBIMatrix dimensions must agree.");
        }
    }

    /**
     * Element-by-element multiplication, data = matrix.*B
     *
     * @param B another matrix
     * @return matrix.*B
     */
    public EBIMatrix arrayTimes(EBIMatrix B) {
        checkMatrixDimensions(B);
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = matrix[i][j] * B.matrix[i][j];
            }
        }
        return X;
    }

    /**
     * Element-by-element multiplication in place, matrix = matrix.*B
     *
     * @param B another matrix
     * @return matrix.*B
     */
    public EBIMatrix arrayTimesEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] *= B.matrix[i][j];
            }
        }
        return this;
    }

    /**
     * Element-by-element right division, data = matrix./B
     *
     * @param B another matrix
     * @return matrix./B
     */
    public EBIMatrix arrayRightDivide(EBIMatrix B) {
        checkMatrixDimensions(B);
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = matrix[i][j] / B.matrix[i][j];
            }
        }
        return X;
    }

    /**
     * Element-by-element right division in place, matrix = matrix./B
     *
     * @param B another matrix
     * @return matrix./B
     */
    public EBIMatrix arrayRightDivideEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] /= B.matrix[i][j];
            }
        }
        return this;
    }

    /**
     * Element-by-element left division, data = matrix.\B
     *
     * @param B another matrix
     * @return matrix.\B
     */
    public EBIMatrix arrayLeftDivide(EBIMatrix B) {
        checkMatrixDimensions(B);
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = B.matrix[i][j] / matrix[i][j];
            }
        }
        return X;
    }

    /**
     * Element-by-element left division in place, matrix = matrix.\B
     *
     * @param B another matrix
     * @return matrix.\B
     */
    public EBIMatrix arrayLeftDivideEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = B.matrix[i][j] / matrix[i][j];
            }
        }
        return this;
    }

    /**
     * Multiply a matrix by a scalar, data = s*matrix
     *
     * @param s scalar
     * @return s*matrix
     */
    public EBIMatrix times(double s) {
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = s * matrix[i][j];
            }
        }
        return X;
    }

    /**
     * Multiply a matrix by a scalar in place, matrix = s*matrix
     *
     * @param s scalar
     * @return replace matrix by s*matrix
     */
    public EBIMatrix timesEquals(double s) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = s * matrix[i][j];
            }
        }
        return this;
    }

    /**
     * Linear algebraic matrix multiplication, matrix * B
     *
     * @param B another matrix
     * @return EBIMatrix product, matrix * B
     * @exception IllegalArgumentException EBIMatrix inner dimensions must
     * agree.
     */
    public EBIMatrix times(EBIMatrix B) {
        if (B.getRowDimension() != columns) {
            throw new IllegalArgumentException("EBIMatrix inner dimensions must agree.");
        }
        int bCols = B.getColumnDimension();
        EBIMatrix X = new EBIMatrix(rows, bCols);
        double[][] C = X.getArray();
        double[][] Barray = B.getArray();
        // Cache-friendly i,k,j loop order: iterate over rows of A, then
        // for each element A[i][k], scatter-multiply across B's row k.
        // This accesses both A and B in row-major order, maximizing L1/L2
        // cache utilization and avoiding the column-copy temporary array.
        for (int i = 0; i < rows; i++) {
            double[] Arowi = matrix[i];
            double[] Crowi = C[i];
            for (int k = 0; k < columns; k++) {
                double aik = Arowi[k];
                if (aik == 0.0) continue;
                double[] Browk = Barray[k];
                for (int j = 0; j < bCols; j++) {
                    Crowi[j] += aik * Browk[j];
                }
            }
        }
        return X;
    }

    /**
     * LU Decomposition
     *
     * @return LUDecomposition
     * @see LUDecomposition
     */
    public LUDecomposition lu() {
        return new LUDecomposition(this);
    }

    /**
     * QR Decomposition
     *
     * @return QRDecomposition
     * @see QRDecomposition
     */
    public QRDecomposition qr() {
        return new QRDecomposition(this);
    }

    /**
     * dataholesky Decomposition
     *
     * @return dataholeskyDecomposition
     * @see CholeskyDecomposition
     */
    public CholeskyDecomposition chol() {
        return new CholeskyDecomposition(this);
    }

    /**
     * Singular Value Decomposition
     *
     * @return SingularValueDecomposition
     * @see SingularValueDecomposition
     */
    public SingularValueDecomposition svd() {
        return new SingularValueDecomposition(this);
    }

    /**
     * Eigenvalue Decomposition
     *
     * @return EigenvalueDecomposition
     * @see EigenvalueDecomposition
     */
    public EigenvalueDecomposition eig() {
        return new EigenvalueDecomposition(this);
    }

    /**
     * Solve matrix*result = B
     *
     * @param B right hand side
     * @return solution if matrix is square, least squares solution otherwise
     */
    public EBIMatrix solve(EBIMatrix B) {
        return (rows == columns ? (new LUDecomposition(this)).solve(B) : (new QRDecomposition(this)).solve(B));
    }

    /**
     * Solve result*matrix = B, which is also matrix'*result' = B'
     *
     * @param B right hand side
     * @return solution if matrix is square, least squares solution otherwise.
     */
    public EBIMatrix solveTranspose(EBIMatrix B) {
        return transpose().solve(B.transpose());
    }

    /**
     * EBIMatrix inverse or pseudoinverse
     *
     * @return inverse(matrix) if matrix is square, pseudoinverse otherwise.
     */
    public EBIMatrix inverse() {
        return solve(identity(rows, rows));
    }

    /**
     * EBIMatrix determinant
     *
     * @return determinant
     */
    public double det() {
        return new LUDecomposition(this).det();
    }

    /**
     * EBIMatrix rank
     *
     * @return effective numerical rank, obtained from SVD.
     */
    public int rank() {
        return new SingularValueDecomposition(this).rank();
    }

    /**
     * EBIMatrix condition (2 norm)
     *
     * @return ratio of largest to smallest singular value.
     */
    public double cond() {
        return new SingularValueDecomposition(this).cond();
    }

    /**
     * EBIMatrix trace.
     *
     * @return sum of the diagonal elements.
     */
    public double trace() {
        double t = 0;
        for (int i = 0; i < min(rows, columns); i++) {
            t += matrix[i][i];
        }
        return t;
    }

    /**
     * Diagonalize this matrix with the Jacobi algorithm.
     *
     * @param nrot dataount of max. rotations
     * @return EBIMatrix m, with m^t * this * m = diagonal
     *
     */
    public EBIMatrix diagonalize(int nrot) {
        EBIMatrix m = duplicate();
        if (m.getRowDimension() != m.getColumnDimension()) {
            LOGGER.debug("EBIMatrix.diagonal: Sizes mismatched");
            return null;
        }
        int n = m.getRowDimension();

        int j, iq, ip, i;

        double tresh, theta, tau, t, sm, s, h, g, c;
        double[] b, z;

        EBIMatrix v = new EBIMatrix(columns, columns);
        List<Double> d = new ArrayList<>(columns);

        b = new double[n + 1];
        z = new double[n + 1];

        for (ip = 0; ip < n; ip++) {
            for (iq = 0; iq < n; iq++) {
                v.matrix[ip][iq] = 0.0;
            }
            v.matrix[ip][ip] = 1.0;
        }

        for (ip = 0; ip < n; ip++) {
            d.add(ip, m.matrix[ip][ip]);
            b[ip] = m.matrix[ip][ip];
            z[ip] = 0.0;
        }

        nrot = 0;
        for (i = 1; i <= 50; i++) {
            sm = 0.0;
            for (ip = 0; ip < n - 1; ip++) {
                for (iq = ip + 1; iq < n; iq++) {
                    sm += abs(m.matrix[ip][iq]);
                }
            }

            // Ready ??
            if (sm == 0.0) {
                return v;
            }

            if (i < 4) {
                tresh = 0.2 * sm / (n * n);
            } else {
                tresh = 0.0;
            }

            for (ip = 0; ip < n - 1; ip++) {
                for (iq = ip + 1; iq < n; iq++) {
                    g = 100.0 * abs(m.matrix[ip][iq]);
                    if ((i > 4) && (abs(d.get(ip) + g) == abs(d.get(ip)))
                            && (abs(d.get(ip) + g) == abs(d.get(ip)))) {
                        m.matrix[ip][iq] = 0.0;
                    } else if (abs(m.matrix[ip][iq]) > tresh) {
                        h = d.get(iq) - d.get(ip);
                        if (abs(h) + g == abs(h)) {
                            t = (m.matrix[ip][iq]) / h;
                        } else {
                            theta = 0.5 * h / (m.matrix[ip][iq]);
                            t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                            if (theta < 0.0) {
                                t = -t;
                            }
                        }
                        c = 1.0 / sqrt(1 + t * t);
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * m.matrix[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d.set(ip, d.get(ip) - h);
                        d.set(iq, d.get(iq) + h);
                        m.matrix[ip][iq] = 0.0;

                        // Case of rotaions 1<=j<p
                        for (j = 0; j < ip; j++) {
                            g = m.matrix[j][ip];
                            h = m.matrix[j][iq];
                            m.matrix[j][ip] = g - s * (h + g * tau);
                            m.matrix[j][iq] = h + s * (g - h * tau);
                        }
                        // Case of rotaions p<j<q
                        for (j = ip + 1; j < iq; j++) {
                            g = m.matrix[ip][j];
                            h = m.matrix[j][iq];
                            m.matrix[ip][j] = g - s * (h + g * tau);
                            m.matrix[j][iq] = h + s * (g - h * tau);
                        }
                        // Case of rotaions q<j<=n
                        for (j = iq + 1; j < n; j++) {
                            g = m.matrix[ip][j];
                            h = m.matrix[iq][j];
                            m.matrix[ip][j] = g - s * (h + g * tau);
                            m.matrix[iq][j] = h + s * (g - h * tau);
                        }
                        for (j = 0; j < n; j++) {
                            g = v.matrix[j][ip];
                            h = v.matrix[j][iq];
                            v.matrix[j][ip] = g - s * (h + g * tau);
                            v.matrix[j][iq] = h + s * (g - h * tau);
                        }
                        ++nrot;
                    }
                }
            }

            for (ip = 0; ip < n; ip++) {
                b[ip] += z[ip];
                d.set(ip, b[ip]);
                z[ip] = 0.0;
            }
        }
        LOGGER.debug("Too many iterations in routine JACOBI");
        return v;
    }

    /**
     * Orthonormalize the vectors of this matrix by Gram-Schmidt.
     *
     * @param S
     * @return
     * @keyword orthonormalization
     * @keyword Gram-Schmidt algorithm
     */
    public EBIMatrix orthonormalize(EBIMatrix S) {
        int p, q, k, i, j;
        double innersum;
        double length;
        //EBIMatrix scr = S.mul(this);
        EBIMatrix result = duplicate();
        for (p = 0; p < columns; p++) // Loops over all vectors
        {
            for (i = 0; i < rows; i++) {
                result.matrix[i][p] = matrix[i][p];
            }

            for (k = 0; k < p; k++) // Substracts the previous vector 
            {
                // First the calculation of the product <phi_p|phi_k>=length
                length = 0;
                for (i = 0; i < rows; i++) // Loops over all vectors
                {
                    innersum = 0;
                    for (j = 0; j < rows; j++) {
                        innersum += result.matrix[j][p] * S.matrix[i][j];
                    }
                    length += result.matrix[i][k] * innersum;
                }

                // Then the substraction of  phi_k*length
                for (q = 0; q < rows; q++) {
                    result.matrix[q][p] -= result.matrix[q][k] * length;
                }
            }

            // Calculates the integral for normalization
            length = 0;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < rows; j++) {
                    length += result.matrix[i][p] * result.matrix[j][p] * S.matrix[i][j];
                }
            }

            length = sqrt(length);

            // Normalizes the vector
            if (length != 0d) {
                for (q = 0; q < rows; q++) {
                    result.matrix[q][p] /= length;
                }
            } else {
                LOGGER.debug("Warning(orthonormalize):" + (p + 1) + ". Vector has length null");
            }
        }
        return result;
    }

    /**
     * Print the matrix to stdout. Line the elements up in columns. Use the
     * format object, and right justify within columns of width characters. Note
     * that is the matrix is to be read back in, you probably will want to use a
     * NumberFormat that is set to US Locale.
     *
     * @param format matrix Formatting object for individual elements.
     * @param width Field width for each column.
     * @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    public void print(NumberFormat format, int width) {
        print(new PrintWriter(System.out, true), format, width);
    }

    /**
     * Print the matrix to stdout. Line the elements up in columns with a
     * Fortran-like 'Fw.d' style format.
     *
     * @param w dataolumn width.
     * @param d Number of digits after the decimal.
     */
    public void print(int w, int d) {
        print(new PrintWriter(System.out, true), w, d);
    }

    /**
     * Print the matrix to the output stream. Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
     *
     * @param output Output stream.
     * @param w dataolumn width.
     * @param d Number of digits after the decimal.
     */
    public void print(PrintWriter output, int w, int d) {
        DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        print(output, format, w + 2);
    }

//    DecimalFormat is a little disappointing coming from Fortran or data's printf.
//    Since it doesn't pad on the left, the elements will come out different
//    widths.  Consequently, we'll pass the desired column width in as an
//    argument and do the extra padding ourselves.
    /**
     * Print the matrix to the output stream. Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters. Note that is the matrix is to be read back in, you probably
     * will want to use a NumberFormat that is set to US Locale.
     *
     * @param output the output stream.
     * @param format matrix formatting object to format the matrix elements
     * @param width dataolumn width.
     * @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    public void print(PrintWriter output, NumberFormat format, int width) {
        output.println();  // start on new line.
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                String s = format.format(matrix[i][j]); // format the number
                int padding = max(1, width - s.length()); // At _least_ 1 space
                for (int k = 0; k < padding; k++) {
                    output.print(' ');
                }
                output.print(s);
            }
            output.println();
        }
        output.println();   // end with blank line.
    }

    /**
     * Return a matrix as a String.
     *
     * @return
     */
    @Override
    public String toString() {
        if ((rows <= 0) || (columns <= 0)) {
            return "[]";
        }
        int i, j;
        DecimalFormat format = new DecimalFormat("00.0000");
        format.setPositivePrefix("+");

        StringBuilder str = new StringBuilder();
        for (i = 0; i < (rows - 1); i++) {
            for (j = 0; j < (columns - 1); j++) {
                if (round(matrix[i][j] * 10000) != 0) {
                    str.append(format.format(matrix[i][j])).append(" ");
                } else {
                    str.append("-------- ");
                }
            }
            if (round(matrix[i][columns - 1] * 10000) != 0) {
                str.append(format.format(matrix[i][columns - 1])).append(NEW_LINE);
            } else {
                str.append("--------").append(NEW_LINE);
            }
        }
        for (j = 0; j < (columns - 1); j++) {
            if (round(matrix[rows - 1][j] * 10000) != 0) {
                str.append(format.format(matrix[rows - 1][j])).append(" ");
            } else {
                str.append("-------- ");
            }
        }
        if (round(matrix[rows - 1][columns - 1] * 10000) != 0) {
            str.append(format.format(matrix[rows - 1][columns - 1]));
        } else {
            str.append("-------- ");
        }
        return str.toString();
    }

    /**
     *
     * @return
     */
    public double contraction() {
        int i, j;
        double result = 0d;
        for (i = 0; i < getRowDimension(); i++) {
            for (j = 0; j < getColumnDimension(); j++) {
                result += matrix[i][j];
            }
        }
        return result;
    }

    /**
     * Similar transformation Ut * M * U
     *
     * @param U
     * @return
     */
    public EBIMatrix similar(EBIMatrix U) {
        EBIMatrix result = new EBIMatrix(U.getColumnDimension(), U.getColumnDimension());
        double sum, innersum;
        for (int i = 0; i < U.getColumnDimension(); i++) {
            for (int j = 0; j < U.getColumnDimension(); j++) {
                sum = 0d;
                for (int k = 0; k < U.getColumnDimension(); k++) {
                    innersum = 0d;
                    for (int l = 0; l < U.getColumnDimension(); l++) {
                        innersum += matrix[k][l] * U.matrix[l][j];
                    }
                    sum += U.matrix[k][i] * innersum;
                }
                result.matrix[i][j] = sum;
            }
        }
        return result;
    }

    /**
     * One norm
     *
     * @return maximum column sum.
     */
    public double norm1() {
        double f = 0;
        for (int j = 0; j < columns; j++) {
            double s = 0;
            for (int i = 0; i < rows; i++) {
                s += abs(matrix[i][j]);
            }
            f = max(f, s);
        }
        return f;
    }

    /**
     * Two norm
     *
     * @return maximum singular value.
     */
    public double norm2() {
        return (new SingularValueDecomposition(this).norm2());
    }

    /**
     * Infinity norm
     *
     * @return maximum row sum.
     */
    public double normInf() {
        double f = 0;
        for (int i = 0; i < rows; i++) {
            double s = 0;
            for (int j = 0; j < columns; j++) {
                s += abs(matrix[i][j]);
            }
            f = max(f, s);
        }
        return f;
    }

    /**
     * Frobenius norm
     *
     * @return sqrt of sum of squares of all elements.
     */
    public double normF() {
        double f = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                f = hypot(f, matrix[i][j]);
            }
        }
        return f;
    }

    /**
     * Unary minus
     *
     * @return -matrix
     */
    public EBIMatrix uminus() {
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = -matrix[i][j];
            }
        }
        return X;
    }

    /**
     * data = matrix + B
     *
     * @param B another matrix
     * @return matrix + B
     */
    public EBIMatrix plus(EBIMatrix B) {
        checkMatrixDimensions(B);
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = matrix[i][j] + B.matrix[i][j];
            }
        }
        return X;
    }

    /**
     * matrix = matrix + B
     *
     * @param B another matrix
     * @return matrix + B
     */
    public EBIMatrix plusEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] += B.matrix[i][j];
            }
        }
        return this;
    }

    /**
     * data = matrix - B
     *
     * @param B another matrix
     * @return matrix - B
     */
    public EBIMatrix minus(EBIMatrix B) {
        checkMatrixDimensions(B);
        EBIMatrix X = new EBIMatrix(rows, columns);
        double[][] C = X.getArray();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                C[i][j] = matrix[i][j] - B.matrix[i][j];
            }
        }
        return X;
    }

    /**
     * matrix = matrix - B
     *
     * @param B another matrix
     * @return matrix - B
     */
    public EBIMatrix minusEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] -= B.matrix[i][j];
            }
        }
        return this;
    }


        // ==================== Inner class CholeskyDecomposition ====================

        /**
     * Cholesky Decomposition.
     * <P>
     * For a symmetric, positive definite matrix A, the Cholesky decomposition is an
     * lower triangular matrix L so that A = L*L'.
     * <P>
     * If the matrix is not symmetric or positive definite, the constructor returns
     * a partial decomposition and sets an internal flag that may be queried by the
     * isSPD() method.
     */
    public static class CholeskyDecomposition implements java.io.Serializable {

        private static final long serialVersionUID = 78619981017L;

        /* ------------------------
         Class variables
         * ------------------------ */
        /**
         * Array for internal storage of decomposition.
         *
         * @serial internal array storage.
         */
        private final double[][] L;
        /**
         * Row and column dimension (square matrix).
         *
         * @serial matrix dimension.
         */
        private final int n;
        /**
         * Symmetric and positive definite flag.
         *
         * @serial is symmetric and positive definite flag.
         */
        private boolean isspd;

        /* ------------------------
         Constructor
         * ------------------------ */
        /**
         * Cholesky algorithm for symmetric and positive definite matrix.
         *
         * @param Arg
         */
        public CholeskyDecomposition(EBIMatrix Arg) {

            // Initialize.
            double[][] A = Arg.getArray();
            n = Arg.getRowDimension();
            L = new double[n][n];
            isspd = (Arg.getColumnDimension() == n);
            // Main loop.
            for (int j = 0; j < n; j++) {
                double[] Lrowj = L[j];
                double d = 0.0;
                for (int k = 0; k < j; k++) {
                    double[] Lrowk = L[k];
                    double s = 0.0;
                    for (int i = 0; i < k; i++) {
                        s += Lrowk[i] * Lrowj[i];
                    }
                    Lrowj[k] = s = (A[j][k] - s) / L[k][k];
                    d += s * s;
                    isspd &= (A[k][j] == A[j][k]);
                }
                d = A[j][j] - d;
                isspd &= (d > 0.0);
                L[j][j] = sqrt(max(d, 0.0));
                for (int k = j + 1; k < n; k++) {
                    L[j][k] = 0.0;
                }
            }
        }

        /* ------------------------
         Temporary, experimental code.
         * ------------------------ *\

         \** Right Triangular Cholesky Decomposition.
         <P>
         For a symmetric, positive definite matrix A, the Right Cholesky
         decomposition is an upper triangular matrix R so that A = R'*R.
         This constructor computes R with the Fortran inspired column oriented
         algorithm used in LINPACK and MATLAB.  In Java, we suspect a row oriented,
         lower triangular decomposition is faster.  We have temporarily included
         this constructor here until timing experiments confirm this suspicion.
         *\

         \** Array for internal storage of right triangular decomposition. **\
         private transient double[][] R;

         \** Cholesky algorithm for symmetric and positive definite matrix.
         @param  A           Square, symmetric matrix.
         @param  rightflag   Actual value ignored.
         @return             Structure to access R and isspd flag.
         *\

         public CholeskyDecomposition (Matrix Arg, int rightflag) {
         // Initialize.
         double[][] A = Arg.getArray();
         n = Arg.getColumnDimension();
         R = new double[n][n];
         isspd = (Arg.getColumnDimension() == n);
         // Main loop.
         for (int j = 0; j < n; j++) {
         double d = 0.0;
         for (int k = 0; k < j; k++) {
         double s = A[k][j];
         for (int i = 0; i < k; i++) {
         s = s - R[i][k]*R[i][j];
         }
         R[k][j] = s = s/R[k][k];
         d = d + s*s;
         isspd = isspd && (A[k][j] == A[j][k]);
         }
         d = A[j][j] - d;
         isspd = isspd && (d > 0.0);
         R[j][j] = Math.sqrt(Math.max(d,0.0));
         for (int k = j+1; k < n; k++) {
         R[k][j] = 0.0;
         }
         }
         }

         \** Return upper triangular factor.
         @return     R
         *\

         public Matrix getR () {
         return new Matrix(R,n,n);
         }

         \* ------------------------
         End of temporary code.
         * ------------------------ */

     /* ------------------------
         Public Methods
         * ------------------------ */
        /**
         * Is the matrix symmetric and positive definite?
         *
         * @return true if A is symmetric and positive definite.
         */
        public boolean isSPD() {
            return isspd;
        }

        /**
         * Return triangular factor.
         *
         * @return L
         */
        public EBIMatrix getL() {
            return new EBIMatrix(L, n, n);
        }

        /**
         * Solve A*X = B
         *
         * @param B A Matrix with as many rows as A and any number of columns.
         * @return X so that L*L'*X = B
         * @exception IllegalArgumentException Matrix row dimensions must agree.
         * @exception RuntimeException Matrix is not symmetric positive definite.
         */
        public EBIMatrix solve(EBIMatrix B) {
            if (B.getRowDimension() != n) {
                throw new IllegalArgumentException("Matrix row dimensions must agree.");
            }
            if (!isspd) {
                throw new RuntimeException("Matrix is not symmetric positive definite.");
            }

            // Copy right hand side.
            double[][] X = B.getArrayCopy();
            int nx = B.getColumnDimension();

            // Solve L*Y = B;
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < nx; j++) {
                    for (int i = 0; i < k; i++) {
                        X[k][j] -= X[i][j] * L[k][i];
                    }
                    X[k][j] /= L[k][k];
                }
            }

            // Solve L'*X = Y;
            for (int k = n - 1; k >= 0; k--) {
                for (int j = 0; j < nx; j++) {
                    for (int i = k + 1; i < n; i++) {
                        X[k][j] -= X[i][j] * L[i][k];
                    }
                    X[k][j] /= L[k][k];
                }
            }

            return new EBIMatrix(X, n, nx);
        }
    }

        // ==================== Inner class EigenvalueDecomposition ====================

        /**
     * Eigenvalues and eigenvectors of a real matrix.
     * <P>
     * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
     * and the eigenvector matrix V is orthogonal. I.e. A =
     * V.times(D.times(V.transpose())) and V.times(V.transpose()) equals the
     * identity matrix.
     * <P>
     * If A is not symmetric, then the eigenvalue matrix D is block diagonal with
     * the real eigenvalues in 1-by-1 blocks and any complex eigenvalues, lambda +
     * i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns of V represent
     * the eigenvectors in the sense that A*V = V*D, i.e. A.times(V) equals
     * V.times(D). The matrix V may be badly conditioned, or even singular, so the
     * validity of the equation A = V*D*inverse(V) depends upon V.cond().
     *
     */
    public static class EigenvalueDecomposition implements java.io.Serializable {

        private static final long serialVersionUID = 17869981017L;
        /* ------------------------
         Class variables
         * ------------------------ */
        /**
         * Row and column dimension (square matrix).
         *
         * @serial matrix dimension.
         */
        private final int n;
        /**
         * Symmetry flag.
         *
         * @serial internal symmetry flag.
         */
        private boolean issymmetric;
        /**
         * Arrays for internal storage of eigenvalues.
         *
         * @serial internal storage of eigenvalues.
         */
        private final double[] d, e;
        /**
         * Array for internal storage of eigenvectors.
         *
         * @serial internal storage of eigenvectors.
         */
        private final double[][] V;
        /**
         * Array for internal storage of nonsymmetric Hessenberg form.
         *
         * @serial internal storage of nonsymmetric Hessenberg form.
         */
        private double[][] H;
        /**
         * Working storage for nonsymmetric algorithm.
         *
         * @serial working storage for nonsymmetric algorithm.
         */
        private double[] ort;

        /* ------------------------
         Private Methods
         * ------------------------ */    private transient double cdivr, cdivi;

        /* ------------------------
        Constructor
         * ------------------------ */
        /**
         * Check for symmetry, then construct the eigenvalue decomposition
         *
         * @param Arg
         */
        public EigenvalueDecomposition(EBIMatrix Arg) {
            double[][] A = Arg.getArray();
            n = Arg.getColumnDimension();
            V = new double[n][n];
            d = new double[n];
            e = new double[n];

            issymmetric = true;
            for (int j = 0; (j < n) && issymmetric; j++) {
                for (int i = 0; (i < n) && issymmetric; i++) {
                    issymmetric = (A[i][j] == A[j][i]);
                }
            }

            if (issymmetric) {
                for (int i = 0; i < n; i++) {
                    arraycopy(A[i], 0, V[i], 0, n);
                }

                // Tridiagonalize.
                tred2();

                // Diagonalize.
                tql2();

            } else {
                H = new double[n][n];
                ort = new double[n];

                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        H[i][j] = A[i][j];
                    }
                }

                // Reduce to Hessenberg form.
                orthes();

                // Reduce Hessenberg to real Schur form.
                hqr2();
            }
        }

        /* ------------------------
        Private Methods
        * ------------------------ */
        // Symmetric Householder reduction to tridiagonal form.
        private void tred2() {
            arraycopy(V[n - 1], 0, d, 0, n);

            // Householder reduction to tridiagonal form.
            for (int i = n - 1; i > 0; i--) {

                // Scale to avoid under/overflow.
                double scale = 0.0;
                double h = 0.0;
                for (int k = 0; k < i; k++) {
                    scale += abs(d[k]);
                }
                if (scale == 0.0) {
                    e[i] = d[i - 1];
                    for (int j = 0; j < i; j++) {
                        d[j] = V[i - 1][j];
                        V[i][j] = 0.0;
                        V[j][i] = 0.0;
                    }
                } else {

                    // Generate Householder vector.
                    for (int k = 0; k < i; k++) {
                        d[k] /= scale;
                        h += d[k] * d[k];
                    }
                    double f = d[i - 1];
                    double g = sqrt(h);
                    if (f > 0) {
                        g = -g;
                    }
                    e[i] = scale * g;
                    h -= f * g;
                    d[i - 1] = f - g;
                    for (int j = 0; j < i; j++) {
                        e[j] = 0.0;
                    }

                    // Apply similarity transformation to remaining columns.
                    for (int j = 0; j < i; j++) {
                        f = d[j];
                        V[j][i] = f;
                        g = e[j] + V[j][j] * f;
                        for (int k = j + 1; k <= i - 1; k++) {
                            g += V[k][j] * d[k];
                            e[k] += V[k][j] * f;
                        }
                        e[j] = g;
                    }
                    f = 0.0;
                    for (int j = 0; j < i; j++) {
                        e[j] /= h;
                        f += e[j] * d[j];
                    }
                    double hh = f / (h + h);
                    for (int j = 0; j < i; j++) {
                        e[j] -= hh * d[j];
                    }
                    for (int j = 0; j < i; j++) {
                        f = d[j];
                        g = e[j];
                        for (int k = j; k <= i - 1; k++) {
                            V[k][j] -= (f * e[k] + g * d[k]);
                        }
                        d[j] = V[i - 1][j];
                        V[i][j] = 0.0;
                    }
                }
                d[i] = h;
            }

            // Accumulate transformations.
            for (int i = 0; i < n - 1; i++) {
                V[n - 1][i] = V[i][i];
                V[i][i] = 1.0;
                double h = d[i + 1];
                if (h != 0.0) {
                    for (int k = 0; k <= i; k++) {
                        d[k] = V[k][i + 1] / h;
                    }
                    for (int j = 0; j <= i; j++) {
                        double g = 0.0;
                        for (int k = 0; k <= i; k++) {
                            g += V[k][i + 1] * V[k][j];
                        }
                        for (int k = 0; k <= i; k++) {
                            V[k][j] -= g * d[k];
                        }
                    }
                }
                for (int k = 0; k <= i; k++) {
                    V[k][i + 1] = 0.0;
                }
            }
            for (int j = 0; j < n; j++) {
                d[j] = V[n - 1][j];
                V[n - 1][j] = 0.0;
            }
            V[n - 1][n - 1] = 1.0;
            e[0] = 0.0;
        }

        // Symmetric tridiagonal QL algorithm.
        private void tql2() {

            //  This is derived from the Algol procedures tql2, by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.
            for (int i = 1; i < n; i++) {
                e[i - 1] = e[i];
            }
            e[n - 1] = 0.0;

            double f = 0.0;
            double tst1 = 0.0;
            double eps = pow(2.0, -52.0);
            for (int l = 0; l < n; l++) {

                // Find small subdiagonal element
                tst1 = max(tst1, abs(d[l]) + abs(e[l]));
                int m = l;
                while (m < n) {
                    if (abs(e[m]) <= eps * tst1) {
                        break;
                    }
                    m++;
                }

                // If m == l, d[l] is an eigenvalue,
                // otherwise, iterate.
                if (m > l) {
                    int iter = 0;
                    do {
                        iter += 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        double g = d[l];
                        double p = (d[l + 1] - g) / (2.0 * e[l]);
                        double r = hypot(p, 1.0);
                        if (p < 0) {
                            r = -r;
                        }
                        d[l] = e[l] / (p + r);
                        d[l + 1] = e[l] * (p + r);
                        double dl1 = d[l + 1];
                        double h = g - d[l];
                        for (int i = l + 2; i < n; i++) {
                            d[i] -= h;
                        }
                        f += h;

                        // Implicit QL transformation.
                        p = d[m];
                        double c = 1.0;
                        double c2 = c;
                        double c3 = c;
                        double el1 = e[l + 1];
                        double s = 0.0;
                        double s2 = 0.0;
                        for (int i = m - 1; i >= l; i--) {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e[i];
                            h = c * p;
                            r = hypot(p, e[i]);
                            e[i + 1] = s * r;
                            s = e[i] / r;
                            c = p / r;
                            p = c * d[i] - s * g;
                            d[i + 1] = h + s * (c * g + s * d[i]);

                            // Accumulate transformation.
                            for (int k = 0; k < n; k++) {
                                h = V[k][i + 1];
                                V[k][i + 1] = s * V[k][i] + c * h;
                                V[k][i] = c * V[k][i] - s * h;
                            }
                        }
                        p = -s * s2 * c3 * el1 * e[l] / dl1;
                        e[l] = s * p;
                        d[l] = c * p;

                        // Check for convergence.
                    } while (abs(e[l]) > eps * tst1);
                }
                d[l] += f;
                e[l] = 0.0;
            }

            // Sort eigenvalues and corresponding vectors.
            for (int i = 0; i < n - 1; i++) {
                int k = i;
                double p = d[i];
                for (int j = i + 1; j < n; j++) {
                    if (d[j] < p) {
                        k = j;
                        p = d[j];
                    }
                }
                if (k != i) {
                    d[k] = d[i];
                    d[i] = p;
                    for (int j = 0; j < n; j++) {
                        p = V[j][i];
                        V[j][i] = V[j][k];
                        V[j][k] = p;
                    }
                }
            }
        }

        // Nonsymmetric reduction to Hessenberg form.
        private void orthes() {

            //  This is derived from the Algol procedures orthes and ortran,
            //  by Martin and Wilkinson, Handbook for Auto. Comp.,
            //  Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutines in EISPACK.
            int low = 0;
            int high = n - 1;

            for (int m = low + 1; m <= high - 1; m++) {

                // Scale column.
                double scale = 0.0;
                for (int i = m; i <= high; i++) {
                    scale += abs(H[i][m - 1]);
                }
                if (scale != 0.0) {

                    // Compute Householder transformation.
                    double h = 0.0;
                    for (int i = high; i >= m; i--) {
                        ort[i] = H[i][m - 1] / scale;
                        h += ort[i] * ort[i];
                    }
                    double g = sqrt(h);
                    if (ort[m] > 0) {
                        g = -g;
                    }
                    h -= ort[m] * g;
                    ort[m] -= g;

                    // Apply Householder similarity transformation
                    // H = (I-u*u'/h)*H*(I-u*u')/h)
                    for (int j = m; j < n; j++) {
                        double f = 0.0;
                        for (int i = high; i >= m; i--) {
                            f += ort[i] * H[i][j];
                        }
                        f /= h;
                        for (int i = m; i <= high; i++) {
                            H[i][j] -= f * ort[i];
                        }
                    }

                    for (int i = 0; i <= high; i++) {
                        double f = 0.0;
                        for (int j = high; j >= m; j--) {
                            f += ort[j] * H[i][j];
                        }
                        f /= h;
                        for (int j = m; j <= high; j++) {
                            H[i][j] -= f * ort[j];
                        }
                    }
                    ort[m] = scale * ort[m];
                    H[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    V[i][j] = (i == j ? 1.0 : 0.0);
                }
            }

            for (int m = high - 1; m >= low + 1; m--) {
                if (H[m][m - 1] != 0.0) {
                    for (int i = m + 1; i <= high; i++) {
                        ort[i] = H[i][m - 1];
                    }
                    for (int j = m; j <= high; j++) {
                        double g = 0.0;
                        for (int i = m; i <= high; i++) {
                            g += ort[i] * V[i][j];
                        }
                        // double division avoids possible underflow
                        g = (g / ort[m]) / H[m][m - 1];
                        for (int i = m; i <= high; i++) {
                            V[i][j] += g * ort[i];
                        }
                    }
                }
            }
        }    // Complex scalar division.

        private void cdiv(double xr, double xi, double yr, double yi) {
            double r, d1;
            if (abs(yr) > abs(yi)) {
                r = yi / yr;
                d1 = yr + r * yi;
                cdivr = (xr + r * xi) / d1;
                cdivi = (xi - r * xr) / d1;
            } else {
                r = yr / yi;
                d1 = yi + r * yr;
                cdivr = (r * xr + xi) / d1;
                cdivi = (r * xi - xr) / d1;
            }
        }

        // Nonsymmetric reduction from Hessenberg to real Schur form.
        private void hqr2() {

            //  This is derived from the Algol procedure hqr2,
            //  by Martin and Wilkinson, Handbook for Auto. Comp.,
            //  Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.
            // Initialize
            int nn = this.n;
            int n1 = nn - 1;
            int low = 0;
            int high = nn - 1;
            double eps = pow(2.0, -52.0);
            double exshift = 0.0;
            double p = 0.0, q = 0.0, r = 0.0, s = 0.0, z = 0.0, t, w, x, y;

            // Store roots isolated by balanc and compute matrix norm
            double norm = 0.0;
            for (int i = 0; i < nn; i++) {
                if (i < low | i > high) {
                    d[i] = H[i][i];
                    e[i] = 0.0;
                }
                for (int j = max(i - 1, 0); j < nn; j++) {
                    norm += abs(H[i][j]);
                }
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n1 >= low) {

                // Look for single small sub-diagonal element
                int l = n1;
                while (l > low) {
                    s = abs(H[l - 1][l - 1]) + abs(H[l][l]);
                    if (s == 0.0) {
                        s = norm;
                    }
                    if (abs(H[l][l - 1]) < eps * s) {
                        break;
                    }
                    l--;
                }

                // Check for convergence
                // One root found
                if (l == n1) {
                    H[n1][n1] += exshift;
                    d[n1] = H[n1][n1];
                    e[n1] = 0.0;
                    n1--;
                    iter = 0;

                    // Two roots found
                } else if (l == n1 - 1) {
                    w = H[n1][n1 - 1] * H[n1 - 1][n1];
                    p = (H[n1 - 1][n1 - 1] - H[n1][n1]) / 2.0;
                    q = p * p + w;
                    z = sqrt(abs(q));
                    H[n1][n1] += exshift;
                    H[n1 - 1][n1 - 1] += exshift;
                    x = H[n1][n1];

                    // Real pair
                    if (q >= 0) {
                        if (p >= 0) {
                            z = p + z;
                        } else {
                            z = p - z;
                        }
                        d[n1 - 1] = x + z;
                        d[n1] = d[n1 - 1];
                        if (z != 0.0) {
                            d[n1] = x - w / z;
                        }
                        e[n1 - 1] = 0.0;
                        e[n1] = 0.0;
                        x = H[n1][n1 - 1];
                        s = abs(x) + abs(z);
                        p = x / s;
                        q = z / s;
                        r = sqrt(p * p + q * q);
                        p /= r;
                        q /= r;

                        // Row modification
                        for (int j = n1 - 1; j < nn; j++) {
                            z = H[n1 - 1][j];
                            H[n1 - 1][j] = q * z + p * H[n1][j];
                            H[n1][j] = q * H[n1][j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= n1; i++) {
                            z = H[i][n1 - 1];
                            H[i][n1 - 1] = q * z + p * H[i][n1];
                            H[i][n1] = q * H[i][n1] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = low; i <= high; i++) {
                            z = V[i][n1 - 1];
                            V[i][n1 - 1] = q * z + p * V[i][n1];
                            V[i][n1] = q * V[i][n1] - p * z;
                        }

                        // Complex pair
                    } else {
                        d[n1 - 1] = x + p;
                        d[n1] = x + p;
                        e[n1 - 1] = z;
                        e[n1] = -z;
                    }
                    n1 -= 2;
                    iter = 0;

                    // No convergence yet
                } else {

                    // Form shift
                    x = H[n1][n1];
                    y = 0.0;
                    w = 0.0;
                    if (l < n1) {
                        y = H[n1 - 1][n1 - 1];
                        w = H[n1][n1 - 1] * H[n1 - 1][n1];
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10) {
                        exshift += x;
                        for (int i = low; i <= n1; i++) {
                            H[i][i] -= x;
                        }
                        s = abs(H[n1][n1 - 1]) + abs(H[n1 - 1][n1 - 2]);
                        x = y = 0.75 * s;
                        w = -0.4375 * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30) {
                        s = (y - x) / 2.0;
                        s = s * s + w;
                        if (s > 0) {
                            s = sqrt(s);
                            if (y < x) {
                                s = -s;
                            }
                            s = x - w / ((y - x) / 2.0 + s);
                            for (int i = low; i <= n1; i++) {
                                H[i][i] -= s;
                            }
                            exshift += s;
                            x = y = w = 0.964;
                        }
                    }

                    iter += 1;   // (Could check iteration count here.)

                    // Look for two consecutive small sub-diagonal elements
                    int m = n1 - 2;
                    while (m >= l) {
                        z = H[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
                        q = H[m + 1][m + 1] - z - r - s;
                        r = H[m + 2][m + 1];
                        s = abs(p) + abs(q) + abs(r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if (m == l) {
                            break;
                        }
                        if (abs(H[m][m - 1]) * (abs(q) + abs(r))
                                < eps * (abs(p) * (abs(H[m - 1][m - 1]) + abs(z)
                                + abs(H[m + 1][m + 1])))) {
                            break;
                        }
                        m--;
                    }

                    for (int i = m + 2; i <= n1; i++) {
                        H[i][i - 2] = 0.0;
                        if (i > m + 2) {
                            H[i][i - 3] = 0.0;
                        }
                    }

                    // double QR step involving rows l:n1 and columns m:n1
                    for (int k = m; k <= n1 - 1; k++) {
                        boolean notlast = (k != n1 - 1);
                        if (k != m) {
                            p = H[k][k - 1];
                            q = H[k + 1][k - 1];
                            r = (notlast ? H[k + 2][k - 1] : 0.0);
                            x = abs(p) + abs(q) + abs(r);
                            if (x != 0.0) {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        if (x == 0.0) {
                            break;
                        }
                        s = sqrt(p * p + q * q + r * r);
                        if (p < 0) {
                            s = -s;
                        }
                        if (s != 0) {
                            if (k != m) {
                                H[k][k - 1] = -s * x;
                            } else if (l != m) {
                                H[k][k - 1] = -H[k][k - 1];
                            }
                            p += s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q /= p;
                            r /= p;

                            // Row modification
                            for (int j = k; j < nn; j++) {
                                p = H[k][j] + q * H[k + 1][j];
                                if (notlast) {
                                    p += r * H[k + 2][j];
                                    H[k + 2][j] -= p * z;
                                }
                                H[k][j] -= p * x;
                                H[k + 1][j] -= p * y;
                            }

                            // Column modification
                            for (int i = 0; i <= min(n1, k + 3); i++) {
                                p = x * H[i][k] + y * H[i][k + 1];
                                if (notlast) {
                                    p += z * H[i][k + 2];
                                    H[i][k + 2] -= p * r;
                                }
                                H[i][k] -= p;
                                H[i][k + 1] -= p * q;
                            }

                            // Accumulate transformations
                            for (int i = low; i <= high; i++) {
                                p = x * V[i][k] + y * V[i][k + 1];
                                if (notlast) {
                                    p += z * V[i][k + 2];
                                    V[i][k + 2] -= p * r;
                                }
                                V[i][k] -= p;
                                V[i][k + 1] -= p * q;
                            }
                        }  // (s != 0)
                    }  // k loop
                }  // check convergence
            }  // while (n1 >= low)

            // Backsubstitute to find vectors of upper triangular form
            if (norm == 0.0) {
                return;
            }

            for (n1 = nn - 1; n1 >= 0; n1--) {
                p = d[n1];
                q = e[n1];

                // Real vector
                if (q == 0) {
                    int l = n1;
                    H[n1][n1] = 1.0;
                    for (int i = n1 - 1; i >= 0; i--) {
                        w = H[i][i] - p;
                        r = 0.0;
                        for (int j = l; j <= n1; j++) {
                            r += H[i][j] * H[j][n1];
                        }
                        if (e[i] < 0.0) {
                            z = w;
                            s = r;
                        } else {
                            l = i;
                            if (e[i] == 0.0) {
                                if (w != 0.0) {
                                    H[i][n1] = -r / w;
                                } else {
                                    H[i][n1] = -r / (eps * norm);
                                }

                                // Solve real equations
                            } else {
                                x = H[i][i + 1];
                                y = H[i + 1][i];
                                q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                                t = (x * s - z * r) / q;
                                H[i][n1] = t;
                                if (abs(x) > abs(z)) {
                                    H[i + 1][n1] = (-r - w * t) / x;
                                } else {
                                    H[i + 1][n1] = (-s - y * t) / z;
                                }
                            }

                            // Overflow control
                            t = abs(H[i][n1]);
                            if ((eps * t) * t > 1) {
                                for (int j = i; j <= n1; j++) {
                                    H[j][n1] /= t;
                                }
                            }
                        }
                    }

                    // Complex vector
                } else if (q < 0) {
                    int l = n1 - 1;

                    // Last vector component imaginary so matrix is triangular
                    if (abs(H[n1][n1 - 1]) > abs(H[n1 - 1][n1])) {
                        H[n1 - 1][n1 - 1] = q / H[n1][n1 - 1];
                        H[n1 - 1][n1] = -(H[n1][n1] - p) / H[n1][n1 - 1];
                    } else {
                        cdiv(0.0, -H[n1 - 1][n1], H[n1 - 1][n1 - 1] - p, q);
                        H[n1 - 1][n1 - 1] = cdivr;
                        H[n1 - 1][n1] = cdivi;
                    }
                    H[n1][n1 - 1] = 0.0;
                    H[n1][n1] = 1.0;
                    for (int i = n1 - 2; i >= 0; i--) {
                        double ra, sa, vr, vi;
                        ra = 0.0;
                        sa = 0.0;
                        for (int j = l; j <= n1; j++) {
                            ra += H[i][j] * H[j][n1 - 1];
                            sa += H[i][j] * H[j][n1];
                        }
                        w = H[i][i] - p;

                        if (e[i] < 0.0) {
                            z = w;
                            r = ra;
                            s = sa;
                        } else {
                            l = i;
                            if (e[i] == 0) {
                                cdiv(-ra, -sa, w, q);
                                H[i][n1 - 1] = cdivr;
                                H[i][n1] = cdivi;
                            } else {

                                // Solve complex equations
                                x = H[i][i + 1];
                                y = H[i + 1][i];
                                vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                                vi = (d[i] - p) * 2.0 * q;
                                if (vr == 0.0 && vi == 0.0) {
                                    vr = eps * norm * (abs(w) + abs(q)
                                            + abs(x) + abs(y) + abs(z));
                                }
                                cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                                H[i][n1 - 1] = cdivr;
                                H[i][n1] = cdivi;
                                if (abs(x) > (abs(z) + abs(q))) {
                                    H[i + 1][n1 - 1] = (-ra - w * H[i][n1 - 1] + q * H[i][n1]) / x;
                                    H[i + 1][n1] = (-sa - w * H[i][n1] - q * H[i][n1 - 1]) / x;
                                } else {
                                    cdiv(-r - y * H[i][n1 - 1], -s - y * H[i][n1], z, q);
                                    H[i + 1][n1 - 1] = cdivr;
                                    H[i + 1][n1] = cdivi;
                                }
                            }

                            // Overflow control
                            t = max(abs(H[i][n1 - 1]), abs(H[i][n1]));
                            if ((eps * t) * t > 1) {
                                for (int j = i; j <= n1; j++) {
                                    H[j][n1 - 1] /= t;
                                    H[j][n1] /= t;
                                }
                            }
                        }
                    }
                }
            }

            // Vectors of isolated roots
            for (int i = 0; i < nn; i++) {
                if (i < low | i > high) {
                    arraycopy(H[i], i, V[i], i, nn - i);
                }
            }

            // Back transformation to get eigenvectors of original matrix
            for (int j = nn - 1; j >= low; j--) {
                for (int i = low; i <= high; i++) {
                    z = 0.0;
                    for (int k = low; k <= min(j, high); k++) {
                        z += V[i][k] * H[k][j];
                    }
                    V[i][j] = z;
                }
            }
        }

        /* ------------------------
        Public Methods
        * ------------------------ */
        /**
         * Return the eigenvector matrix
         *
         * @return V
         */
        public EBIMatrix getV() {
            return new EBIMatrix(V, n, n);
        }

        /**
         * Return the real parts of the eigenvalues
         *
         * @return real(diag(D))
         */
        public double[] getRealEigenvalues() {
            return d;
        }

        /**
         * Return the imaginary parts of the eigenvalues
         *
         * @return imag(diag(D))
         */
        public double[] getImagEigenvalues() {
            return e;
        }

        /**
         * Return the block diagonal eigenvalue matrix
         *
         * @return D
         */
        public EBIMatrix getD() {
            EBIMatrix X = new EBIMatrix(n, n);
            double[][] D = X.getArray();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    D[i][j] = 0.0;
                }
                D[i][i] = d[i];
                if (e[i] > 0) {
                    D[i][i + 1] = e[i];
                } else if (e[i] < 0) {
                    D[i][i - 1] = e[i];
                }
            }
            return X;
        }
    }

        // ==================== Inner class LUDecomposition ====================

        /**
     * LU Decomposition.
     * <P>
     * For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n unit
     * lower triangular matrix L, an n-by-n upper triangular matrix U, and a
     * permutation vector piv of length m so that A(piv,:) = L*U. If m < n, then L
     * is m-by-m and U is m-by-n. <P>
     * The LU decompostion with pivoting always exists, even if the matrix is
     * singular, so the constructor will never fail. The primary use of the LU
     * decomposition is in the solution of square systems of simultaneous linear
     * equations. This will fail if isNonsingular() returns false.
     */
    public static class LUDecomposition implements java.io.Serializable {

        private static final long serialVersionUID = 19978681017L;

        /* ------------------------
         Class variables
         * ------------------------ */
        /**
         * Array for internal storage of decomposition.
         *
         * @serial internal array storage.
         */
        private final double[][] LU;

        /**
         * Row and column dimensions, and pivot sign.
         *
         * @serial column dimension.
         * @serial row dimension.
         * @serial pivot sign.
         */
        private final int m, n;
        private int pivsign;

        /**
         * Internal storage of pivot vector.
         *
         * @serial pivot vector.
         */
        private final int[] piv;

        /* ------------------------
         Constructor
         * ------------------------ */
        /**
         * LU Decomposition
         *
         * @param A Rectangular matrix
         */
        public LUDecomposition(EBIMatrix A) {

            // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
            LU = A.getArrayCopy();
            m = A.getRowDimension();
            n = A.getColumnDimension();
            piv = new int[m];
            for (int i = 0; i < m; i++) {
                piv[i] = i;
            }
            pivsign = 1;
            double[] LUrowi;
            double[] LUcolj = new double[m];

            // Outer loop.
            for (int j = 0; j < n; j++) {

                // Make a duplicate of the j-th column to localize references.
                for (int i = 0; i < m; i++) {
                    LUcolj[i] = LU[i][j];
                }

                // Apply previous transformations.
                for (int i = 0; i < m; i++) {
                    LUrowi = LU[i];

                    // Most of the time is spent in the following dot product.
                    int kmax = min(i, j);
                    double s = 0.0;
                    for (int k = 0; k < kmax; k++) {
                        s += LUrowi[k] * LUcolj[k];
                    }

                    LUrowi[j] = LUcolj[i] -= s;
                }

                // Find pivot and exchange if necessary.
                int p = j;
                for (int i = j + 1; i < m; i++) {
                    if (abs(LUcolj[i]) > abs(LUcolj[p])) {
                        p = i;
                    }
                }
                if (p != j) {
                    for (int k = 0; k < n; k++) {
                        double t = LU[p][k];
                        LU[p][k] = LU[j][k];
                        LU[j][k] = t;
                    }
                    int k = piv[p];
                    piv[p] = piv[j];
                    piv[j] = k;
                    pivsign = -pivsign;
                }

                // Compute multipliers.
                if (j < m && LU[j][j] != 0.0) {
                    for (int i = j + 1; i < m; i++) {
                        LU[i][j] /= LU[j][j];
                    }
                }
            }
        }

        /* ------------------------
         Temporary, experimental code.
         ------------------------ *\

         \** LU Decomposition, computed by Gaussian elimination.
         <P>
         This constructor computes L and U with the "daxpy"-based elimination
         algorithm used in LINPACK and MATLAB.  In Java, we suspect the dot-product,
         Crout algorithm will be faster.  We have temporarily included this
         constructor until timing experiments confirm this suspicion.
         <P>
         @param  A             Rectangular matrix
         @param  linpackflag   Use Gaussian elimination.  Actual value ignored.
         @return               Structure to access L, U and piv.
         *\

         public LUDecomposition (Matrix A, int linpackflag) {
         // Initialize.
         LU = A.getArrayCopy();
         m = A.getRowDimension();
         n = A.getColumnDimension();
         piv = new int[m];
         for (int i = 0; i < m; i++) {
         piv[i] = i;
         }
         pivsign = 1;
         // Main loop.
         for (int k = 0; k < n; k++) {
         // Find pivot.
         int p = k;
         for (int i = k+1; i < m; i++) {
         if (Math.abs(LU[i][k]) > Math.abs(LU[p][k])) {
         p = i;
         }
         }
         // Exchange if necessary.
         if (p != k) {
         for (int j = 0; j < n; j++) {
         double t = LU[p][j]; LU[p][j] = LU[k][j]; LU[k][j] = t;
         }
         int t = piv[p]; piv[p] = piv[k]; piv[k] = t;
         pivsign = -pivsign;
         }
         // Compute multipliers and eliminate k-th column.
         if (LU[k][k] != 0.0) {
         for (int i = k+1; i < m; i++) {
         LU[i][k] /= LU[k][k];
         for (int j = k+1; j < n; j++) {
         LU[i][j] -= LU[i][k]*LU[k][j];
         }
         }
         }
         }
         }

         \* ------------------------
         End of temporary code.
         * ------------------------ */

     /* ------------------------
         Public Methods
         * ------------------------ */
        /**
         * Is the matrix nonsingular?
         *
         * @return true if U, and hence A, is nonsingular.
         */
        public boolean isNonsingular() {
            for (int j = 0; j < n; j++) {
                if (LU[j][j] == 0) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Return lower triangular factor
         *
         * @return L
         */
        public EBIMatrix getL() {
            EBIMatrix X = new EBIMatrix(m, n);
            double[][] L = X.getArray();
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (i > j) {
                        L[i][j] = LU[i][j];
                    } else if (i == j) {
                        L[i][j] = 1.0;
                    } else {
                        L[i][j] = 0.0;
                    }
                }
            }
            return X;
        }

        /**
         * Return upper triangular factor
         *
         * @return U
         */
        public EBIMatrix getU() {
            EBIMatrix X = new EBIMatrix(n, n);
            double[][] U = X.getArray();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i <= j) {
                        U[i][j] = LU[i][j];
                    } else {
                        U[i][j] = 0.0;
                    }
                }
            }
            return X;
        }

        /**
         * Return pivot permutation vector
         *
         * @return piv
         */
        public int[] getPivot() {
            int[] p = new int[m];
            arraycopy(piv, 0, p, 0, m);
            return p;
        }

        /**
         * Return pivot permutation vector as a one-dimensional double array
         *
         * @return (double) piv
         */
        public double[] getDoublePivot() {
            double[] vals = new double[m];
            for (int i = 0; i < m; i++) {
                vals[i] = valueOf(piv[i]);
            }
            return vals;
        }

        /**
         * Determinant
         *
         * @return det(A)
         * @exception IllegalArgumentException Matrix must be square
         */
        public double det() {
            if (m != n) {
                throw new IllegalArgumentException("Matrix must be square.");
            }
            double d = valueOf(pivsign);
            for (int j = 0; j < n; j++) {
                d *= LU[j][j];
            }
            return d;
        }

        /**
         * Solve A*X = B
         *
         * @param B A Matrix with as many rows as A and any number of columns.
         * @return X so that L*U*X = B(piv,:)
         * @exception IllegalArgumentException Matrix row dimensions must agree.
         * @exception RuntimeException Matrix is singular.
         */
        public EBIMatrix solve(EBIMatrix B) {
            if (B.getRowDimension() != m) {
                throw new IllegalArgumentException("Matrix row dimensions must agree.");
            }
            if (!this.isNonsingular()) {
                throw new RuntimeException("Matrix is singular.");
            }

            // Copy right hand side with pivoting
            int nx = B.getColumnDimension();
            EBIMatrix Xmat = B.getMatrix(piv, 0, nx - 1);
            double[][] X = Xmat.getArray();

            // Solve L*Y = B(piv,:)
            for (int k = 0; k < n; k++) {
                for (int i = k + 1; i < n; i++) {
                    for (int j = 0; j < nx; j++) {
                        X[i][j] -= X[k][j] * LU[i][k];
                    }
                }
            }
            // Solve U*X = Y;
            for (int k = n - 1; k >= 0; k--) {
                for (int j = 0; j < nx; j++) {
                    X[k][j] /= LU[k][k];
                }
                for (int i = 0; i < k; i++) {
                    for (int j = 0; j < nx; j++) {
                        X[i][j] -= X[k][j] * LU[i][k];
                    }
                }
            }
            return Xmat;
        }
    }

        // ==================== Inner class QRDecomposition ====================

        /**
     * QR Decomposition.
     * <P>
     * For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
     * orthogonal matrix Q and an n-by-n upper triangular matrix R so that A = Q*R.
     * <P>
     * The QR decompostion always exists, even if the matrix does not have full
     * rank, so the constructor will never fail. The primary use of the QR
     * decomposition is in the least squares solution of nonsquare systems of
     * simultaneous linear equations. This will fail if isFullRank() returns false.
     */
    public static class QRDecomposition implements java.io.Serializable {

        private static final long serialVersionUID = 199810878617L;
        private final static ILoggingTool LOGGER
                = createLoggingTool(QRDecomposition.class);

        /* ------------------------
         Class variables
         * ------------------------ */
        /**
         * Array for internal storage of decomposition.
         *
         * @serial internal array storage.
         */
        private final double[][] QR;

        /**
         * Row and column dimensions.
         *
         * @serial column dimension.
         * @serial row dimension.
         */
        private final int m, n;

        /**
         * Array for internal storage of diagonal of R.
         *
         * @serial diagonal of R.
         */
        private final double[] Rdiag;

        /* ------------------------
         Constructor
         * ------------------------ */
        /**
         * QR Decomposition, computed by Householder reflections.
         *
         * @param A Rectangular matrix
         */
        public QRDecomposition(EBIMatrix A) {
            // Initialize.
            QR = A.getArrayCopy();
            m = A.getRowDimension();
            n = A.getColumnDimension();
            Rdiag = new double[n];

            // Main loop.
            for (int k = 0; k < n; k++) {
                // Compute 2-norm of k-th column without under/overflow.
                double nrm = 0.0;
                for (int i = k; i < m; i++) {
                    nrm = hypot(nrm, QR[i][k]);
                }

                if (nrm != 0.0) {
                    // Form k-th Householder vector.
                    if (QR[k][k] < 0) {
                        nrm = -nrm;
                    }
                    for (int i = k; i < m; i++) {
                        QR[i][k] /= nrm;
                    }
                    QR[k][k] += 1.0;

                    // Apply transformation to remaining columns.
                    for (int j = k + 1; j < n; j++) {
                        double s = 0.0;
                        for (int i = k; i < m; i++) {
                            s += QR[i][k] * QR[i][j];
                        }
                        s = -s / QR[k][k];
                        for (int i = k; i < m; i++) {
                            QR[i][j] += s * QR[i][k];
                        }
                    }
                }
                Rdiag[k] = -nrm;
            }
        }

        /* ------------------------
         Public Methods
         * ------------------------ */
        /**
         * Is the matrix full rank?
         *
         * @return true if R, and hence A, has full rank.
         */
        public boolean isFullRank() {
            for (int j = 0; j < n; j++) {
                if (Rdiag[j] == 0) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Return the Householder vectors
         *
         * @return Lower trapezoidal matrix whose columns define the reflections
         */
        public EBIMatrix getH() {
            EBIMatrix X = new EBIMatrix(m, n);
            double[][] H = X.getArray();
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (i >= j) {
                        H[i][j] = QR[i][j];
                    } else {
                        H[i][j] = 0.0;
                    }
                }
            }
            return X;
        }

        /**
         * Return the upper triangular factor
         *
         * @return R
         */
        public EBIMatrix getR() {
            EBIMatrix X = new EBIMatrix(n, n);
            double[][] R = X.getArray();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i < j) {
                        R[i][j] = QR[i][j];
                    } else if (i == j) {
                        R[i][j] = Rdiag[i];
                    } else {
                        R[i][j] = 0.0;
                    }
                }
            }
            return X;
        }

        /**
         * Generate and return the (economy-sized) orthogonal factor
         *
         * @return Q
         */
        public EBIMatrix getQ() {
            EBIMatrix X = new EBIMatrix(m, n);
            double[][] Q = X.getArray();
            for (int k = n - 1; k >= 0; k--) {
                for (int i = 0; i < m; i++) {
                    Q[i][k] = 0.0;
                }
                Q[k][k] = 1.0;
                for (int j = k; j < n; j++) {
                    if (QR[k][k] != 0) {
                        double s = 0.0;
                        for (int i = k; i < m; i++) {
                            s += QR[i][k] * Q[i][j];
                        }
                        s = -s / QR[k][k];
                        for (int i = k; i < m; i++) {
                            Q[i][j] += s * QR[i][k];
                        }
                    }
                }
            }
            return X;
        }

        /**
         * Least squares solution of A*X = B
         *
         * @param B A EBIMatrix with as many rows as A and any number of columns.
         * @return X that minimizes the two norm of Q*R*X-B.
         * @exception IllegalArgumentException EBIMatrix row dimensions must agree.
         * @exception RuntimeException EBIMatrix is rank deficient.
         */
        public EBIMatrix solve(EBIMatrix B) {
            if (B.getRowDimension() != m) {
                throw new IllegalArgumentException("Matrix row dimensions must agree.");
            }
            if (!this.isFullRank()) {
                throw new RuntimeException("Matrix is rank deficient.");
            }

            // Copy right hand side
            int nx = B.getColumnDimension();
            double[][] X = B.getArrayCopy();

            // Compute Y = transpose(Q)*B
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < nx; j++) {
                    double s = 0.0;
                    for (int i = k; i < m; i++) {
                        s += QR[i][k] * X[i][j];
                    }
                    s = -s / QR[k][k];
                    for (int i = k; i < m; i++) {
                        X[i][j] += s * QR[i][k];
                    }
                }
            }
            // Solve R*X = Y;
            for (int k = n - 1; k >= 0; k--) {
                for (int j = 0; j < nx; j++) {
                    X[k][j] /= Rdiag[k];
                }
                for (int i = 0; i < k; i++) {
                    for (int j = 0; j < nx; j++) {
                        X[i][j] -= X[k][j] * QR[i][k];
                    }
                }
            }
            return (new EBIMatrix(X, n, nx).getMatrix(0, n - 1, 0, nx - 1));
        }
    }

        // ==================== Inner class SingularValueDecomposition ====================

        /**
     * Singular Value Decomposition.
     * <P>
     * For an m-by-n matrix A with m >= n, the singular value decomposition is an
     * m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and an n-by-n
     * orthogonal matrix V so that A = U*S*V'.
     * <P>
     * The singular values, sigma[k] = S[k][k], are ordered so that sigma[0] >=
     * sigma[1] >= ... >= sigma[n-1].
     * <P>
     * The singular value decompostion always exists, so the constructor will never
     * fail. The matrix condition number and the effective numerical rank can be
     * computed from this decomposition.
     */
    public static class SingularValueDecomposition implements java.io.Serializable {

        private static final long serialVersionUID = 19981017786L;
        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(SingularValueDecomposition.class);
        /* ------------------------
         Class variables
         * ------------------------ */
        /**
         * Arrays for internal storage of U and V.
         *
         * @serial internal storage of U.
         * @serial internal storage of V.
         */
        private double[][] U, V;
        /**
         * Array for internal storage of singular values.
         *
         * @serial internal storage of singular values.
         */
        private double[] s;
        /**
         * Row and column dimensions.
         *
         * @serial row dimension.
         * @serial column dimension.
         */
        private int m, n;

        /* ------------------------
         Constructor
         * ------------------------ */
        /**
         * Construct the singular value decomposition
         *
         * @param Arg
         */
        public SingularValueDecomposition(EBIMatrix Arg) {

            // Derived from LINPACK code.
            // Initialize.
            double[][] A = Arg.getArrayCopy();
            m = Arg.getRowDimension();
            n = Arg.getColumnDimension();

            /* Apparently the failing cases are only a proper subset of (m<n), 
             so let's not throw error.  Correct fix to come later?
             if (m<n) {
             throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
             */
            int nu = min(m, n);
            s = new double[min(m + 1, n)];
            U = new double[m][nu];
            V = new double[n][n];
            double[] e = new double[n];
            double[] work = new double[m];
            boolean wantu = true;
            boolean wantv = true;

            // Reduce A to bidiagonal form, storing the diagonal elements
            // in s and the super-diagonal elements in e.
            int nct = min(m - 1, n);
            int nrt = max(0, min(n - 2, m));
            for (int k = 0; k < max(nct, nrt); k++) {
                if (k < nct) {

                    // Compute the transformation for the k-th column and
                    // place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    s[k] = 0.0;
                    for (int i = k; i < m; i++) {
                        s[k] = hypot(s[k], A[i][k]);
                    }
                    if (s[k] != 0.0) {
                        if (A[k][k] < 0.0) {
                            s[k] = -s[k];
                        }
                        for (int i = k; i < m; i++) {
                            A[i][k] /= s[k];
                        }
                        A[k][k] += 1.0;
                    }
                    s[k] = -s[k];
                }
                for (int j = k + 1; j < n; j++) {
                    if ((k < nct) && (s[k] != 0.0)) {

                        // Apply the transformation.
                        double t = 0.0;
                        for (int i = k; i < m; i++) {
                            t += A[i][k] * A[i][j];
                        }
                        t = -t / A[k][k];
                        for (int i = k; i < m; i++) {
                            A[i][j] += t * A[i][k];
                        }
                    }

                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.
                    e[j] = A[k][j];
                }
                if (wantu && (k < nct)) {

                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (int i = k; i < m; i++) {
                        U[i][k] = A[i][k];
                    }
                }
                if (k < nrt) {

                    // Compute the k-th row transformation and place the
                    // k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0.0;
                    for (int i = k + 1; i < n; i++) {
                        e[k] = hypot(e[k], e[i]);
                    }
                    if (e[k] != 0.0) {
                        if (e[k + 1] < 0.0) {
                            e[k] = -e[k];
                        }
                        for (int i = k + 1; i < n; i++) {
                            e[i] /= e[k];
                        }
                        e[k + 1] += 1.0;
                    }
                    e[k] = -e[k];
                    if ((k + 1 < m) && (e[k] != 0.0)) {

                        // Apply the transformation.
                        for (int i = k + 1; i < m; i++) {
                            work[i] = 0.0;
                        }
                        for (int j = k + 1; j < n; j++) {
                            for (int i = k + 1; i < m; i++) {
                                work[i] += e[j] * A[i][j];
                            }
                        }
                        for (int j = k + 1; j < n; j++) {
                            double t = -e[j] / e[k + 1];
                            for (int i = k + 1; i < m; i++) {
                                A[i][j] += t * work[i];
                            }
                        }
                    }
                    if (wantv) {

                        // Place the transformation in V for subsequent
                        // back multiplication.
                        for (int i = k + 1; i < n; i++) {
                            V[i][k] = e[i];
                        }
                    }
                }
            }

            // Set up the final bidiagonal matrix or order p.
            int p = min(n, m + 1);
            if (nct < n) {
                s[nct] = A[nct][nct];
            }
            if (m < p) {
                s[p - 1] = 0.0;
            }
            if (nrt + 1 < p) {
                e[nrt] = A[nrt][p - 1];
            }
            e[p - 1] = 0.0;

            // If required, generate U.
            if (wantu) {
                for (int j = nct; j < nu; j++) {
                    for (int i = 0; i < m; i++) {
                        U[i][j] = 0.0;
                    }
                    U[j][j] = 1.0;
                }
                for (int k = nct - 1; k >= 0; k--) {
                    if (s[k] != 0.0) {
                        for (int j = k + 1; j < nu; j++) {
                            double t = 0.0;
                            for (int i = k; i < m; i++) {
                                t += U[i][k] * U[i][j];
                            }
                            t = -t / U[k][k];
                            for (int i = k; i < m; i++) {
                                U[i][j] += t * U[i][k];
                            }
                        }
                        for (int i = k; i < m; i++) {
                            U[i][k] = -U[i][k];
                        }
                        U[k][k] = 1.0 + U[k][k];
                        for (int i = 0; i < k - 1; i++) {
                            U[i][k] = 0.0;
                        }
                    } else {
                        for (int i = 0; i < m; i++) {
                            U[i][k] = 0.0;
                        }
                        U[k][k] = 1.0;
                    }
                }
            }

            // If required, generate V.
            if (wantv) {
                for (int k = n - 1; k >= 0; k--) {
                    if ((k < nrt) && (e[k] != 0.0)) {
                        for (int j = k + 1; j < nu; j++) {
                            double t = 0.0;
                            for (int i = k + 1; i < n; i++) {
                                t += V[i][k] * V[i][j];
                            }
                            t = -t / V[k + 1][k];
                            for (int i = k + 1; i < n; i++) {
                                V[i][j] += t * V[i][k];
                            }
                        }
                    }
                    for (int i = 0; i < n; i++) {
                        V[i][k] = 0.0;
                    }
                    V[k][k] = 1.0;
                }
            }

            // Main iteration loop for the singular values.
            int pp = p - 1;
            int iter = 0;
            double eps = pow(2.0, -52.0);
            double tiny = pow(2.0, -966.0);
            while (p > 0) {
                int k, kase;

                // Here is where a test for too many iterations would go.
                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.
                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and
                //              s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).
                for (k = p - 2; k >= -1; k--) {
                    if (k == -1) {
                        break;
                    }
                    if (abs(e[k])
                            <= tiny + eps * (abs(s[k]) + abs(s[k + 1]))) {
                        e[k] = 0.0;
                        break;
                    }
                }
                if (k == p - 2) {
                    kase = 4;
                } else {
                    int ks;
                    for (ks = p - 1; ks >= k; ks--) {
                        if (ks == k) {
                            break;
                        }
                        double t = (ks != p ? abs(e[ks]) : 0.)
                                + (ks != k + 1 ? abs(e[ks - 1]) : 0.);
                        if (abs(s[ks]) <= tiny + eps * t) {
                            s[ks] = 0.0;
                            break;
                        }
                    }
                    if (ks == k) {
                        kase = 3;
                    } else if (ks == p - 1) {
                        kase = 1;
                    } else {
                        kase = 2;
                        k = ks;
                    }
                }
                k++;

                // Perform the task indicated by kase.
                switch (kase) {

                    // Deflate negligible s(p).
                    case 1: {
                        double f = e[p - 2];
                        e[p - 2] = 0.0;
                        for (int j = p - 2; j >= k; j--) {
                            double t = hypot(s[j], f);
                            double cs = s[j] / t;
                            double sn = f / t;
                            s[j] = t;
                            if (j != k) {
                                f = -sn * e[j - 1];
                                e[j - 1] = cs * e[j - 1];
                            }
                            if (wantv) {
                                for (int i = 0; i < n; i++) {
                                    t = cs * V[i][j] + sn * V[i][p - 1];
                                    V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                                    V[i][j] = t;
                                }
                            }
                        }
                    }
                    break;

                    // Split at negligible s(k).
                    case 2: {
                        double f = e[k - 1];
                        e[k - 1] = 0.0;
                        for (int j = k; j < p; j++) {
                            double t = hypot(s[j], f);
                            double cs = s[j] / t;
                            double sn = f / t;
                            s[j] = t;
                            f = -sn * e[j];
                            e[j] = cs * e[j];
                            if (wantu) {
                                for (int i = 0; i < m; i++) {
                                    t = cs * U[i][j] + sn * U[i][k - 1];
                                    U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                                    U[i][j] = t;
                                }
                            }
                        }
                    }
                    break;

                    // Perform one qr step.
                    case 3: {

                        // Calculate the shift.
                        double scale = max(max(max(max(abs(s[p - 1]), abs(s[p - 2])), abs(e[p - 2])),
                                abs(s[k])), abs(e[k]));
                        double sp = s[p - 1] / scale;
                        double spm1 = s[p - 2] / scale;
                        double epm1 = e[p - 2] / scale;
                        double sk = s[k] / scale;
                        double ek = e[k] / scale;
                        double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                        double c = (sp * epm1) * (sp * epm1);
                        double shift = 0.0;
                        if ((b != 0.0) | (c != 0.0)) {
                            shift = sqrt(b * b + c);
                            if (b < 0.0) {
                                shift = -shift;
                            }
                            shift = c / (b + shift);
                        }
                        double f = (sk + sp) * (sk - sp) + shift;
                        double g = sk * ek;

                        // Chase zeros.
                        for (int j = k; j < p - 1; j++) {
                            double t = hypot(f, g);
                            double cs = f / t;
                            double sn = g / t;
                            if (j != k) {
                                e[j - 1] = t;
                            }
                            f = cs * s[j] + sn * e[j];
                            e[j] = cs * e[j] - sn * s[j];
                            g = sn * s[j + 1];
                            s[j + 1] = cs * s[j + 1];
                            if (wantv) {
                                for (int i = 0; i < n; i++) {
                                    t = cs * V[i][j] + sn * V[i][j + 1];
                                    V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                                    V[i][j] = t;
                                }
                            }
                            t = hypot(f, g);
                            cs = f / t;
                            sn = g / t;
                            s[j] = t;
                            f = cs * e[j] + sn * s[j + 1];
                            s[j + 1] = -sn * e[j] + cs * s[j + 1];
                            g = sn * e[j + 1];
                            e[j + 1] = cs * e[j + 1];
                            if (wantu && (j < m - 1)) {
                                for (int i = 0; i < m; i++) {
                                    t = cs * U[i][j] + sn * U[i][j + 1];
                                    U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                                    U[i][j] = t;
                                }
                            }
                        }
                        e[p - 2] = f;
                        iter += 1;
                    }
                    break;

                    // Convergence.
                    case 4: {

                        // Make the singular values positive.
                        if (s[k] <= 0.0) {
                            s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                            if (wantv) {
                                for (int i = 0; i <= pp; i++) {
                                    V[i][k] = -V[i][k];
                                }
                            }
                        }

                        // Order the singular values.
                        while (k < pp) {
                            if (s[k] >= s[k + 1]) {
                                break;
                            }
                            double t = s[k];
                            s[k] = s[k + 1];
                            s[k + 1] = t;
                            if (wantv && (k < n - 1)) {
                                for (int i = 0; i < n; i++) {
                                    t = V[i][k + 1];
                                    V[i][k + 1] = V[i][k];
                                    V[i][k] = t;
                                }
                            }
                            if (wantu && (k < m - 1)) {
                                for (int i = 0; i < m; i++) {
                                    t = U[i][k + 1];
                                    U[i][k + 1] = U[i][k];
                                    U[i][k] = t;
                                }
                            }
                            k++;
                        }
                        iter = 0;
                        p--;
                    }
                    break;
                }
            }
        }

        /* ------------------------
         Public Methods
         * ------------------------ */
        /**
         * Return the left singular vectors
         *
         * @return U
         */
        public EBIMatrix getU() {
            return new EBIMatrix(U, m, min(m + 1, n));
        }

        /**
         * Return the right singular vectors
         *
         * @return V
         */
        public EBIMatrix getV() {
            return new EBIMatrix(V, n, n);
        }

        /**
         * Return the one-dimensional array of singular values
         *
         * @return diagonal of S.
         */
        public double[] getSingularValues() {
            return s;
        }

        /**
         * Return the diagonal matrix of singular values
         *
         * @return S
         */
        public EBIMatrix getS() {
            EBIMatrix X = new EBIMatrix(n, n);
            double[][] S = X.getArray();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    S[i][j] = 0.0;
                }
                S[i][i] = this.s[i];
            }
            return X;
        }

        /**
         * Two norm
         *
         * @return max(S)
         */
        public double norm2() {
            return s[0];
        }

        /**
         * Two norm condition number
         *
         * @return max(S)/min(S)
         */
        public double cond() {
            return s[0] / s[min(m, n) - 1];
        }

        /**
         * Effective numerical matrix rank
         *
         * @return Number of nonnegligible singular values.
         */
        public int rank() {
            double eps = pow(2.0, -52.0);
            double tol = max(m, n) * s[0] * eps;
            int r = 0;
            for (int i = 0; i < s.length; i++) {
                if (s[i] > tol) {
                    r++;
                }
            }
            return r;
        }
    }

}
