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
package uk.ac.ebi.reactionblast.tools;

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
import static java.lang.System.out;
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
import uk.ac.ebi.reactionblast.tools.matrix.CholeskyDecomposition;
import uk.ac.ebi.reactionblast.tools.matrix.EigenvalueDecomposition;
import uk.ac.ebi.reactionblast.tools.matrix.LUDecomposition;
import static uk.ac.ebi.reactionblast.tools.matrix.Maths.hypot;
import uk.ac.ebi.reactionblast.tools.matrix.QRDecomposition;
import uk.ac.ebi.reactionblast.tools.matrix.SingularValueDecomposition;

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
    public synchronized void initMatrix(double v) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = v;
            }
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
    public synchronized double getValue(int i, int j) {

        double val = -1.0d;
        if (i <= rows && j <= columns) {
            val = matrix[i][j];
        } else {

            LOGGER.debug("Error: Array of out bound");
        }
        return val;
    }

    /**
     * Make a deep duplicate of a matrix
     *
     * @return
     */
    public synchronized EBIMatrix duplicate() {
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
    public synchronized Object clone() throws CloneNotSupportedException {
        return this.duplicate();
    }

    /**
     * data duplicate the internal two-dimensional array.
     *
     * @return Two-dimensional array duplicate of matrix elements.
     */
    public synchronized double[][] getArrayCopy() {
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
    public synchronized double[] getColumnPackedCopy() {
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
    public synchronized double[] getRowPackedCopy() {
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
    public synchronized double[][] getArray() {
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
    public synchronized EBIMatrix getMatrix(int rowStart, int rowEnd, int colStart, int colEnd) {
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
    public synchronized EBIMatrix transpose() {
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
    public synchronized boolean setValue(int row, int col, double value) {

        double val = value;
        boolean flag = false;

        if (row <= rows && col <= columns) {
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
    public synchronized int getRowDimension() {
        return this.rows;
    }

    /**
     * Get column dimension.
     *
     * @return columns, the number of columns.
     */
    public synchronized int getColumnDimension() {
        return this.columns;
    }

    /**
     *
     * @return
     */
    public synchronized List<Double> getDiagonalElements() {

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

            out.println("Row =/= Columns");
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
    public synchronized EBIMatrix getMatrix(int[] r, int[] c) {
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
    public synchronized EBIMatrix getMatrix(int rowStart, int rowEnd, int[] c) {
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
    public synchronized EBIMatrix getMatrix(int[] r, int colStart, int colEnd) {
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
    public synchronized void set(int i, int j, double s) {
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
    public synchronized void setMatrix(int rowStart, int rowEnd, int colStart, int colEnd, EBIMatrix X) {
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
    public synchronized void setMatrix(int[] r, int[] c, EBIMatrix X) {
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
    public synchronized void setMatrix(int[] r, int colStart, int colEnd, EBIMatrix X) {
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
    public synchronized void setMatrix(int rowStart, int rowEnd, int[] c, EBIMatrix X) {
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
    public synchronized boolean is_element_max_in_column(int iPos, int jPos) {

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
    public synchronized boolean is_element_min_in_column(int iPos, int jPos) {

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
    public synchronized void reSizeMatrix(int RowSize, int ColSize) {

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
    public synchronized boolean is_element_max_in_row(int iPos, int jPos) {

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
    public synchronized boolean is_element_min_in_row(int iPos, int jPos) {

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
    public synchronized void swapColumns(int coloumn1, int coloumn2) {

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
    public synchronized void swapRows(int row1, int row2) throws CDKException {

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
    public synchronized void pivot(int row, int col) {

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
    public synchronized EBIMatrix normalize(EBIMatrix S) {
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
                out.println("Warning(orthonormalize):" + (p + 1) + ". Vector has length null");
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
    public synchronized EBIMatrix mul(double a) {
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
    public synchronized List<Double> mul(List<Double> a) {
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
    public synchronized EBIMatrix mul(EBIMatrix b) {
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
    private synchronized void checkMatrixDimensions(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayTimes(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayTimesEquals(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayRightDivide(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayRightDivideEquals(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayLeftDivide(EBIMatrix B) {
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
    public synchronized EBIMatrix arrayLeftDivideEquals(EBIMatrix B) {
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
    public synchronized EBIMatrix times(double s) {
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
    public synchronized EBIMatrix timesEquals(double s) {
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
    public synchronized EBIMatrix times(EBIMatrix B) {
        if (B.getRowDimension() != columns) {
            throw new IllegalArgumentException("EBIMatrix inner dimensions must agree.");
        }
        EBIMatrix X = new EBIMatrix(rows, B.getColumnDimension());
        double[][] C = X.getArray();
        double[] Bcolj = new double[columns];
        for (int j = 0; j < B.getColumnDimension(); j++) {
            for (int k = 0; k < columns; k++) {
                Bcolj[k] = B.matrix[k][j];
            }
            for (int i = 0; i < rows; i++) {
                double[] Arowi = matrix[i];
                double s = 0;
                for (int k = 0; k < columns; k++) {
                    s += Arowi[k] * Bcolj[k];
                }
                C[i][j] = s;
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
    public synchronized LUDecomposition lu() {
        return new LUDecomposition(this);
    }

    /**
     * QR Decomposition
     *
     * @return QRDecomposition
     * @see QRDecomposition
     */
    public synchronized QRDecomposition qr() {
        return new QRDecomposition(this);
    }

    /**
     * dataholesky Decomposition
     *
     * @return dataholeskyDecomposition
     * @see CholeskyDecomposition
     */
    public synchronized CholeskyDecomposition chol() {
        return new CholeskyDecomposition(this);
    }

    /**
     * Singular Value Decomposition
     *
     * @return SingularValueDecomposition
     * @see SingularValueDecomposition
     */
    public synchronized SingularValueDecomposition svd() {
        return new SingularValueDecomposition(this);
    }

    /**
     * Eigenvalue Decomposition
     *
     * @return EigenvalueDecomposition
     * @see EigenvalueDecomposition
     */
    public synchronized EigenvalueDecomposition eig() {
        return new EigenvalueDecomposition(this);
    }

    /**
     * Solve matrix*result = B
     *
     * @param B right hand side
     * @return solution if matrix is square, least squares solution otherwise
     */
    public synchronized EBIMatrix solve(EBIMatrix B) {
        return (rows == columns ? (new LUDecomposition(this)).solve(B) : (new QRDecomposition(this)).solve(B));
    }

    /**
     * Solve result*matrix = B, which is also matrix'*result' = B'
     *
     * @param B right hand side
     * @return solution if matrix is square, least squares solution otherwise.
     */
    public synchronized EBIMatrix solveTranspose(EBIMatrix B) {
        return transpose().solve(B.transpose());
    }

    /**
     * EBIMatrix inverse or pseudoinverse
     *
     * @return inverse(matrix) if matrix is square, pseudoinverse otherwise.
     */
    public synchronized EBIMatrix inverse() {
        return solve(identity(rows, rows));
    }

    /**
     * EBIMatrix determinant
     *
     * @return determinant
     */
    public synchronized double det() {
        return new LUDecomposition(this).det();
    }

    /**
     * EBIMatrix rank
     *
     * @return effective numerical rank, obtained from SVD.
     */
    public synchronized int rank() {
        return new SingularValueDecomposition(this).rank();
    }

    /**
     * EBIMatrix condition (2 norm)
     *
     * @return ratio of largest to smallest singular value.
     */
    public synchronized double cond() {
        return new SingularValueDecomposition(this).cond();
    }

    /**
     * EBIMatrix trace.
     *
     * @return sum of the diagonal elements.
     */
    public synchronized double trace() {
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
    public synchronized EBIMatrix diagonalize(int nrot) {
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
        out.println("Too many iterations in routine JACOBI");
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
    public synchronized EBIMatrix orthonormalize(EBIMatrix S) {
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
                out.println("Warning(orthonormalize):" + (p + 1) + ". Vector has length null");
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
    public synchronized void print(NumberFormat format, int width) {
        print(new PrintWriter(out, true), format, width);
    }

    /**
     * Print the matrix to stdout. Line the elements up in columns with a
     * Fortran-like 'Fw.d' style format.
     *
     * @param w dataolumn width.
     * @param d Number of digits after the decimal.
     */
    public synchronized void print(int w, int d) {
        print(new PrintWriter(out, true), w, d);
    }

    /**
     * Print the matrix to the output stream. Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
     *
     * @param output Output stream.
     * @param w dataolumn width.
     * @param d Number of digits after the decimal.
     */
    public synchronized void print(PrintWriter output, int w, int d) {
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
    public synchronized void print(PrintWriter output, NumberFormat format, int width) {
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
    public synchronized String toString() {
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
    public synchronized double norm1() {
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
    public synchronized double norm2() {
        return (new SingularValueDecomposition(this).norm2());
    }

    /**
     * Infinity norm
     *
     * @return maximum row sum.
     */
    public synchronized double normInf() {
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
    public synchronized double normF() {
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
    public synchronized EBIMatrix uminus() {
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
    public synchronized EBIMatrix plus(EBIMatrix B) {
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
    public synchronized EBIMatrix plusEquals(EBIMatrix B) {
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
    public synchronized EBIMatrix minus(EBIMatrix B) {
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
    public synchronized EBIMatrix minusEquals(EBIMatrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] -= B.matrix[i][j];
            }
        }
        return this;
    }

}
