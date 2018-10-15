package uk.ac.ebi.reactionblast.tools.matrix;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.System.arraycopy;
import uk.ac.ebi.reactionblast.tools.EBIMatrix;
import static uk.ac.ebi.reactionblast.tools.matrix.Maths.hypot;

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
public class EigenvalueDecomposition implements java.io.Serializable {

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
        for (int j = 0; (j < n) & issymmetric; j++) {
            for (int i = 0; (i < n) & issymmetric; i++) {
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
                            if (vr == 0.0 & vi == 0.0) {
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
