/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.calculator;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import uk.ac.ebi.centres.SignCalculator;

/**
 * @author John May
 * @param <A>
 */
public abstract class AbstractSignCalculator<A> implements SignCalculator<A> {

    private static final int x = 0;
    private static final int y = 1;
    private static final int z = 2;

    /**
     *
     * @param atom
     * @return
     */
    public abstract double getX(A atom);

    /**
     *
     * @param atom
     * @return
     */
    public abstract double getY(A atom);

    double[] crossproduct(double[] v1, double[] v2) {
        return new double[]{(v1[y] * v2[z]) - (v2[y] * v1[y]),
            (v1[z] * v2[x]) - (v2[z] * v1[x]),
            (v1[x] * v2[y]) - (v2[x] * v1[y])};
    }

    double magnitude(double[] vector) {
        return sqrt(vector[x] * vector[x]
                + vector[y] * vector[y]
                + vector[z] * vector[z]);

    }

    abstract double[] toVector(A base, A atom);

    /**
     * Copy pasta - http://www.roseindia.net/tutorial/java/core/finddeterminant.html
     */
    double determinant(double[][] arr) {
        double result = 0;
        if (arr.length == 1) {
            result = arr[0][0];
            return result;
        }
        if (arr.length == 2) {
            result = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
            return result;
        }
        for (int i = 0; i < arr[0].length; i++) {
            double temp[][] = new double[arr.length - 1][arr[0].length - 1];

            for (int j = 1; j < arr.length; j++) {
                for (int k = 0; k < arr[0].length; k++) {

                    if (k < i) {
                        temp[j - 1][k] = arr[j][k];
                    } else if (k > i) {
                        temp[j - 1][k - 1] = arr[j][k];
                    }
                }
            }
            result += arr[0][i] * pow(-1, i) * determinant(temp);
        }
        return result;
    }
}
