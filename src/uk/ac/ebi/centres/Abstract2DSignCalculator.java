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
package uk.ac.ebi.centres;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.signum;

/**
 * @author John May
 * @param <A>
 */
public abstract class Abstract2DSignCalculator<A> implements SignCalculator<A> {

    @Override
    public int getSign(Ligand<A> a1, Ligand<A> a2, Ligand<A> a3, Ligand<A> a4) {

        // unspecified
        if (a1.getDepth() == 0 && a2.getDepth() == 0
                && a3.getDepth() == 0 && a4.getDepth() == 0) {
            return 0;
        }

        double[][] matrix = new double[][]{{getX(a1.getAtom()), getY(a1.getAtom()), 1, a1.getDepth()},
            {getX(a2.getAtom()), getY(a2.getAtom()), 1, a2.getDepth()},
            {getX(a3.getAtom()), getY(a3.getAtom()), 1, a3.getDepth()},
            {getX(a4.getAtom()), getY(a4.getAtom()), 1, a4.getDepth()},};


        return (int) abs(determinant(matrix));


    }

    @Override
    public int getSign(A a1, A a2, A a3) {
        double[][] matrix = new double[][]{{getX(a1), getY(a1), 1},
            {getX(a2), getY(a2), 1},
            {getX(a3), getY(a3), 1}};
        return (int) signum(determinant(matrix));
    }

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

    /**
     * Copy pasta - http://www.roseindia.net/tutorial/java/core/finddeterminant.html
     */
    private double determinant(double[][] arr) {
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
