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

import static java.lang.Math.signum;
import uk.ac.ebi.centres.Ligand;

/**
 * @author John May
 * @param <A>
 */
public abstract class ThreeDimensionalSignCalculator<A>
        extends AbstractSignCalculator<A> {

    /**
     *
     * @param atom
     * @return
     */
    public abstract double getZ(A atom);

    @Override
    double[] toVector(A base, A atom) {
        return new double[]{getX(atom) - getX(base),
            getY(atom) - getY(base),
            getZ(atom) - getZ(base)};
    }

    @Override
    public int getSign(Ligand<A> a1, Ligand<A> a2, Ligand<A> a3, Ligand<A> a4) {

        // unspecified
        if (a1.getDepth() == 0 && a2.getDepth() == 0
                && a3.getDepth() == 0 && a4.getDepth() == 0) {
            return 0;
        }

        double[][] matrix = new double[][]{{getX(a1.getAtom()), getY(a1.getAtom()), getZ(a1.getAtom()), 1},
            {getX(a2.getAtom()), getY(a2.getAtom()), getZ(a2.getAtom()), 1},
            {getX(a3.getAtom()), getY(a3.getAtom()), getZ(a3.getAtom()), 1},
            {getX(a4.getAtom()), getY(a4.getAtom()), getZ(a4.getAtom()), 1},};

        return (int) signum(determinant(matrix));


    }

    @Override
    public int getSign(A a1, A a2, A a3) {
        double[][] matrix = new double[][]{{getX(a1), getY(a1), getZ(a1)},
            {getX(a2), getY(a2), getZ(a2)},
            {getX(a3), getY(a3), getZ(a3)}};


        // checking the size of the sign doesn't work for 3D as it does for 2D
        // instead we used the magnitude of the cross-product.
        double magnitude = magnitude(crossproduct(toVector(a2, a1), toVector(a2, a3)));
        return magnitude < 0.2 ? 0 : (int) signum(determinant(matrix));
    }
}
