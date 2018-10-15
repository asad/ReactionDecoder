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
package uk.ac.ebi.centres.cdk;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.calculator.ThreeDimensionalSignCalculator;

/**
 * @author John May
 */
public class CDK3DSignCalculator extends ThreeDimensionalSignCalculator<IAtom> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CDK3DSignCalculator.class);

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public double getX(IAtom atom) {
        return atom.getPoint3d().x;
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public double getY(IAtom atom) {
        return atom.getPoint3d().y;
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public double getZ(IAtom atom) {
        return atom.getPoint3d().z;
    }
}
