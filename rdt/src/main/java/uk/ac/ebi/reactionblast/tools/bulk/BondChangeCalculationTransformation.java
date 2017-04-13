/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

package uk.ac.ebi.reactionblast.tools.bulk;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;
import static uk.ac.ebi.reactionblast.interfaces.ITransformation.TargetType.REACTION;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;

/**
 *
 * @author asad
 */
public class BondChangeCalculationTransformation implements ITransformation<IReaction> {

    private static final Logger LOG = getLogger(BondChangeCalculationTransformation.class.getName());

    /**
     *
     * @return
     */
    @Override
    public TargetType getTargetType() {
        return REACTION;
    }

    /**
     *
     * @param reaction
     * @return
     */
    @Override
    public IReaction transform(IReaction reaction) {
        try {
            boolean generate2D = false;
            boolean generate3D = false;
            BondChangeCalculator bcc =
                    new BondChangeCalculator(reaction, generate2D, generate3D);
            return bcc.getReactionWithCompressUnChangedHydrogens();
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
}
