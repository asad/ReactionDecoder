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
package uk.ac.ebi.reactionblast.tools.descriptors;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class LipinskiRuleOfFiveDescriptor {

    /**
     * *********Default Lipinski Rule***********
     */
    private static double xlogPvalueLipinski = 5.0;
    private static int acceptorsLipinski = 10;
    private static int donorsLipinski = 5;
    private static double mwvalueLipinski = 500.0;
    private static double rotatablebondsLipinski = 10.0;
     private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(LipinskiRuleOfFiveDescriptor.class);
    /**
     * *********Value from the User***********
     */
    private double xlogPvalue = -1.0;
    private int acceptors = -1;
    private int donors = -1;
    private double mwvalue = -1.0;
    private double rotatablebonds = -1.0;

    /**
     * *******************************************
     */
    public LipinskiRuleOfFiveDescriptor() {
    }

    /**
     *
     * @param xlogPvalue overwrite default getXlogP cutoff for lipinski rule
     * @param acceptors overwrite default Hydrogen bond acceptors cutoff for
     * lipinski rule
     * @param donors overwrite default Hydrogen bond donors cutoff for lipinski
     * rule
     * @param mwvalue overwrite default molecular weight cutoff for lipinski
     * rule
     * @param rotatablebonds overwrite default rotatablebonds cutoff for
     * lipinski rule
     */
    public void resetRuleOfFiveDescriptor(double xlogPvalue, int acceptors, int donors, double mwvalue, double rotatablebonds) {
        xlogPvalueLipinski = xlogPvalue;
        acceptorsLipinski = acceptors;
        donorsLipinski = donors;
        mwvalueLipinski = mwvalue;
        rotatablebondsLipinski = rotatablebonds;
    }

    /**
     *
     * @param molecule
     * @return total number of lipinski rule failures (Max. is five)
     * @throws Exception *
     */
    public int calculate(IAtomContainer molecule) throws Exception {

        int lipinskifailures = 0;

        CDKMolecularDescriptor cmd = new CDKMolecularDescriptor(molecule);
        this.xlogPvalue = cmd.getXlogP(true);
        this.acceptors = cmd.getHBondAcceptors(true);
        this.donors = cmd.getHBondDoners(true);
        this.mwvalue = cmd.getMolecularWeight();
        this.rotatablebonds = cmd.getRotatableBondsCountDescriptor(true, false);

        if (xlogPvalue > xlogPvalueLipinski) {
            lipinskifailures += 1;
        }
        if (acceptors > acceptorsLipinski) {
            lipinskifailures += 1;
        }
        if (donors > donorsLipinski) {
            lipinskifailures += 1;
        }
        if (mwvalue > mwvalueLipinski) {
            lipinskifailures += 1;
        }
        if (rotatablebonds > rotatablebondsLipinski) {
            lipinskifailures += 1;
        }

        return lipinskifailures;

    }

    /**
     *
     * @param xlogPvalue
     * @param acceptors
     * @param donors
     * @param mwvalue
     * @param rotatablebonds
     * @return total number of lipinski rule failures (Max. is five)
     */
    public int calculate(double xlogPvalue, int acceptors, int donors, double mwvalue, double rotatablebonds) {

        int lipinskifailures = 0;

        if (xlogPvalue > xlogPvalueLipinski) {
            lipinskifailures += 1;
        }
        if (acceptors > acceptorsLipinski) {
            lipinskifailures += 1;
        }
        if (donors > donorsLipinski) {
            lipinskifailures += 1;
        }
        if (mwvalue > mwvalueLipinski) {
            lipinskifailures += 1;
        }
        if (rotatablebonds > rotatablebondsLipinski) {
            lipinskifailures += 1;
        }

        return lipinskifailures;

    }
}
