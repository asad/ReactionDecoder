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
package uk.ac.ebi.reactionblast.mechanism;

import java.io.Serializable;
import static java.lang.System.getProperty;

import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class MappingSolution implements Serializable {

    private static final long serialVersionUID = 1678787866L;

    private final IMappingAlgorithm algorithmID;
    private final double bondEnergySum;
    private final double energyDelta;
    private final int totalBondChanges;
    private final int totalFragmentChanges;
    private final int totalStereoChanges;
    private final int smallestFragmentCount;
    private final int totalCarbonBondChanges;
    private final IReaction reaction;
    private final Reactor reactor;
    private final int totalChanges;
    private boolean chosen;
    private final BondChangeCalculator bondChangeCalculator;
    private boolean generate3D;
    private boolean generate2D;

    /**
     *
     * @param bondChangeCalculator
     * @param ma
     * @param reactor
     * @param localScore
     * @param bondEnergyChange
     * @param totalCarbonBondChanges
     * @param totalBondChanges
     * @param totalFragmentChanges
     * @param totalStereoChanges
     * @param smallestFragmentCount
     * @param reaction
     * @param energyDelta
     */
    public MappingSolution(BondChangeCalculator bondChangeCalculator, IMappingAlgorithm ma, IReaction reaction, Reactor reactor, double bondEnergyChange, int totalCarbonBondChanges, int totalBondChanges, int totalFragmentChanges, int totalStereoChanges, int smallestFragmentCount, int localScore, double energyDelta) {
        this.algorithmID = ma;
        this.bondEnergySum = bondEnergyChange;
        this.totalBondChanges = totalBondChanges;
        this.totalFragmentChanges = totalFragmentChanges;
        this.totalStereoChanges = totalStereoChanges;
        this.smallestFragmentCount = smallestFragmentCount;
        this.reaction = reaction;
        this.reactor = reactor;
        this.totalChanges = localScore;
        this.bondChangeCalculator = bondChangeCalculator;
        this.totalCarbonBondChanges = totalCarbonBondChanges;
        this.chosen = false;
        this.generate2D = false;
        this.generate3D = false;
        this.energyDelta = energyDelta;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        String property = getProperty("line.separator");
        sb.append(property);
        sb.append("-----------------------------------");
        sb.append(property);
        sb.append("Chosen Algorithm is= ").append(this.getAlgorithmID().description());
        sb.append(", ");
        sb.append(property).append("Scores=" + "(Chaos Delta:");
        sb.append(this.getTotalFragmentChanges()).append(", Sigma:");
        sb.append(this.getTotalBondChanges()).append(", Energy:");
        sb.append(this.getBondEnergySum()).append(")");
        sb.append(property);
        sb.append("-----------------------------------");
        sb.append(property);
        sb.append("MappingSolution{" + "algorithmID=").append(algorithmID).append(", bondEnergyChange=").append(bondEnergySum).append(", totalBondChanges=").append(totalBondChanges).append(", totalFragmentChanges=").append(totalFragmentChanges).append(", totalStereoChanges=").append(totalStereoChanges).append(", smallestFragmentCount=").append(smallestFragmentCount) //                .append(", reaction=").append(reaction)
                //                .append(", reactor=").append(reactor)
                .append(", totalChanges=").append(totalChanges).append(", chosen=").append(chosen) //                .append(", bondChangeCalculator=").append(bondChangeCalculator)
                .append(", generate3D=").append(generate3D).append(", generate2D=").append(generate2D).append('}');
        sb.append(property);
        return sb.toString();
    }

    /**
     * @return the algorithmID
     */
    public IMappingAlgorithm getAlgorithmID() {
        return algorithmID;
    }

    /**
     * @return the bondEnergySum
     */
    public double getBondEnergySum() {
        return bondEnergySum;
    }

    /**
     * @return the totalBondChanges
     */
    public int getTotalBondChanges() {
        return totalBondChanges;
    }

    /**
     * @return the totalFragmentChanges
     */
    public int getTotalFragmentChanges() {
        return totalFragmentChanges;
    }

    /**
     * @return the totalStereoChanges
     */
    public int getTotalStereoChanges() {
        return totalStereoChanges;
    }

    /**
     * @return the smallestFragmentCount
     */
    public int getSmallestFragmentCount() {
        return smallestFragmentCount;
    }

    /**
     * @return the reaction
     */
    public IReaction getReaction() {
        return reaction;
    }

    /**
     * @return the reactor
     */
    public Reactor getReactor() {
        return reactor;
    }

    /**
     * @return the totalChanges
     */
    public int getTotalChanges() {
        return totalChanges;
    }

    /**
     * @return the chosen
     */
    public boolean isChosen() {
        return chosen;
    }

    /**
     * @return the bondChangeCalculator
     */
    public BondChangeCalculator getBondChangeCalculator() {
        return bondChangeCalculator;
    }

    /**
     * @param chosen the chosen to set
     */
    void setChosen(boolean chosen) {
        this.chosen = chosen;
    }

    /**
     * @return the generate3D
     */
    boolean isGenerate3D() {
        return generate3D;
    }

    /**
     * @param generate3D the generate3D to set
     */
    void setGenerate3D(boolean generate3D) {
        this.generate3D = generate3D;
    }

    /**
     * @return the generate2D
     */
    boolean isGenerate2D() {
        return generate2D;
    }

    /**
     * @param generate2D the generate2D to set
     */
    void setGenerate2D(boolean generate2D) {
        this.generate2D = generate2D;
    }

    /**
     *
     * @return
     */
    public int getTotalCarbonBondChanges() {
        return this.totalCarbonBondChanges;
    }

    /**
     * @return the energyDelta
     */
    public double getEnergyDelta() {
        return energyDelta;
    }
}
