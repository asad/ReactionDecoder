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

import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import static org.openscience.cdk.interfaces.IAtomType.Hybridization.SP3;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.DOWN_INVERTED;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_INVERTED;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_OR_DOWN;
import static org.openscience.cdk.interfaces.IBond.Stereo.UP_OR_DOWN_INVERTED;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import uk.ac.ebi.centres.Centre;
import uk.ac.ebi.centres.CentreProvider;
import uk.ac.ebi.centres.ConnectionTable;
import uk.ac.ebi.centres.DescriptorManager;
import uk.ac.ebi.centres.graph.ConnectionTableDigraph;
import uk.ac.ebi.centres.ligand.PlanarCentre;
import uk.ac.ebi.centres.ligand.TetrahedralCentre;

/**
 * @author John May
 */
public class CDKCentreProvider implements CentreProvider<IAtom> {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CDKCentreProvider.class);

    private final IAtomContainer container;
    private final ConnectionTable<IAtom> table;
    private IAtomContainer cyclicFragments;

    /**
     *
     * @param container
     */
    public CDKCentreProvider(IAtomContainer container) {
        this.container = container;
        this.table = new CDKConnectionTable(container);
    }

    /**
     *
     * @return
     */
    @Override
    public Integer getAtomCount() {
        return container.getAtomCount();
    }

    /**
     *
     * @param manager
     * @return
     */
    @Override
    public Collection<Centre<IAtom>> getCentres(DescriptorManager<IAtom> manager) {

        List<Centre<IAtom>> centres = new ArrayList<>(container.getAtomCount());

        // tetrahedral centres
        for (IAtom atom : container.atoms()) {

            // might need refinement
            if (SP3.equals(atom.getHybridization())
                    && container.getConnectedBondsCount(atom) > 2
                    && atom.getFormalNeighbourCount() == 4
                    && hasStereoBonds(container, atom)) {
                TetrahedralCentre<IAtom> centre = new TetrahedralCentre<>(manager.getDescriptor(atom), atom);
                centre.setProvider(new ConnectionTableDigraph<>(centre, manager, table));
                centres.add(centre);
            }
        }

        // planar centres
        for (IBond bond : container.bonds()) {
            if (DOUBLE.equals(bond.getOrder())
                    && container.getConnectedBondsCount(bond.getAtom(0)) > 1
                    && container.getConnectedBondsCount(bond.getAtom(1)) > 1
                    && bond.getFlag(ISAROMATIC) == FALSE
                    && onlyConnectedToSingleBonds(bond, container)
                    && !getCyclicFragments().contains(bond)
                    && !hasVariableBond(container, bond.getAtom(0))
                    && !hasVariableBond(container, bond.getAtom(1))) {
                PlanarCentre<IAtom> centre = new PlanarCentre<>(bond.getAtom(0), bond.getAtom(1),
                        manager.getDescriptor(bond.getAtom(0), bond.getAtom(1)));
                centre.setProvider(new ConnectionTableDigraph<>(centre, manager, table));
                centres.add(centre);

            }
        }

        return centres;

    }

    /**
     * stops tandem double bonds being provided C=C=C \
     *
     * being provided. see. unit test of 2-iminoethen-1-ol (testIminoethenol)
     *
     * @param bond
     * @param container
     * @return
     */
    private boolean onlyConnectedToSingleBonds(IBond bond, IAtomContainer container) {
        return onlyConnectedToSingleBonds(bond, bond.getAtom(0), container)
                && onlyConnectedToSingleBonds(bond, bond.getAtom(1), container);
    }

    private boolean onlyConnectedToSingleBonds(IBond bond, IAtom atom, IAtomContainer container) {
        for (IBond connected : container.getConnectedBondsList(atom)) {
            if (!SINGLE.equals(connected.getOrder()) && !connected.equals(bond)) {
                return FALSE;
            }
        }
        return TRUE;
    }

    private IAtomContainer getCyclicFragments() {
        if (cyclicFragments == null) {
            cyclicFragments = new SpanningTree(container).getCyclicFragmentsContainer();
        }
        return cyclicFragments;
    }

    private boolean hasVariableBond(IAtomContainer container, IAtom atom) {
        for (IBond bond : container.getConnectedBondsList(atom)) {
            IBond.Stereo stereo = bond.getStereo();
            if (UP_OR_DOWN.equals(stereo)
                    || UP_OR_DOWN_INVERTED.equals(stereo)) {
                return TRUE;
            }
        }
        return FALSE;
    }

    private boolean hasStereoBonds(IAtomContainer container, IAtom atom) {
        for (IBond bond : container.getConnectedBondsList(atom)) {
            IBond.Stereo stereo = bond.getStereo();
            if (UP.equals(stereo)
                    || DOWN.equals(stereo)
                    || UP_INVERTED.equals(stereo)
                    || DOWN_INVERTED.equals(stereo)) {
                return TRUE;
            }
        }
        return FALSE;
    }
}
