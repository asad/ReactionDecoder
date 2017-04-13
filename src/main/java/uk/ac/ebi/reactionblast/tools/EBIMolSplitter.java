/*
 * Copyright (C) 2007-2017 Syed Asad Rahman <asad @ ebi.ac.uk>.
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

 /*
 * EBIMolSplitter.java
 *
 * Created on 24 January 2007, 14:07
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast.tools;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.VISITED;
import static org.openscience.cdk.graph.PathTools.breadthFirstSearch;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IElectronContainer;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class EBIMolSplitter {

    private static final Logger LOG = getLogger(EBIMolSplitter.class.getName());

    /**
     * Check whether a set of atoms in an atomcontainer is connected
     *
     * @param atomContainer The GraphAtomContainer to be check for connectedness
     * @return true if the GraphAtomContainer is connected
     */
    public static boolean isConnected(IAtomContainer atomContainer) {
        boolean flag = false;

        IAtomContainer ac = atomContainer.getBuilder().newInstance(IAtomContainer.class);
        IAtom atom = null;
        IAtomContainer molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
        List<IAtom> sphere = new ArrayList<>();
        for (int f = 0; f < atomContainer.getAtomCount(); f++) {
            atom = atomContainer.getAtom(f);
            atom.setFlag(VISITED, false);
            ac.addAtom(atomContainer.getAtom(f));
        }

        Iterator<IBond> bonds = atomContainer.bonds().iterator();
        while (bonds.hasNext()) {
            IBond bond = bonds.next();
            bond.setFlag(VISITED, false);
            ac.addBond(bond);
        }
        atom = ac.getAtom(0);
        sphere.add(atom);
        atom.setFlag(VISITED, true);
        breadthFirstSearch(ac, sphere, molecule);
        if (molecule.getAtomCount() == atomContainer.getAtomCount()) {
            flag = true;
        }
        return flag;
    }

    /**
     * Partitions the atoms in an GraphAtomContainer into covalently connected
     * components.
     *
     * @param atomContainer The GraphAtomContainer to be partitioned into
     * connected components, i.e. molecules
     * @return A MoleculeSet.
     *
     *
     */
    public static IAtomContainerSet splitMolecules(IAtomContainer atomContainer) {
        IAtomContainer ac = atomContainer.getBuilder().newInstance(IAtomContainer.class);
        IAtom atom;
        IElectronContainer eContainer;
        IAtomContainer molecule;
        IAtomContainerSet molecules = atomContainer.getBuilder().newInstance(IAtomContainerSet.class);
        List<IAtom> sphere = new ArrayList<>();
        for (int f = 0; f < atomContainer.getAtomCount(); f++) {
            atom = atomContainer.getAtom(f);
            atom.setFlag(VISITED, false);
            ac.addAtom(atom);
        }
        Iterator<IElectronContainer> eContainers = atomContainer.electronContainers().iterator();
        while (eContainers.hasNext()) {
            eContainer = eContainers.next();
            eContainer.setFlag(VISITED, false);
            ac.addElectronContainer(eContainer);
        }
        while (ac.getAtomCount() > 0) {
            atom = ac.getAtom(0);
            molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            sphere.clear();
            sphere.add(atom);
            atom.setFlag(VISITED, true);
            breadthFirstSearch(ac, sphere, molecule);
            molecules.addAtomContainer(molecule);
            ac.remove(molecule);
        }
        return molecules;
    }

    private EBIMolSplitter() {
    }
}
