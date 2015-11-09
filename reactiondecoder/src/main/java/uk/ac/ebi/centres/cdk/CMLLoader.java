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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.ChemFile;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.config.Isotopes;
import uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator;

/**
 * @author John May
 */
public class CMLLoader {

    public static IAtomContainer loadCML(InputStream in) {
        CMLReader reader = new CMLReader(in);
        try {
            IChemFile chemfile = reader.read(new ChemFile());
            Iterator<IAtomContainer> iterator = ChemFileManipulator.getAllAtomContainers(chemfile).iterator();
            if (iterator.hasNext()) {
                IAtomContainer container = iterator.next();
                // due to a bug need to reconfigure molecule
                for (IAtom atom : container.atoms()) {
                    atom.setAtomicNumber(PeriodicTable.getAtomicNumber(atom.getSymbol()));
                    if (!atom.getSymbol().equals("R")) {
                        atom.setMassNumber(Isotopes.getInstance().getMajorIsotope(atom.getSymbol()).getMassNumber());
                    }
                }
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(container);
                return container;
            }
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException ex) {
            Logger.getLogger(CMLLoader.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            if (reader != null) {
                reader.close();
            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return new AtomContainer();
    }
}
