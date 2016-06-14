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

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.config.Isotopes.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.ChemFile;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties;
import static org.openscience.cdk.tools.manipulator.ChemFileManipulator.getAllAtomContainers;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getAtomicNumber;

/**
 * @author John May
 */
public class CMLLoader {

    private static final Logger LOG = getLogger(CMLLoader.class.getName());

    /**
     *
     * @param in
     * @return
     */
    public static IAtomContainer loadCML(InputStream in) {
        CMLReader reader = new CMLReader(in);
        try {
            IChemFile chemfile = reader.read(new ChemFile());
            Iterator<IAtomContainer> iterator = getAllAtomContainers(chemfile).iterator();
            if (iterator.hasNext()) {
                IAtomContainer container = iterator.next();
                // due to a bug need to reconfigure molecule
                for (IAtom atom : container.atoms()) {
                    atom.setAtomicNumber(getAtomicNumber(atom.getSymbol()));
                    if (!atom.getSymbol().equals("R")) {
                        atom.setMassNumber(getInstance().getMajorIsotope(atom.getSymbol()).getMassNumber());
                    }
                }
                percieveAtomTypesAndConfigureUnsetProperties(container);
                return container;
            }
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException ex) {
            getLogger(CMLLoader.class.getName()).log(SEVERE, null, ex);
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

    private CMLLoader() {
    }
}
