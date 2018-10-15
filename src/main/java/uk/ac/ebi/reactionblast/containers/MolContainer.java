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
package uk.ac.ebi.reactionblast.containers;

//~--- non-JDK imports --------------------------------------------------------
import java.io.IOException;
import static java.util.Collections.synchronizedSortedMap;
import static java.util.Collections.unmodifiableMap;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.Substructure;
import uk.ac.ebi.reactionblast.interfaces.IMolContainer;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

//~--- classes ----------------------------------------------------------------
/**
 * @RCSfile: atomMapperTool.java,v
 *
 * @Author: Syed Asad Rahman
 * @Date: 2004/06/3
 * @Revision: 1.10
 *
 * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
 *
 * @Contact: asad@ebi.ac.uk
 *
 * @This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version. All we ask is that proper credit is given for our
 * work, which includes - but is not limited to - adding the above copyright
 * notice to the beginning of your source code files, and to any copyright
 * notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *
 */
public class MolContainer implements IMolContainer {

    /*
     * Singleton Pattern Implementation
     */
    private static MolContainer _instance = null;
    private static Map<String, IAtomContainer> molContainer = null;
    private final static ILoggingTool LOGGER
            = createLoggingTool(MolContainer.class);

    /**
     *
     * @return
     */
    public static synchronized MolContainer getInstance() {
        if (_instance == null) {
            _instance = new MolContainer();
        }

        return _instance;
    }

    //~--- constructors -------------------------------------------------------
    private MolContainer() {
        molContainer = synchronizedSortedMap(new TreeMap<>());
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        molContainer.clear();
        molContainer = synchronizedSortedMap(new TreeMap<String, IAtomContainer>());
    }

    /**
     *
     * @param key
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Erase(String key) throws IOException {
        molContainer.remove(key);
    }

    /**
     *
     * @param key
     * @throws java.io.IOException
     */
    @Override
    public synchronized void put(String key, IAtomContainer Value) throws IOException {
        try {
            molContainer.put(key, Value);
        } catch (Exception e) {
            LOGGER.debug(e);
        }
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param key
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized IAtomContainer getAtomContainer(String key)
            throws IOException {
        return molContainer.get(key);
    }

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized Map<String, IAtomContainer> getAtomContainerMap() throws IOException {
        return unmodifiableMap(molContainer);
    }

    /**
     *
     * @param key
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isKeyPresent(String key) throws IOException {
        return molContainer.containsKey(key);
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param key
     * @param value
     * @throws java.io.IOException
     */
    @Override
    public synchronized void add(String key, IAtomContainer value)
            throws IOException {
        molContainer.put(key, value);
    }

    /**
     *
     * @param _queryMol mol to be compared
     * @param _targetMol clean mol from the container
     * @param removeHydrogen
     * @return
     * @throws Exception
     */
    public synchronized boolean isIdentical(IAtomContainer _queryMol,
            IAtomContainer _targetMol,
            boolean removeHydrogen) throws Exception {

        _targetMol = cloneWithIDs(_targetMol);
        if (_queryMol.getAtomCount() == 1 && _targetMol.getAtomCount() == 1) {
            IAtom a = _queryMol.atoms().iterator().next();
            IAtom b = _targetMol.atoms().iterator().next();
            return a.getSymbol().equalsIgnoreCase(b.getSymbol())
                    && Objects.equals(a.getFormalCharge(), b.getFormalCharge());
        }
        return isSubgraphIdentical(_queryMol, _targetMol, removeHydrogen);
    }

    private synchronized boolean isSubgraphIdentical(IAtomContainer _mol,
            IAtomContainer _rMol,
            boolean removeHydrogen) throws CDKException, IOException {
//        System.out.println("Graph matching");

        IAtomContainer mol1 = _mol;
        IAtomContainer mol2 = _rMol;

        Substructure mcs = new Substructure(mol1, mol2, true, true, true, false);
        mcs.setChemFilters(false, false, false);
        return mcs.isSubgraph() && !mcs.isStereoMisMatch()
                && mol1.getAtomCount() == mol2.getAtomCount();
    }

    /**
     *
     * @param key
     * @param mol
     * @return
     * @throws Exception
     */
    @Override
    public synchronized boolean compareAtomContainer(String key, IAtomContainer mol) throws Exception {
        mol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
        try {
            boolean flag = molContainer.containsKey(key);
            if (flag && mol.getAtomCount() > 0) {
                IAtomContainer molFromContainer = molContainer.get(key);
                return isIdentical(mol, molFromContainer, true);
            }
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return false;
    }

    @Override
    public synchronized String getMoleculeID(IAtomContainer mol) throws Exception {
        IAtomContainer queryMol = removeHydrogensExceptSingleAndPreserveAtomID(mol);
        percieveAtomTypesAndConfigureAtoms(queryMol);
        CDKHydrogenAdder instance = CDKHydrogenAdder.getInstance(queryMol.getBuilder());
        for (IAtom atom : queryMol.atoms()) {
            try {
                instance.addImplicitHydrogens(queryMol, atom);
            } catch (CDKException e) {
                LOGGER.error("WARNING: Error in adding H to the molecule");
            }
        }

        aromatizeMolecule(queryMol);

        for (Map.Entry<String, IAtomContainer> map : molContainer.entrySet()) {
            String key = map.getKey();
            IAtomContainer tMol = map.getValue();
            if (isIdentical(queryMol, tMol, true)) {
                return key;
            }
        }
        //System.LOGGER.debug("Error: Unable to Find AtomContainer ID!!!");
        return null;
    }

    @Override
    public synchronized boolean isValuePresent(IAtomContainer Value) throws IOException {
        return molContainer.containsValue(Value);
    }

    /**
     *
     * @return
     */
    public synchronized boolean isEmpty() {
        return molContainer.isEmpty();
    }

    @Override
    public void write() throws IOException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
