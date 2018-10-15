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
import java.util.BitSet;
import static java.util.Collections.synchronizedSortedMap;
import static java.util.Collections.unmodifiableMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.interfaces.IFingerPrintContainer;

//~--- classes ----------------------------------------------------------------
/**
 * @RCSfile: FingerPrintContainer.java,v
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
public class FingerPrintContainer implements IFingerPrintContainer {

    /*
     * Singleton Pattern Implementation
     */
    private static FingerPrintContainer _instance = null;
    private static Map<String, BitSet> FingerPrintMap = null;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(FingerPrintContainer.class);

    /**
     *
     * @return
     */
    public static synchronized FingerPrintContainer getInstance() {
        if (_instance == null) {
            _instance = new FingerPrintContainer();
        }

        return _instance;
    }

    //~--- constructors -------------------------------------------------------
    private FingerPrintContainer() {
        FingerPrintMap = synchronizedSortedMap(new TreeMap<>());
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        FingerPrintMap.clear();
        FingerPrintMap = synchronizedSortedMap(new TreeMap<>());
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Erase(String Key) throws IOException {
        FingerPrintMap.remove(Key);
    }

    /**
     *
     * @return
     */
    public synchronized Integer getCount() {
        return FingerPrintMap.size();
    }

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    @Override
    public synchronized void put(String Key, BitSet Value) throws IOException {
        try {
            FingerPrintMap.put(Key, Value);
        } catch (Exception e) {
            LOGGER.error(SEVERE, null, e);
        }
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized BitSet getFingerPrint(String Key) throws IOException {
        BitSet value = FingerPrintMap.get(Key);
        return value;
    }

    /**
     *
     * @param bitset
     * @return
     * @throws java.io.IOException
     */
    @Override
    public synchronized String getMoleculeID(BitSet bitset)
            throws IOException {
        String Key = null;
        for (Map.Entry<String, BitSet> map : FingerPrintMap.entrySet()) {
            String key = map.getKey();
            try {
                if (getTanimotoSimilarity(map.getValue(), bitset) == 1.0) {
                    Key = key;
                    break;
                }
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        //System.LOGGER.debug("Error: Unable to Find AtomContainer ID!!!");
        return Key;
    }

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized Map<String, BitSet> getFingerPrintMap() throws IOException {
        return unmodifiableMap(FingerPrintMap);
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isKeyPresent(String Key) throws IOException {
        return FingerPrintMap.containsKey(Key);
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    @Override
    public synchronized void setValue(String Key, BitSet Value)
            throws IOException {
//        System.out.println("KEY " + Key + " val: " + Value.cardinality());
        FingerPrintMap.put(Key, Value);
//        System.out.println("FingerPrintMap " + FingerPrintMap.size() + " val: " + Value.cardinality());
    }

    /**
     *
     * @return
     */
    public synchronized Set<String> getCompoundIDSet() {
        return FingerPrintMap.keySet();
    }

    /**
     *
     * @param value
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isValuePresent(BitSet value) throws IOException {
        for (BitSet bitset : FingerPrintMap.values()) {
            try {
                if (getTanimotoSimilarity(value, bitset) == 1.0) {
                    return true;
                }
            } catch (Exception ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        return false;
    }

    /**
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isEmpty() throws IOException {
        return FingerPrintMap.isEmpty();
    }

    @Override
    public void write() throws IOException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
