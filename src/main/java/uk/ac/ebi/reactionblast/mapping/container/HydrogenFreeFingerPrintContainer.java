/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.container;

//~--- non-JDK imports --------------------------------------------------------
import java.io.IOException;
import java.io.Serializable;
import java.util.BitSet;
import static java.util.Collections.synchronizedSortedMap;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.fingerprints.tools.Similarity.getTanimotoSimilarity;
import uk.ac.ebi.reactionblast.interfaces.IFingerPrintContainer;

//~--- classes ----------------------------------------------------------------
/**
 * @RCSfile: atomMapperTool.java,v
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
public class HydrogenFreeFingerPrintContainer implements IFingerPrintContainer, Serializable {

    //define the FINGER_SIZE of the fingerprint
    //NOTE: this should be a multiple of 64 and preferably not 1024 or 2048
    //as for these values we often get the random numbers for one-atom or
    //two-atom paths the same!
    private static final int FINGER_SIZE = 64 * 30;
    //depth search is set to 6, if not given explicitly as a parameter to Fingerprinter
    private static final int DEPTH_SEARCH = 8;
    private static final long serialVersionUID = 987987606869669691L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(HydrogenFreeFingerPrintContainer.class);

    /**
     *
     * @return
     */
    public static int getFingerPrintSize() {
        return FINGER_SIZE;
    }

    /**
     *
     * @return
     */
    public static int getFingerPrintDepth() {
        return DEPTH_SEARCH;
    }
    private final Map<String, BitSet> fingerPrintMap;

    //~--- constructors -------------------------------------------------------
    /**
     * HydrogenFreeFingerPrintContainer container
     */
    public HydrogenFreeFingerPrintContainer() {
        fingerPrintMap = synchronizedSortedMap(new TreeMap<>());
    }

    @Override
    public String toString() {
        return "HydrogenFreeFingerPrintContainer{" + "fingerPrintMap=" + fingerPrintMap + '}';
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        fingerPrintMap.clear();
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Erase(String Key) throws IOException {
        fingerPrintMap.remove(Key);
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
            fingerPrintMap.put(Key, Value);
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
        return fingerPrintMap.get(Key);
    }

    /**
     *
     * @param bitset
     * @return
     * @throws java.io.IOException
     */
    @Override
    public synchronized String getMoleculeID(BitSet bitset) throws IOException {
        String Key = null;
        for (Map.Entry<String, BitSet> map : fingerPrintMap.entrySet()) {
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
    synchronized public Map<String, BitSet> getFingerPrintMap()
            throws IOException {
        return synchronizedSortedMap(new TreeMap<>(fingerPrintMap));
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isKeyPresent(String Key) throws IOException {
        return fingerPrintMap.containsKey(Key);
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    @Override
    public synchronized void setValue(String Key, BitSet Value) throws IOException {
        fingerPrintMap.put(Key, Value);
    }

    /**
     *
     * @param value
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized boolean isValuePresent(BitSet value) throws IOException {
        for (BitSet bitset : fingerPrintMap.values()) {
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
     *
     * @return
     */
    public int getSize() {
        return fingerPrintMap.size();
    }

    @Override
    public boolean isEmpty() throws IOException {
        return fingerPrintMap.isEmpty();
    }

    @Override
    public void write() throws IOException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
