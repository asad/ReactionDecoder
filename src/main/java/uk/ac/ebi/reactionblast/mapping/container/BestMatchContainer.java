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

//~--- JDK imports ------------------------------------------------------------
import java.io.IOException;
import java.io.Serializable;
import static java.util.Collections.synchronizedMap;
import java.util.HashMap;
import java.util.Map;
import static java.util.logging.Level.SEVERE;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import uk.ac.ebi.reactionblast.mapping.container.helper.Key;
import uk.ac.ebi.reactionblast.mapping.interfaces.BestMatch;
import uk.ac.ebi.reactionblast.mapping.interfaces.IKey;

//~--- classes --------------------------------------
/**
 *
 *
 * @Author: Syed Asad Rahman: asad@ebi.ac.uk
 * @Copyright (C) 2004-2004 The Atom Mapper Tool (AMT) project
 * @Date: 2004/06/3
 * @Revision: 1.10
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
 */
public class BestMatchContainer extends BestMatch implements Serializable {

    private final static ILoggingTool LOGGER
            = createLoggingTool(BestMatchContainer.class);
    private static final long serialVersionUID = 10947239472059259L;
    private final Map<IKey, AtomAtomMapping> mcsAtomMap;
    private final Map<IKey, Integer> fragmentCount;
    private final Map<IKey, Double> bondBreakingEnergy;
    private final Map<IKey, Double> stereoScore;
    private final Map<IKey, Double> similarity;

    //~--- constructors -------------------------------------------------------
    /**
     *
     */
    public BestMatchContainer() {
        mcsAtomMap = synchronizedMap(new HashMap<>());
        fragmentCount = synchronizedMap(new HashMap<>());
        bondBreakingEnergy = synchronizedMap(new HashMap<>());
        stereoScore = synchronizedMap(new HashMap<>());
        similarity = synchronizedMap(new HashMap<>());
        // System.out.println("FingerPrint Map Created");
    }

    @Override
    public String toString() {
        return "BestMatchContainer{" + "mcsAtomMap=" + mcsAtomMap + ", fragmentCount=" + fragmentCount + ", bondBreakingEnergy=" + bondBreakingEnergy + ", stereoScore=" + stereoScore + ", similarity=" + similarity + '}';
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        mcsAtomMap.clear();
        fragmentCount.clear();
        bondBreakingEnergy.clear();
        stereoScore.clear();
        similarity.clear();
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Erase(int indexI, int indexJ) throws IOException {
        IKey key = new Key(indexI, indexJ);
        mcsAtomMap.remove(key);
        fragmentCount.remove(key);
        bondBreakingEnergy.remove(key);
        stereoScore.remove(key);
        similarity.remove(key);
    }

    //~--- get methods --------------------------------------------------------
    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     * @throws IOException
     */
    @Override
    public synchronized AtomAtomMapping getAtomMatch(int indexI, int indexJ)
            throws IOException {
        IKey key = new Key(indexI, indexJ);
        if (mcsAtomMap.containsKey(key)) {
            return mcsAtomMap.get(key);
        } else {
            try {
                throw new CDKException("Key not found:" + key + " in " + mcsAtomMap.keySet());
            } catch (CDKException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }
        return null;
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    @Override
    public synchronized double getBondEnergy(int indexI, int indexJ) {
        IKey key = new Key(indexI, indexJ);
        return bondBreakingEnergy.containsKey(key) ? bondBreakingEnergy.get(key) : 0.;
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    @Override
    public synchronized int getTotalFragmentCount(int indexI, int indexJ) {
        IKey key = new Key(indexI, indexJ);
        return fragmentCount.containsKey(key) ? fragmentCount.get(key) : 0;
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param value
     */
    @Override
    public synchronized void setBondEnergy(int indexI, int indexJ, double value) {
        IKey key = new Key(indexI, indexJ);
        bondBreakingEnergy.put(key, value);
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param value
     */
    @Override
    public synchronized void setTotalFragmentCount(int indexI, int indexJ, Integer value) {
        IKey key = new Key(indexI, indexJ);
        fragmentCount.put(key, value);
    }

    //~--- set methods --------------------------------------------------------
    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param matchingAtoms
     * @throws IOException
     */
    @Override
    public synchronized void putBestMapping(int indexI, int indexJ, AtomAtomMapping matchingAtoms)
            throws IOException {
        IKey key = new Key(indexI, indexJ);
        mcsAtomMap.put(key, matchingAtoms);
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     * @throws IOException
     */
    @Override
    public synchronized boolean containsKey(int indexI, int indexJ)
            throws IOException {
        IKey key = new Key(indexI, indexJ);
        return mcsAtomMap.containsKey(key);
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param stereoVal
     */
    @Override
    public synchronized void setStereoScore(int indexI, int indexJ, double stereoVal) {
        IKey key = new Key(indexI, indexJ);
        stereoScore.put(key, stereoVal);
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    @Override
    public synchronized double getStereoScore(int indexI, int indexJ) {
        IKey key = new Key(indexI, indexJ);
        return stereoScore.containsKey(key) == true ? stereoScore.get(key) : 0.0;
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @param _simVal
     */
    @Override
    public synchronized void setGraphSimilarity(int indexI, int indexJ, double _simVal) {
        IKey key = new Key(indexI, indexJ);
        similarity.put(key, _simVal);
    }

    /**
     * String IKey = substrateIndex + "_" + productIndex;
     *
     * @param indexI
     * @param indexJ
     * @return
     */
    @Override
    public synchronized double getGraphSimilarity(int indexI, int indexJ) {
        IKey key = new Key(indexI, indexJ);
        return similarity.containsKey(key) == true ? similarity.get(key) : 0.0;
    }
}
