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
import java.util.ArrayList;
import static java.util.Collections.synchronizedMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.exception.CDKException;
import uk.ac.ebi.reactionblast.mapping.container.helper.MolMapping;

/**
 * @RCSfile: atomMapperTool.java,v
 *
 * @Author: Syed Asad Rahman
 * @Date: 2009/06/3
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
public class MoleculeMoleculeMapping implements Serializable {

    private static final long serialVersionUID = 1094750239472059259L;

    //~--- fields -------------------------------------------------------------
    private final Map<String, List<MolMapping>> reactant_product_mapping_map;

    //~--- constructors -------------------------------------------------------
    /**
     *
     */
    public MoleculeMoleculeMapping() {
        reactant_product_mapping_map = synchronizedMap(new HashMap<>());
    }

    @Override
    public String toString() {
        return "MoleculeMoleculeMapping{" + "reactant_product_mapping_map=" + reactant_product_mapping_map + '}';
    }

    /**
     *
     * @throws java.io.IOException
     */
    public synchronized void Clear() throws IOException {
        reactant_product_mapping_map.clear();
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    public synchronized void Erase(String Key) throws IOException {
        reactant_product_mapping_map.remove(Key);
    }

    /**
     *
     * @param Key
     * @return
     * @throws java.io.IOException
     */
    public synchronized boolean isPresent(String Key)
            throws IOException {
        return reactant_product_mapping_map.containsKey(Key);
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    public synchronized void setMolMappings(String Key, List<MolMapping> Value) throws
            IOException {
        reactant_product_mapping_map.put(Key, Value);
        // Stores Reaction ID and RPAIR ID as Value in ArrayList
    }

    /**
     *
     * @param RID
     * @return
     * @throws CDKException
     */
    public synchronized List<MolMapping> getMolMappings(String RID) throws CDKException {
        return reactant_product_mapping_map.containsKey(RID) == true ? reactant_product_mapping_map.get(RID) : null;
    }

    /**
     *
     * @return Reaction count with RPAIR
     */
    public synchronized long getCount() {
        return reactant_product_mapping_map.size();
    }

    /**
     *
     * @return
     */
    public synchronized Set<String> getKeySet() {
        return reactant_product_mapping_map.keySet();
    }

    /**
     *
     * @return
     */
    public synchronized Set<Map.Entry<String, List<MolMapping>>> getEntrySet() {
        return reactant_product_mapping_map.entrySet();
    }

    /**
     *
     * @param reactionID
     * @param rName
     * @param pName
     * @return
     */
    public synchronized List<MolMapping> getMapping(String reactionID, String rName, String pName) {
        List<MolMapping> mMap = reactant_product_mapping_map.get(reactionID);
        List<MolMapping> mappedMap = new ArrayList<>();
        for (MolMapping map : mMap) {
            if ((map.getTarget().equalsIgnoreCase(rName) && map.getQuery().equalsIgnoreCase(pName))
                    || (map.getTarget().equalsIgnoreCase(pName) && map.getQuery().equalsIgnoreCase(rName))) {
                mappedMap.add(map);
            }
        }
        return mappedMap;
    }
}
