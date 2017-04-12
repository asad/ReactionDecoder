/**
 *
 * Copyright (C) 2009-2015 Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this program; if not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Class that cleans redundant mappings from the solution set.
 * <OL>
 *
 * <lI>1: Stereo match, bond type, ring etc,
 * <lI>2: Fragment size,
 * <lI>3: Bond breaking energy
 *
 * </OL>
 *
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class PostFilter {

    /**
     *
     * Creates a new instance of Post Filter and removes redundant mapping(s).
     *
     * @param mappings
     * @return Filtered non-redundant mappings
     */
    public synchronized static List<Map<Integer, Integer>> filter(List<List<Integer>> mappings) {
        List<Map<Integer, Integer>> final_MAPPINGS = new ArrayList<>();
        if (mappings != null && !mappings.isEmpty()) {
            List<Map<Integer, Integer>> removeRedundantMapping = removeRedundantMapping(mappings);
            final_MAPPINGS.addAll(removeRedundantMapping);
        }
        return final_MAPPINGS;
    }

    private synchronized static boolean hasMap(Map<Integer, Integer> newMap, List<Map<Integer, Integer>> nonRedundantMapping) {
        for (Map<Integer, Integer> storedMap : nonRedundantMapping) {
            if (storedMap.equals(newMap)) {
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param mapping_org
     * @return
     */
    private synchronized static List<Map<Integer, Integer>> removeRedundantMapping(List<List<Integer>> mapping_org) {
        List<Map<Integer, Integer>> nonRedundantMapping = Collections.synchronizedList(new ArrayList<Map<Integer, Integer>>());
        for (List<Integer> M : mapping_org) {
            Map<Integer, Integer> newMap = getMappingMapFromList(M);
            if (!hasMap(newMap, nonRedundantMapping)) {
                nonRedundantMapping.add(newMap);
            }
        }
        return nonRedundantMapping;
    }

    private synchronized static Map<Integer, Integer> getMappingMapFromList(List<Integer> list) {
        Map<Integer, Integer> newMap = Collections.synchronizedSortedMap(new TreeMap<Integer, Integer>());
        for (int index = 0; index < list.size(); index += 2) {
            newMap.put(list.get(index), list.get(index + 1));
        }
        return newMap;
    }
}
