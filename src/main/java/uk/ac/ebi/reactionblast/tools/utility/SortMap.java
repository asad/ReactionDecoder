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
package uk.ac.ebi.reactionblast.tools.utility;

import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class SortMap {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(SortMap.class);

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Object, Double> valueInAscendingOrder(Map<Object, Double> map) {
        List<Map.Entry<Object, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Object, Double> entry, Map.Entry<Object, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Aceending Order Sort)
        );
        // LOGGER.info(list);
        Map<Object, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Object, Double> valueInDescendingOrder(Map<Object, Double> map) {
        List<Map.Entry<Object, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Object, Double> entry, Map.Entry<Object, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() < entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Decending Order Sort)
        );
        // LOGGER.info(list);
        Map<Object, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Integer, Integer> intValueInAscendingOrder(Map<Integer, Integer> map) {
        List<Map.Entry<Integer, Integer>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Integer, Integer> entry, Map.Entry<Integer, Integer> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Aceending Order Sort)
        );
        // LOGGER.info(list);
        Map<Integer, Integer> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Integer, Integer> intValueInDescendingOrder(Map<Integer, Integer> map) {
        List<Map.Entry<Integer, Integer>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Integer, Integer> entry, Map.Entry<Integer, Integer> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() < entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Decending Order Sort)
        );
        // LOGGER.info(list);
        Map<Integer, Integer> result = new LinkedHashMap<>();
        for (Map.Entry<Integer, Integer> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Double, Double> doubleValueInAscendingOrder(Map<Double, Double> map) {
        List<Map.Entry<Double, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Double, Double> entry, Map.Entry<Double, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Aceending Order Sort)
        );
        // LOGGER.info(list);
        Map<Double, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    /**
     *
     * @param map
     * @return
     */
    public static synchronized Map<Double, Double> doubleValueInDescendingOrder(Map<Double, Double> map) {
        List<Map.Entry<Double, Double>> list = new LinkedList<>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        sort(list, (Map.Entry<Double, Double> entry, Map.Entry<Double, Double> entry1) -> (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() < entry1.getValue() ? 1 : -1)) // Return 0 for bond match, -1 for less than and +1 for more then (Decending Order Sort)
        );
        // LOGGER.info(list);
        Map<Double, Double> result = new LinkedHashMap<>();
        list.forEach((entry) -> {
            result.put(entry.getKey(), entry.getValue());
        });
        return result;
    }

    private SortMap() {
    }
}
