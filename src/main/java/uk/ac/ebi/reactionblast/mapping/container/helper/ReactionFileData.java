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
package uk.ac.ebi.reactionblast.mapping.container.helper;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.synchronizedList;
import static java.util.Collections.unmodifiableCollection;
import java.util.List;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionFileData extends Object implements Serializable {

    private static final long serialVersionUID = 193790837047304701L;
    private List<String> _data = null;

    /**
     *
     */
    public ReactionFileData() {
        _data = synchronizedList(new ArrayList<>());
    }

    /**
     * @param data
     * @return true/false
     */
    public boolean addData(String data) {
        return _data.add(data.trim());
    }

    /**
     * @param index
     * @param data
     */
    public void addData(int index, String data) {
        _data.add(index, data.trim());
    }

    /**
     * @param data
     * @return true/false
     */
    public boolean isPresent(String data) {
        return _data.contains(data.trim());
    }

    /**
     * @param index
     * @return true/false
     */
    public String getValue(int index) {
        return _data.get(index);
    }

    /**
     * @return true/false
     */
    public Collection<String> getValues() {
        return unmodifiableCollection(_data);
    }

    /**
     * @return true/false
     */
    public boolean isEmpty() {
        return _data.isEmpty();
    }

    /**
     *
     * @return data size
     */
    public int getSize() {
        return _data.size();
    }

    /**
     * Print the data
     *
     * @return
     */
    @Override
    public String toString() {
        int index = 0;
        String data = "";
        for (String val : _data) {
            data += "Index: " + index++ + " Data: " + val;
        }
        return data;
    }
}
