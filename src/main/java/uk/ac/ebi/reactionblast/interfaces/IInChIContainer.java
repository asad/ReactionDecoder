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

package uk.ac.ebi.reactionblast.interfaces;

import java.io.IOException;
import java.util.Map;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public interface IInChIContainer {

    /**
     *
     * @throws java.io.IOException
     */
    void Clear() throws IOException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    void Erase(String Key) throws IOException;

    /**
     *
     * @return
     * @throws CloneNotSupportedException
     */
    Object clone() throws CloneNotSupportedException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    String getInChI(String Key) throws IOException;

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    Map<String, String> getInChIMap() throws IOException;

    /**
     *
     * @param Value
     * @return
     * @throws java.io.IOException
     */
    String getMoleculeID(String Value) throws IOException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    boolean isKeyPresent(String Key) throws IOException;

    /**
     *
     * @param Value
     * @throws java.io.IOException
     * @return
     */
    boolean isValuePresent(String Value) throws IOException;

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    void put(String Key, String Value) throws IOException;

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    void setValue(String Key, String Value) throws IOException;

    /**
     *
     * @throws java.io.IOException
     */
    void write() throws IOException;

}
