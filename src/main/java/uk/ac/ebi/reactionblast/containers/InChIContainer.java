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
import java.util.TreeMap;

import uk.ac.ebi.reactionblast.interfaces.IInChIContainer;

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
public class InChIContainer implements IInChIContainer, Cloneable {

    /* Singleton Pattern Implementation */
    private static InChIContainer _instance = null;
    private static Map<String, String> InChIMap = null;

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    public static int getCount() throws IOException {
        return InChIMap.size();
    }

    /**
     *
     * @return
     */
    public static synchronized InChIContainer getInstance() {
        if (_instance == null) {
            _instance = new InChIContainer();
        }
        return _instance;
    }

    //~--- constructors -------------------------------------------------------
    private InChIContainer() {
        InChIMap = synchronizedSortedMap(new TreeMap<>());
    }

    //~--- methods ------------------------------------------------------------
    /**
     *
     * @throws java.io.IOException
     */
    @Override
    public synchronized void Clear() throws IOException {
        InChIMap.clear();
        _instance = null;
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    @Override
    public void Erase(String Key) throws IOException {
        InChIMap.remove(Key);
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();

        // that'll teach 'em
    }

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    @Override
    synchronized public void put(String Key, String Value) throws IOException {
        try {
            if (Value != null) {
                InChIMap.put(Key, Value);
            }
        } catch (Exception e) {
            e.printStackTrace();
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
    synchronized public String getInChI(String Key)
            throws IOException {
        String value = InChIMap.get(Key);
        return value == null ? "" : value;
    }

    /**
     *
     * @param Value
     * @return
     * @throws java.io.IOException
     */
    @Override
    public synchronized String getMoleculeID(String Value) throws IOException {
        String Key = "Key Not Found";
        for (Map.Entry<String, String> map : InChIMap.entrySet()) {
            if (map.getValue().equals(Value)) {
                return Key;
            }
        }//System.LOGGER.debug("Error: Unable to Find AtomContainer ID!!!");
        return Key;
    }

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    @Override
    public synchronized Map<String, String> getInChIMap() throws IOException {
        return unmodifiableMap(InChIMap);
    }

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    @Override
    synchronized public boolean isKeyPresent(String Key) throws IOException {
        boolean flag = InChIMap.containsKey(Key);

        return flag;
    }

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    @Override
    synchronized public void setValue(String Key, String Value)
            throws IOException {
        InChIMap.put(Key, Value);
    }

    /**
     *
     * @param Value
     * @throws java.io.IOException
     * @return
     */
    @Override
    synchronized public boolean isValuePresent(String Value) throws IOException {
        boolean flag = InChIMap.containsValue(Value);
        return flag;
    }

    @Override
    public void write() throws IOException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
//~ Formatted by Jindent --- http://www.jindent.com

