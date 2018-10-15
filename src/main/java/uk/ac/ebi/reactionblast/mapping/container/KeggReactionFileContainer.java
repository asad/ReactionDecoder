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

 /*
 * ReactionContainer.java
 *
 * Created on 19 January 2006, 10:48
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast.mapping.container;

//~--- JDK imports ------------------------------------------------------------
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import uk.ac.ebi.reactionblast.mapping.container.helper.ReactionFileData;

//~--- classes ----------------------------------------------------------------
/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class KeggReactionFileContainer implements Serializable {

    private static final long serialVersionUID = 1094750239472059259L;
    private static TreeMap<String, Map<String, ReactionFileData>> reactionMap = null;
    private static KeggReactionFileContainer ref = null;
    //~--- fields -------------------------------------------------------------
    private final static String entry = "ENTRY";
    private final static String name = "NAME";
    private final static String equat = "EQUATION";
    private final static String defin = "DEFINITION";
    private final static String path = "PATHWAY";
    private final static String rpair = "RPAIR";
    private final static String ortho = "ORTHOLOGY";
    private final static String enzyme = "ENZYME";
    private final static String end = "///";
    private final static String comment = "COMMENT";
    private final static String remark = "REMARK";

    /**
     * Creates a new instance of ReactionContainer
     *
     * @return
     * @throws java.lang.Exception
     */
    public static synchronized KeggReactionFileContainer getInstance() throws Exception {
        if (ref == null) {

            // it's ok, we can call this constructor
            ref = new KeggReactionFileContainer();
        }

        return ref;
    }

    //~--- constructors -------------------------------------------------------
    /**
     *
     * @throws Exception
     */
    protected KeggReactionFileContainer() throws Exception {
        reactionMap = new TreeMap<>();
    }

    //~--- methods ------------------------------------------------------------
    /**
     * Clear the container
     */
    synchronized public void Clear() {
        reactionMap.clear();
        ref = null;
    }

    /**
     *
     * @param _recID reaction ID
     * @param Key {NAME/EQUATION etc}
     * @return True/False
     */
    public synchronized boolean containsKey(String _recID, String Key) {
        Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);
        return DataMap.containsKey(Key) == true;
    }

    //~--- get methods --------------------------------------------------------
    /**
     *
     * @return Reaction Ids in the container
     */
    public synchronized Set<String> getKeySet() {
        return reactionMap.keySet();
    }

    /**
     *
     * @param _recID reaction ID
     * @return Defination of the reaction
     */
    public synchronized String getDefination(String _recID) {
        String _data = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(defin)) {
                _data = DataMap.get(defin).getValue(0);
            }
        }

        return _data.trim();
    }

    /**
     *
     * @param _recID
     * @return entry ID of the reaction
     */
    public synchronized String getEntry(String _recID) {
        String _data = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(entry)) {
                _data = DataMap.get(entry).getValue(0);
            }
        }

        return _data.trim();
    }

    /**
     *
     * @param _recID
     * @return enzymes of this reaction
     */
    public synchronized List<String> getEnzyme(String _recID) {
        List<String> ecData = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(enzyme)) {
                ecData = new ArrayList<>();
                ReactionFileData eData = DataMap.get(enzyme);
                for (String data : eData.getValues()) {
//                    System.out.println("RID " + _recID + " EC: " + data);
                    ecData.add(data);
                }
            }
        }

        return ecData;
    }

    /**
     *
     * @param _recID
     * @return reaction equation
     */
    public synchronized String getEquation(String _recID) {
        String _data = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(equat)) {
                _data = DataMap.get(equat).getValue(0);
            }
        }

        return _data.trim();
    }

    /**
     *
     * @return
     */
    synchronized Set<String> getKeys() {
        return reactionMap.keySet();
    }

//    /**
//     *
//     * @return
//     */
//    synchronized public TreeMap<String, Map<String, ReactionFileData>> getMap() {
//        return reactionMap;
//    }
    /**
     *
     * @param _recID reaction ID
     * @return Name of the reaction
     */
    synchronized public String getName(String _recID) {
        String _data = null;
        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(name)) {
                _data = DataMap.get(name).getValue(0).trim();
            }
        }
        return _data;
    }

    /**
     *
     * @param _recID
     * @return Pathway Names
     */
    synchronized public List<String> getPathway(String _recID) {
        List<String> _data = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(path)) {
                _data = new ArrayList<>();
                for (String data : DataMap.get(path).getValues()) {
                    _data.add(data);
                }
            }
        }
        return _data;
    }

    /**
     *
     * @param _recID reaction ID
     * @return Reaction crosslink(s)
     */
    synchronized public List<String> getRemark(String _recID) {
        List<String> _data = null;

        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);

            if (DataMap.containsKey(remark)) {
                _data = new ArrayList<>();
                for (String data : DataMap.get(remark).getValues()) {
                    _data.add(data);
                }
            }
        }

        return _data;
    }

    /**
     *
     * @param _recID reaction ID
     * @return RPAIR(s)
     */
    public synchronized String getRpair(String _recID) {
        String _data = null;
        if (reactionMap.containsKey(_recID)) {
            Map<String, ReactionFileData> DataMap = reactionMap.get(_recID);
            if (DataMap.containsKey(rpair)) {
                _data = DataMap.get(rpair).getValue(0);
            }
        }

        return _data.trim();
    }

    /**
     *
     * @return reaction count in the container
     */
    public synchronized Integer getCount() {
        return reactionMap.size();
    }

    /**
     *
     * @param key
     * @param DataMap
     */
    public synchronized void put(String key, Map<String, ReactionFileData> DataMap) {
        reactionMap.put(key, DataMap);
    }
}
