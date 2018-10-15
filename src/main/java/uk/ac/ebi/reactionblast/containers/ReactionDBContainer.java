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

import java.io.Serializable;
import java.util.Collection;
import static java.util.Collections.synchronizedMap;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ReactionDBContainer implements Serializable {

    private static Map<String, ReactionInfoCollector> reactionsFingerprints = null;
    private static ReactionDBContainer ref = null;
    private static final long serialVersionUID = 19998987876L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(ReactionDBContainer.class);

    /**
     * Creates a new instance of CompoundContainer
     *
     * @return
     * @throws java.lang.Exception
     */
    public static synchronized ReactionDBContainer getInstance()
            throws Exception {
        if (ref == null) {

            // it's ok, we can call this constructor
            ref = new ReactionDBContainer();
        }

        return ref;
    }

    //~--- constructors -------------------------------------------------------
    /**
     *
     * @throws Exception
     */
    protected ReactionDBContainer() throws Exception {
        reactionsFingerprints = synchronizedMap(new HashMap<>());

    }

    /**
     *
     * @param Key Reaction ID
     * @return true or false
     */
    synchronized public boolean containsKey(String Key) {

        return reactionsFingerprints.containsKey(Key);
    }

    /**
     *
     * @return
     */
    synchronized public int size() {
        return reactionsFingerprints.size();
    }

    /**
     *
     * @return
     */
    synchronized public boolean isEmpty() {
        return reactionsFingerprints.isEmpty();
    }

    /**
     *
     * @param value
     * @return
     */
    synchronized public boolean containsValue(ReactionInfoCollector value) {
        return reactionsFingerprints.containsValue(value);
    }

    /**
     *
     * @param key
     * @return
     */
    synchronized public ReactionInfoCollector get(String key) {
        return reactionsFingerprints.get(key);
    }

    /**
     *
     * @param key
     * @param value
     * @return
     */
    synchronized public ReactionInfoCollector put(String key, ReactionInfoCollector value) {
        return reactionsFingerprints.put(key, value);
    }

    /**
     *
     * @param key
     * @return
     */
    synchronized public ReactionInfoCollector remove(String key) {
        return reactionsFingerprints.remove(key);
    }

    /**
     *
     * @param m
     */
    synchronized public void putAll(Map<String, ReactionInfoCollector> m) {
        reactionsFingerprints.putAll(m);
    }

    /**
     *
     */
    synchronized public void clear() {
        reactionsFingerprints.clear();
        ref = null;
    }

    /**
     *
     * @return
     */
    synchronized public Set<String> keySet() {
        return reactionsFingerprints.keySet();
    }

    /**
     *
     * @return
     */
    synchronized public Collection<ReactionInfoCollector> values() {
        return reactionsFingerprints.values();
    }

    /**
     *
     * @return
     */
    synchronized public Set<Entry<String, ReactionInfoCollector>> entrySet() {
        return reactionsFingerprints.entrySet();
    }
}
