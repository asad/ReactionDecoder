/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping.cache;

import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 * @param <K>
 * @param <V>
 */
//Thread Safe Cache
public class ThreadSafeCache<K, V> implements Cache<K, V> {
    //Shared resource that needs protection

    private final Map<K, V> map;

    //Single instance kept
    private static final ThreadSafeCache SC = new ThreadSafeCache();

    //Access method
    public static ThreadSafeCache getInstance() {
        return SC;
    }

    //Private constructor to prevent instantiation
    private ThreadSafeCache() {
        //ConcurrentHashMap takes care of 
        //synchronization in a a multi-threaded 
        //environment
        map = new ConcurrentHashMap<>();
    }

    //Now this is thread safe
    @Override
    public void put(K key, V value) {
        map.put(key, value);
    }

    //Now this is thread safe
    @Override
    public V get(K key) {
        return map.get(key);
    }

    //Now this is thread safe
    /**
     * Check if Key Present
     *
     * @param key
     * @return
     */
    public boolean containsKey(K key) {
        return map.containsKey(key);
    }

    // CLEANUP method
    public void cleanup() {
        map.clear();
        Thread.yield();
    }

    /**
     * Size of the map
     *
     * @return
     */
    public int size() {
        return map.size();
    }

    /**
     * Size of the map
     *
     * @return
     */
    public Set<K> keySet() {
        return map.keySet();
    }
}
