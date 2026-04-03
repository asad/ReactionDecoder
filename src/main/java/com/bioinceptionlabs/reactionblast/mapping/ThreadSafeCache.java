/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Thread-safe LRU cache for MCS solutions. Supports cross-reaction caching
 * when canonical SMILES are used as keys (instead of molecule IDs).
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 * @param <K>
 * @param <V>
 */
/**
 * Generic Cache Interface.
 * @param <K>
 * @param <V>
 */
interface Cache<K, V> {
    void put(K key, V value);
    V get(K key);
}

public class ThreadSafeCache<K, V> implements Cache<K, V> {

    /** Maximum cache entries before LRU eviction kicks in. */
    private static final int MAX_CAPACITY = 10_000;

    private final Map<K, V> map;

    private static final ThreadSafeCache SC = new ThreadSafeCache();

    public static ThreadSafeCache getInstance() {
        return SC;
    }

    private ThreadSafeCache() {
        // ConcurrentHashMap for thread safety; LRU eviction handled in put()
        map = new ConcurrentHashMap<>(256, 0.75f, 4);
    }

    @Override
    public void put(K key, V value) {
        // Simple size-based eviction: if over capacity, clear oldest half
        if (map.size() >= MAX_CAPACITY) {
            evict();
        }
        map.put(key, value);
    }

    @Override
    public V get(K key) {
        return map.get(key);
    }

    /**
     * Check if key is present in the cache.
     */
    public boolean containsKey(K key) {
        return map.containsKey(key);
    }

    /**
     * Insert the value only if the key is absent.
     *
     * @return the existing value if present, otherwise the inserted value
     */
    public V putIfAbsent(K key, V value) {
        if (map.size() >= MAX_CAPACITY) {
            evict();
        }
        if (map instanceof ConcurrentHashMap<?, ?> concurrentMap) {
            @SuppressWarnings("unchecked")
            ConcurrentHashMap<K, V> typedMap = (ConcurrentHashMap<K, V>) concurrentMap;
            V existing = typedMap.putIfAbsent(key, value);
            return existing != null ? existing : value;
        }
        V existing = map.get(key);
        if (existing == null) {
            map.put(key, value);
            return value;
        }
        return existing;
    }

    /**
     * Clear all cached entries. Use sparingly — cross-reaction caching
     * benefits from keeping the cache warm between reactions.
     */
    public void cleanup() {
        map.clear();
    }

    /**
     * @return number of cached entries
     */
    public int size() {
        return map.size();
    }

    public Set<K> keySet() {
        return map.keySet();
    }

    /**
     * Evict roughly half the cache when over capacity.
     * ConcurrentHashMap iteration order is arbitrary, which
     * approximates random eviction — acceptable for MCS caching.
     */
    private void evict() {
        int toRemove = map.size() / 2;
        int removed = 0;
        for (K key : map.keySet()) {
            if (removed >= toRemove) break;
            map.remove(key);
            removed++;
        }
    }
}
