/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import java.lang.ref.SoftReference;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Thread-safe cache for MCS solutions backed by {@link SoftReference} values.
 * <p>
 * Under normal heap pressure the cache behaves like a regular map — entries
 * remain reachable and provide O(1) MCS reuse across reactions with identical
 * molecule pairs.  When the JVM is low on memory the GC is free to reclaim
 * any soft-referenced value; a subsequent {@link #get} for that key simply
 * returns {@code null} and the caller falls through to a fresh MCS computation.
 * <p>
 * A hard capacity limit ({@link #MAX_CAPACITY}) prevents unbounded growth of
 * the key set itself; when reached, approximately half the entries are evicted.
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 * @param <K> key type (typically a canonical SMILES pair key)
 * @param <V> value type (typically {@code MCSSolution})
 */
interface Cache<K, V> {
    void put(K key, V value);
    V get(K key);
}

public class ThreadSafeCache<K, V> implements Cache<K, V> {

    /** Maximum number of key entries before random eviction kicks in. */
    private static final int MAX_CAPACITY = 500;

    private final ConcurrentHashMap<K, SoftReference<V>> map;

    @SuppressWarnings("rawtypes")
    private static final ThreadSafeCache SC = new ThreadSafeCache();

    @SuppressWarnings("unchecked")
    public static <K, V> ThreadSafeCache<K, V> getInstance() {
        return SC;
    }

    private ThreadSafeCache() {
        map = new ConcurrentHashMap<>(256, 0.75f, 4);
    }

    @Override
    public void put(K key, V value) {
        if (map.size() >= MAX_CAPACITY) {
            evict();
        }
        map.put(key, new SoftReference<>(value));
    }

    @Override
    public V get(K key) {
        SoftReference<V> ref = map.get(key);
        if (ref == null) {
            return null;
        }
        V value = ref.get();
        if (value == null) {
            // Referent was GC'd — remove the stale key
            map.remove(key);
        }
        return value;
    }

    /**
     * Check if key is present and its referent is still alive.
     */
    public boolean containsKey(K key) {
        SoftReference<V> ref = map.get(key);
        if (ref == null) {
            return false;
        }
        if (ref.get() == null) {
            map.remove(key);
            return false;
        }
        return true;
    }

    /**
     * Insert the value only if the key is absent (or its referent was GC'd).
     *
     * @return the existing live value if present, otherwise the newly inserted value
     */
    public V putIfAbsent(K key, V value) {
        while (true) {
            SoftReference<V> existingRef = map.get(key);
            if (existingRef != null) {
                V existing = existingRef.get();
                if (existing != null) {
                    return existing;
                }
                // Stale reference — remove and retry
                map.remove(key, existingRef);
            }
            if (map.size() >= MAX_CAPACITY) {
                evict();
            }
            SoftReference<V> newRef = new SoftReference<>(value);
            SoftReference<V> prev = map.putIfAbsent(key, newRef);
            if (prev == null) {
                return value;
            }
            V prevValue = prev.get();
            if (prevValue != null) {
                return prevValue;
            }
            // Another thread inserted a stale reference — retry
            map.remove(key, prev);
        }
    }

    /**
     * Clear all cached entries.
     */
    public void cleanup() {
        map.clear();
    }

    /**
     * @return approximate number of key entries (some may have GC'd referents)
     */
    public int size() {
        return map.size();
    }

    public Set<K> keySet() {
        return map.keySet();
    }

    /**
     * Evict roughly half the entries when over capacity.
     * Also purges any keys whose soft references have been cleared by GC.
     */
    private void evict() {
        // First pass: remove stale (GC'd) entries
        map.entrySet().removeIf(e -> e.getValue().get() == null);
        // If still over capacity, remove half
        if (map.size() >= MAX_CAPACITY) {
            int toRemove = map.size() / 2;
            int removed = 0;
            for (K key : map.keySet()) {
                if (removed >= toRemove) break;
                map.remove(key);
                removed++;
            }
        }
    }
}
