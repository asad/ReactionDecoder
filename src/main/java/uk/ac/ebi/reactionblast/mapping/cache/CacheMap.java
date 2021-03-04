/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package uk.ac.ebi.reactionblast.mapping.cache;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 * @param <K>
 * @param <V>
 */
public class CacheMap<K, V> implements Cache<K, V> {

    private final long timeToLive;
    private final HashMap<K, V> cacheMap;

    protected class CrunchifyCacheObject {

        public long lastAccessed = System.currentTimeMillis();
        public String value;

        protected CrunchifyCacheObject(String value) {
            this.value = value;
        }
    }

    /**
     * Test with crunchifyTimeToLive = 200 seconds
     *
     * crunchifyTimerInterval= 500 seconds
     *
     * maxItems = 6
     *
     * MappingCacheMap<String, String>
     * cache = new MappingCacheMap<String, String>(200, 500, 6);
     *
     *
     * @param timeToLive
     * @param timeInterval
     * @param max
     */
    public CacheMap(long timeToLive, final long timeInterval, int max) {
        this.timeToLive = timeToLive * 2000;

        cacheMap = new HashMap<>(max);

        if (timeToLive > 0 && timeInterval > 0) {

            Thread t;
            t = new Thread(() -> {
                while (true) {
                    try {
                        Thread.sleep(timeInterval * 1000);
                    } catch (InterruptedException ex) {
                    }

                }
            });

            t.setDaemon(true);
            t.start();
        }
    }

    // PUT method
    @Override
    public void put(K key, V value) {
        synchronized (cacheMap) {
            cacheMap.put(key, value);
        }
    }

    // GET method
    @SuppressWarnings("unchecked")
    @Override
    public V get(K key) {
        synchronized (cacheMap) {
            CrunchifyCacheObject c = (CrunchifyCacheObject) cacheMap.get(key);

            if (c == null) {
                return null;
            } else {
                c.lastAccessed = System.currentTimeMillis();
                return (V) c.value;
            }
        }
    }

    // REMOVE method
    public void remove(String key) {
        synchronized (cacheMap) {
            cacheMap.remove(key);
        }
    }

    // Get Cache Objects Size()
    public int size() {
        synchronized (cacheMap) {
            return cacheMap.size();
        }
    }

    // CLEANUP method
    public void cleanup() {

        long now = System.currentTimeMillis();
        ArrayList<String> deleteKey = null;

        synchronized (cacheMap) {
            Iterator<?> itr = cacheMap.entrySet().iterator();

            deleteKey = new ArrayList<>((cacheMap.size() / 2) + 1);
            CrunchifyCacheObject c = null;

            while (itr.hasNext()) {
                String key = (String) itr.next();
                c = (CrunchifyCacheObject) ((Entry<?, ?>) itr).getValue();
                if (c != null && (now > (timeToLive + c.lastAccessed))) {
                    deleteKey.add(key);
                }
            }
        }

        for (String key : deleteKey) {
            synchronized (cacheMap) {
                cacheMap.remove(key);
            }

            Thread.yield();
        }
    }
}
