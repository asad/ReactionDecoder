/*
 * Copyright (c) 2018-2026. BioInception Labs Pvt. Ltd.
 */
package com.bioinceptionlabs.reactionblast.mapping;

/**
 *
 * @author Syed Asad Rahman <asad.rahman at bioinceptionlabs.com>
 * @param <K>
 * @param <V>
 */
//Generic Cache Interface
public interface Cache<K, V> {

    public void put(K key, V value);

    public V get(K key);
}
