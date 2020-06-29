/*
 * Copyright (c) 2018-2020. BioInception Labs Pvt. Ltd.
 */
package uk.ac.ebi.reactionblast.mapping.cache;

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
