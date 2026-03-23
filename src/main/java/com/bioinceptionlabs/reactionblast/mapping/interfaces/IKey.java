/*
 * Copyright (C) 2003-2020 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mapping.interfaces;

import java.util.Comparator;
import com.bioinceptionlabs.reactionblast.mapping.container.helper.Key;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public interface IKey extends Comparable<Key>, Comparator<Key> {

    @Override
    int compareTo(Key t);

    @Override
    boolean equals(Object o);

    /**
     *
     * @return
     */
    @Override
    int hashCode();

    /**
     *
     * @return
     */
    @Override
    String toString();

    @Override
    int compare(Key t, Key t1);
}
