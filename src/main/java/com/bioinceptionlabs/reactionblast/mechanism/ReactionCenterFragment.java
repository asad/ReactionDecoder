/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.mechanism;

import java.io.Serializable;

import com.bioinceptionlabs.reactionblast.mechanism.EnumSubstrateProduct;

/**
 *
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ReactionCenterFragment implements Serializable {

    private static final long serialVersionUID = 9879878799977781L;

    private final String signature;
    private final int level;
    private final EnumSubstrateProduct rpf;

    /**
     *
     * @param signature
     * @param level
     * @param rpf
     */
    public ReactionCenterFragment(String signature, int level, EnumSubstrateProduct rpf) {
        this.signature = signature;
        this.level = level;
        this.rpf = rpf;
    }

    @Override
    public String toString() {
        return "ReactionCenterFragment{" + "signature=" + signature + ", level=" + level + ", rpf=" + rpf + '}';
    }

    /**
     *
     * @return
     */
    public int getLevel() {
        return level;
    }

    /**
     *
     * @return
     */
    public EnumSubstrateProduct getReactantProductInfo() {
        return rpf;
    }

    /**
     *
     * @return
     */
    public String getSignature() {
        return signature;
    }
}
