/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.mcsplus;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.BondMatcher;

public class MCSPlusBuilder {

    private IAtomContainer ac1;
    private IAtomContainer ac2;
    private AtomMatcher am;
    private BondMatcher bm;

    public MCSPlusBuilder() {
    }

    public MCSPlusBuilder setQuery(IAtomContainer ac1) {
        this.ac1 = ac1;
        return this;
    }

    public MCSPlusBuilder setTarget(IAtomContainer ac2) {
        this.ac2 = ac2;
        return this;
    }

    public MCSPlusBuilder setAtomMatcher(AtomMatcher am) {
        this.am = am;
        return this;
    }

    public MCSPlusBuilder setBondMatcher(BondMatcher bm) {
        this.bm = bm;
        return this;
    }

    public MCSPlus createMCSPlus() {
        return new MCSPlus(ac1, ac2, am, bm);
    }

}
