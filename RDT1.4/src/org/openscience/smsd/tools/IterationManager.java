/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.tools;

import java.io.Serializable;

/**
 * Class that handles execution time of the MCS search.
 *
 *
 * 
 * 
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class IterationManager implements Serializable {

    private static final long serialVersionUID = 396239639826981L;
    private int max;
    private int counter;
    private int coverage;
    private final int limit;

    /**
     * Constructor for storing execution time
     */
    public IterationManager() {
        this(Integer.MAX_VALUE);
    }

    /**
     * Constructor for storing execution time
     *
     * @param maxIteration
     */
    public IterationManager(int maxIteration) {
        this.counter = 0;
        this.coverage = 250;
        this.max = maxIteration;
        this.limit = this.max * this.coverage;
    }

    /**
     * Returns Number of iterations
     *
     * @return Number of iterations
     */
    public synchronized int getCounter() {
        return counter;
    }

    /**
     * increment the counter
     *
     *
     */
    public synchronized void increment() {
        counter++;
    }

    /**
     * decrement the counter
     *
     *
     */
    public synchronized void decrement() {
        counter--;
    }

    public synchronized boolean isMaxIteration() {
        return getCounter() > limit;
    }

    /**
     * @return the coverage
     */
    public synchronized int getCoverage() {
        return coverage;
    }

    /**
     * @param coverage the coverage to set
     */
    public synchronized void setCoverage(int coverage) {
        this.coverage = coverage;
    }

    /**
     * Returns max allowed iterations (upper limit)
     *
     * @return
     */
    public int getIterationLimit() {
        return limit;
    }
}
