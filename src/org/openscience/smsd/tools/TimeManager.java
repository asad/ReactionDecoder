/* Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
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

import static java.lang.System.currentTimeMillis;
import java.text.SimpleDateFormat;
import static java.util.TimeZone.getTimeZone;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;

/**
 * Class that handles execution time of the MCS search.
 *
 * long diffSeconds = time / 1000; 
 * long diffMinutes = time / (60 * 1000); 
 * long diffHours = time / (60 * 60 * 1000); 
 * long diffDays = time / (24 * 60 * 60 * 1000);
 *
 * 
 * 
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.tools.TimeManagerTest")
public class TimeManager {
    private static final Logger LOG = getLogger(TimeManager.class.getName());

    private final double startTime;
    private final SimpleDateFormat dateFormat;

    /**
     * Constructor for storing execution time
     */
    @TestMethod("testTimeManager")
    public TimeManager() {

        dateFormat = new SimpleDateFormat("HH:mm:ss");
        dateFormat.setTimeZone(getTimeZone("GMT"));
        startTime = currentTimeMillis();
    }

    /**
     * Returns Elapsed Time In Hours
     *
     * @return Elapsed Time In Hours
     */
    @TestMethod("testGetElapsedTimeInHours")
    public synchronized double getElapsedTimeInHours() {
        double currentTime = currentTimeMillis();
        return (currentTime - startTime) / (60 * 60 * 1000);


    }

    /**
     * Returns Elapsed Time In Minutes
     *
     * @return Elapsed Time In Minutes
     */
    @TestMethod("testGetElapsedTimeInMinutes")
    public synchronized double getElapsedTimeInMinutes() {
        double currentTime = currentTimeMillis();
        return (currentTime - startTime) / (60 * 1000);

    }

    /**
     * Return Elapsed Time In Seconds
     *
     * @return Elapsed Time In Seconds
     */
    @TestMethod("testGetElapsedTimeInSeconds")
    public synchronized double getElapsedTimeInSeconds() {
        double currentTime = currentTimeMillis();
        return ((currentTime - startTime) / 1000);

    }

    /**
     * Returns Elapsed Time In Mill Seconds
     *
     * @return Elapsed Time In Mill Seconds
     */
    @TestMethod("testGetElapsedTimeInMilliSeconds")
    public synchronized double getElapsedTimeInMilliSeconds() {
        double currentTime = currentTimeMillis();
        return (currentTime - startTime);

    }
}
