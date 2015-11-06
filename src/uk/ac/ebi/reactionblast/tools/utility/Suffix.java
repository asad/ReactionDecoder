/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
/*
 * Suffix.java
 *
 * Created on January 22, 2006, 10:04 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package uk.ac.ebi.reactionblast.tools.utility;

//~--- JDK imports ------------------------------------------------------------
import java.io.IOException;
import java.util.Calendar;
import java.util.GregorianCalendar;
import org.openscience.cdk.math.RandomNumbersTool;

//~--- classes ----------------------------------------------------------------
/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public final class Suffix {

    private static String suffix = "";
    private static String timeSuffix = "";
    private static String randonNumberSuffix = "";
    private static Suffix ref = null;

    //~--- constructors -------------------------------------------------------
    protected Suffix() throws IOException {
        Calendar cal = new GregorianCalendar();
        int ms = cal.get(Calendar.YEAR);
        timeSuffix = String.valueOf(ms);
        ms = cal.get(Calendar.MONTH);
        timeSuffix = timeSuffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.DATE);
        timeSuffix = timeSuffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.HOUR);
        timeSuffix = timeSuffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.MINUTE);
        timeSuffix = timeSuffix.concat(String.valueOf(ms));
        ms = cal.get(Calendar.MILLISECOND);
        timeSuffix = timeSuffix.concat(String.valueOf(ms));

        randonNumberSuffix = String.valueOf(RandomNumbersTool.randomInt(1, 1000));
        suffix = timeSuffix + randonNumberSuffix;
    }

    //~--- get methods --------------------------------------------------------
    /**
     * Creates a new instance of Suffix
     *
     * @return
     * @throws IOException
     */
    public synchronized static Suffix getInstance() throws IOException {
        if (ref == null) {

            // it's ok, we can call this constructor
            ref = new Suffix();
        }

        return ref;
    }

    /**
     *
     * @return time + randomnumber
     */
    public String getSuffix() {
        return suffix;
    }

    /**
     *
     * @return
     */
    public String getTimeSuffix() {
        return timeSuffix;
    }

    /**
     *
     * @return
     */
    public String getRandonNumberSuffix() {
        return randonNumberSuffix;
    }
}
