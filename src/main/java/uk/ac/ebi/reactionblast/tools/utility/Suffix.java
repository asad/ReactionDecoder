/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
import static java.lang.String.valueOf;
import java.util.Calendar;
import static java.util.Calendar.DATE;
import static java.util.Calendar.HOUR;
import static java.util.Calendar.MILLISECOND;
import static java.util.Calendar.MINUTE;
import static java.util.Calendar.MONTH;
import static java.util.Calendar.YEAR;
import java.util.GregorianCalendar;
import static org.openscience.cdk.math.RandomNumbersTool.randomInt;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

//~--- classes ----------------------------------------------------------------
/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Suffix {

    private static String suffix = "";
    private static String timeSuffix = "";
    private static String randonNumberSuffix = "";
    private static Suffix ref = null;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Suffix.class);

    //~--- get methods --------------------------------------------------------
    /**
     * Creates a new instance of Suffix
     *
     * @return
     * @throws IOException
     */
    public static synchronized Suffix getInstance() throws IOException {
        if (ref == null) {
            
            // it's ok, we can call this constructor
            ref = new Suffix();
        }
        
        return ref;
    }

    //~--- constructors -------------------------------------------------------

    /**
     *
     * @throws IOException
     */

    protected Suffix() throws IOException {
        Calendar cal = new GregorianCalendar();
        int ms = cal.get(YEAR);
        timeSuffix = valueOf(ms);
        ms = cal.get(MONTH);
        timeSuffix = timeSuffix.concat(valueOf(ms));
        ms = cal.get(DATE);
        timeSuffix = timeSuffix.concat(valueOf(ms));
        ms = cal.get(HOUR);
        timeSuffix = timeSuffix.concat(valueOf(ms));
        ms = cal.get(MINUTE);
        timeSuffix = timeSuffix.concat(valueOf(ms));
        ms = cal.get(MILLISECOND);
        timeSuffix = timeSuffix.concat(valueOf(ms));

        randonNumberSuffix = valueOf(randomInt(1, 1000));
        suffix = timeSuffix + randonNumberSuffix;
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
