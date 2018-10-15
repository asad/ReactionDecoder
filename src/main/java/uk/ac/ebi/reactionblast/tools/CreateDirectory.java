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
 * CreateDirectory.java
 *
 * Created on 13 January 2006, 09:41
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package uk.ac.ebi.reactionblast.tools;

//~--- JDK imports ------------------------------------------------------------
import java.io.File;
import java.io.IOException;
import static java.lang.Thread.sleep;

import uk.ac.ebi.reactionblast.tools.utility.Suffix;
import static uk.ac.ebi.reactionblast.tools.utility.Suffix.getInstance;

//~--- classes ----------------------------------------------------------------
/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class CreateDirectory {

    private static final int MKDIR_RETRY_SLEEP_MILLIS = 10;
    private static String suffix = null;

    /**
     *
     * @return
     */
    public static String getSuffix() {
        return suffix;
    }
    //~--- fields -------------------------------------------------------------
    private File dir;

    //~--- methods ------------------------------------------------------------
    private void createDir() throws IOException {
        if (dir == null) {
            throw new IOException("dir attribute is required");
        }

        if (dir.isFile()) {
            throw new IOException("Unable to create directory as a file "
                    + "already exists with that name: "
                    + dir.getAbsolutePath());
        }

        if (!dir.exists()) {
            boolean result = mkdirs(dir);

            if (!result) {
                String msg
                        = "Directory " + dir.getAbsolutePath()
                        + " creation was not successful for an unknown reason";

                throw new IOException(msg);
            }

            //System.out.println("Created dir: " + dir.getAbsolutePath());
            // path=dir.getAbsolutePath();
        }
    }

    /**
     * Creates a new instance of CreateDirectory
     */
    public void CreateDirectory() {
    }

    private boolean mkdirs(File f) {
        boolean status = true;

        if (!f.mkdirs()) {
            try {
                sleep(MKDIR_RETRY_SLEEP_MILLIS);
                status = f.mkdirs();
            } catch (InterruptedException ex) {
                status = f.mkdirs();
            }
        }

        return status;
    }

    private void setSuffix() throws IOException {
        Suffix init = getInstance();
        suffix = init.getSuffix();
    }
    //~--- get methods --------------------------------------------------------

    //~--- set methods --------------------------------------------------------
    /**
     *
     * @param dirname
     * @param SUFFIX
     * @return
     * @throws java.io.IOException
     */
    public File createDirectory(String dirname, boolean SUFFIX) throws IOException {
        if (SUFFIX) {
            setSuffix();
            String dirSuffix = getSuffix();
            dir = new File(dirname.concat(dirSuffix));
            createDir();
        } else {
            dir = new File(dirname);
            createDir();
        }
        return dir;
    }

    /**
     *
     * @param parentDirName
     * @param childDirName
     * @param SUFFIX
     * @return
     * @throws IOException
     */
    public File createDirectory(String parentDirName, String childDirName, boolean SUFFIX)
            throws IOException {
        if (SUFFIX) {
            setSuffix();
            String dirSuffix = getSuffix();
            childDirName = childDirName.concat(".");
            childDirName = childDirName.concat(dirSuffix);
            dir = new File(parentDirName, childDirName);
            createDir();
        } else {
            dir = new File(parentDirName, childDirName);
            createDir();
        }
        return dir;
    }

}
