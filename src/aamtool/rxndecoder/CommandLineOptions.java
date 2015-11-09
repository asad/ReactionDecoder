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
package aamtool.rxndecoder;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.apache.commons.cli.OptionBuilder.create;
import org.apache.commons.cli.Options;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CommandLineOptions {
    private static final Logger LOG = getLogger(CommandLineOptions.class.getName());

    /**
     *
     */
    public CommandLineOptions() {
    }

    /**
     *
     * @return
     */
    protected Options createAAMOptions() {
        Options optionsAAM = new Options();
        optionsAAM.addOption("h", "help", false, "Help page for command usage");

        optionsAAM.addOption(create("Q"));

        optionsAAM.addOption(create("q"));

        optionsAAM.addOption(create("j"));

        optionsAAM.addOption("g", "image", false, "create png of the mapping");
        optionsAAM.addOption("m", "mappings", false, "Report all Mappings");
        optionsAAM.addOption("p", "patterns", false, "Report all Mol Mol Pair (RPAIR type)");

        optionsAAM.addOption(create("f"));
        optionsAAM.addOption(create("x"));
        return optionsAAM;
    }

    /**
     *
     * @return
     */
    protected Options createTransformationOptions() {
        Options options = new Options();
        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption(create("Q"));

        options.addOption(create("q"));

        options.addOption(create("j"));
        options.addOption(create("c"));

        options.addOption("r", "recursive", false, "Resursive Matching Filter");
        options.addOption(create("x"));
        options.addOption(create("f"));

        return options;
    }

    /**
     *
     * @return
     */
    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(create("Q"));

        optionsCompare.addOption(create("q"));

        optionsCompare.addOption(create("T"));

        optionsCompare.addOption(create("t"));

        optionsCompare.addOption(create("j"));

        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("m", "mappings", false, "Report all Mappings");
        optionsCompare.addOption(create("x"));
        optionsCompare.addOption(create("f"));

        return optionsCompare;
    }
}
