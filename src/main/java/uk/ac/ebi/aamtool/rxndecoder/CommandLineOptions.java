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
package uk.ac.ebi.aamtool.rxndecoder;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CommandLineOptions {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(CommandLineOptions.class);

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

        optionsAAM.addOption(Option.builder("Q").build());

        optionsAAM.addOption(Option.builder("q").build());

        optionsAAM.addOption(Option.builder("j").build());

        optionsAAM.addOption("g", "image", false, "create png of the mapping");
        optionsAAM.addOption("m", "mappings", false, "Report all Mappings");
        optionsAAM.addOption("p", "patterns", false, "Report all Mol Mol Pair (RPAIR type)");

        optionsAAM.addOption(Option.builder("f").build());
        optionsAAM.addOption(Option.builder("x").build());
        return optionsAAM;
    }

    /**
     *
     * @return
     */
    protected Options createTransformationOptions() {
        Options options = new Options();
        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption(Option.builder("Q").build());

        options.addOption(Option.builder("q").build());

        options.addOption(Option.builder("j").build());
        options.addOption(Option.builder("c").build());

        options.addOption("r", "recursive", false, "Resursive Matching Filter");
        options.addOption(Option.builder("x").build());
        options.addOption(Option.builder("f").build());

        return options;
    }

    /**
     *
     * @return
     */
    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(Option.builder("Q").build());

        optionsCompare.addOption(Option.builder("q").build());

        optionsCompare.addOption(Option.builder("T").build());

        optionsCompare.addOption(Option.builder("t").build());

        optionsCompare.addOption(Option.builder("j").build());

        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("m", "mappings", false, "Report all Mappings");
        optionsCompare.addOption(Option.builder("x").build());
        optionsCompare.addOption(Option.builder("f").build());

        return optionsCompare;
    }
}
