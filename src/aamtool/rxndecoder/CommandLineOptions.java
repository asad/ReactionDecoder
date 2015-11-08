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
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CommandLineOptions {

    public CommandLineOptions() {
    }

    protected Options createAAMOptions() {
        Options optionsAAM = new Options();
        optionsAAM.addOption("h", "help", false, "Help page for command usage");

        optionsAAM.addOption(
                OptionBuilder.create("Q"));

        optionsAAM.addOption(
                OptionBuilder.create("q"));

        optionsAAM.addOption(
                OptionBuilder.create("j"));

        optionsAAM.addOption("g", "image", false, "create png of the mapping");
        optionsAAM.addOption("m", "mappings", false, "Report all Mappings");
        optionsAAM.addOption("p", "patterns", false, "Report all Mol Mol Pair (RPAIR type)");

        optionsAAM.addOption(
                OptionBuilder.create("f"));
        optionsAAM.addOption(
                OptionBuilder.create("x"));
        return optionsAAM;
    }

    protected Options createTransformationOptions() {
        Options options = new Options();
        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption(
                OptionBuilder.create("Q"));

        options.addOption(
                OptionBuilder.create("q"));

        options.addOption(
                OptionBuilder.create("j"));
        options.addOption(
                OptionBuilder.create("c"));

        options.addOption("r", "recursive", false, "Resursive Matching Filter");
        options.addOption(
                OptionBuilder.create("x"));
        options.addOption(
                OptionBuilder.create("f"));

        return options;
    }

    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(
                OptionBuilder.create("Q"));

        optionsCompare.addOption(
                OptionBuilder.create("q"));

        optionsCompare.addOption(
                OptionBuilder.create("T"));

        optionsCompare.addOption(
                OptionBuilder.create("t"));

        optionsCompare.addOption(
                OptionBuilder.create("j"));

        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("m", "mappings", false, "Report all Mappings");
        optionsCompare.addOption(
                OptionBuilder.create("x"));
        optionsCompare.addOption(
                OptionBuilder.create("f"));

        return optionsCompare;
    }
    private static final Logger LOG = Logger.getLogger(CommandLineOptions.class.getName());
}
