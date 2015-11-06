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
                OptionBuilder.hasArg().withDescription("Query Type (RXN/SMI)").withArgName("formatQ").create("Q"));

        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Task (AAM)").withArgName("job").create("j"));

        optionsAAM.addOption("g", "image", false, "create png of the mapping");
        optionsAAM.addOption("m", "mappings", false, "Report all Mappings");
        optionsAAM.addOption("p", "patterns", false, "Report all Mol Mol Pair (RPAIR type)");

        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Output Format (TEXT/XML)").withArgName("formatO").create("f"));
        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Job Preffix").withArgName("preffix").create("x"));
        return optionsAAM;
    }

    protected Options createTransformationOptions() {
        Options options = new Options();
        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption(
                OptionBuilder.hasArg().withDescription("Query Type (RXN/SMI/MOL2/SDF/CML)").withArgName("formatQ").create("Q"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Task (TRANSFORM)").withArgName("job").create("j"));
        options.addOption(
                OptionBuilder.hasArg().withDescription("Hits (max:100)").withArgName("hit").create("c"));

        options.addOption("r", "recursive", false, "Resursive Matching Filter");
        options.addOption(
                OptionBuilder.hasArg().withDescription("Job Preffix").withArgName("preffix").create("x"));
        options.addOption(
                OptionBuilder.hasArg().withDescription("Output Format (TEXT/XML)").withArgName("formatO").create("f"));

        return options;
    }

    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query Type (RXN/SMI/RID)").withArgName("formatQ").create("Q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Target Type (RXN/SMI/RID)").withArgName("formatT").create("T"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Target").withArgName("target").create("t"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Task (COMPARE)").withArgName("job").create("j"));

        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("m", "mappings", false, "Report all Mappings");
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Job Preffix").withArgName("preffix").create("x"));
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Output Format (TEXT/XML/BOTH)").withArgName("formatO").create("f"));

        return optionsCompare;
    }
}
