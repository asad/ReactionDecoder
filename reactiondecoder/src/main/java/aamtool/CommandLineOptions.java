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
package aamtool;

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
        optionsAAM.addOption("m", "mappings", false, "Report all mappings");
        optionsAAM.addOption("u", "premap", false, "use user defined mappings");
        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Job prefix").withArgName("prefix").create("p"));

        optionsAAM.addOption(
                OptionBuilder.hasArg().withDescription("Output format (TEXT/XML)").withArgName("formatO").create("f"));

        return optionsAAM;
    }

    protected Options createAnnotateOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query type (RXN/SMI)").withArgName("formatQ").create("Q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Task (ANNOTATE)").withArgName("job").create("j"));
        optionsCompare.addOption("u", "premap", false, "use user defined mappings");
        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Job prefix").withArgName("prefix").create("p"));
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Output format (TEXT/XML/BOTH)").withArgName("formatO").create("f"));
        optionsCompare.addOption("x", "patterns", false, "Report all matched molecular pairs (RPAIR type)");
        return optionsCompare;
    }

    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");
        optionsCompare.addOption("u", "premap", false, "use user defined mappings");
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query type (RXN/SMI)").withArgName("formatQ").create("Q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Target type (RXN/SMI)").withArgName("formatT").create("T"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Target").withArgName("target").create("t"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Task (COMPARE)").withArgName("job").create("j"));

        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Job prefix").withArgName("prefix").create("p"));
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Output format (TEXT/XML/BOTH)").withArgName("formatO").create("f"));
        optionsCompare.addOption("x", "patterns", false, "Report all matched molecular pairs (RPAIR type)");
        return optionsCompare;
    }

    protected Options createSimilarityOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Query type (RXN/SMI)").withArgName("formatQ").create("Q"));
        optionsCompare.addOption("u", "premap", false, "use user defined mappings");
        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Input Reaction List").withArgName("queries").create("q"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Task (SIMILARITY)").withArgName("job").create("j"));

        optionsCompare.addOption(
                OptionBuilder.hasArg().withDescription("Job prefix").withArgName("prefix").create("p"));
        return optionsCompare;
    }

    protected Options createSearchOptions() {
        Options options = new Options();
        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption(
                OptionBuilder.hasArg().withDescription("Query Type (RXN/SMI/RID)").withArgName("formatQ").create("Q"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Query").withArgName("query").create("q"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Task (SEARCH)").withArgName("job").create("j"));
        options.addOption(
                OptionBuilder.hasArg().withDescription("Similarity Type (BOND/CENTRE/STRUCTURE)").withArgName("sim").create("s"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Hits (max:100)").withArgName("hit").create("c"));

        options.addOption(
                OptionBuilder.hasArg().withDescription("Output Format (TEXT/XML)").withArgName("formatO").create("f"));
        options.addOption(
                OptionBuilder.hasArg().withDescription("Job Preffix").withArgName("preffix").create("x"));
        return options;
    }
}
