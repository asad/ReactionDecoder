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
package uk.ac.ebi.aamtool;

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
        optionsAAM.addOption("Q", "formatQ", true, "Query Type (RXN/SMI)");
        optionsAAM.addOption("q", "query", true, "Query");
        optionsAAM.addOption("j", "job", true, "Task (AAM)");
        optionsAAM.addOption("g", "image", false, "create png of the mapping");
        optionsAAM.addOption("m", "mappings", false, "Report all mappings");
        optionsAAM.addOption("u", "premap", false, "use user defined mappings");
        optionsAAM.addOption("p", "prefix", true, "Job prefix");
        optionsAAM.addOption("f", "formatO", true, "Output format (TEXT/XML)");

        return optionsAAM;
    }

    /**
     *
     * @return
     */
    protected Options createAnnotateOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");
        optionsCompare.addOption("Q", "formatQ", true, "Query Type (RXN/SMI)");
        optionsCompare.addOption("q", "query", true, "Query");
        optionsCompare.addOption("j", "job", true, "Task (ANNOTATE)");
        optionsCompare.addOption("u", "premap", false, "use user defined mappings");
        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("p", "prefix", true, "Job prefix");
        optionsCompare.addOption("f", "formatO", true, "Output format (TEXT/XML/BOTH)");
        optionsCompare.addOption("x", "patterns", false, "Report all matched molecular pairs (RPAIR type)");
        return optionsCompare;
    }

    /**
     *
     * @return
     */
    protected Options createCompareOptions() {
        Options optionsCompare = new Options();
        optionsCompare.addOption("h", "help", false, "Help page for command usage");
        optionsCompare.addOption("u", "premap", false, "use user defined mappings");
        optionsCompare.addOption("Q", "formatQ", true, "Query Type (RXN/SMI)");
        optionsCompare.addOption("q", "query", true, "Query");
        optionsCompare.addOption("T", "formatT", true, "Target Type (RXN/SMI)");
        optionsCompare.addOption("t", "target", true, "Target");
        optionsCompare.addOption("j", "job", true, "Task (COMPARE)");
        optionsCompare.addOption("g", "image", false, "create png of the mapping");
        optionsCompare.addOption("p", "prefix", true, "Job prefix");
        optionsCompare.addOption("f", "formatO", true, "Output format (TEXT/XML/BOTH)");
        optionsCompare.addOption("x", "patterns", false, "Report all matched molecular pairs (RPAIR type)");
        return optionsCompare;
    }

}
