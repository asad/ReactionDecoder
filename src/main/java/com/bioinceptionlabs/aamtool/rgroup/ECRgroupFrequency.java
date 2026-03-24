/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.aamtool.rgroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import com.bioinceptionlabs.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, BioInception.
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ECRgroupFrequency {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(ECRgroupFrequency.class);

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        ECRgroupFrequency ecRgroupFrequency = new ECRgroupFrequency(args);
    }

    /**
     *
     * @param args
     */
    public ECRgroupFrequency(String[] args) {
        LOGGER.debug("------------------------------------------------------");
        Map<String, MultiReactionContainer> reactionMap = new TreeMap<>();
        for (String dir : args) {
            File f = new File(dir);
            if (f.isDirectory()) {
                File[] files = f.listFiles();
                for (File ec : files) {
                    String ecNumber = ec.getName();
                    if (ec.isDirectory()) {
                        File[] listReactionFiles = ec.listFiles();
                        for (File reactionFile : listReactionFiles) {
                            String reactionName = reactionFile.getName().split("\\.")[0];
                            try (FileReader fileReader = new FileReader(reactionFile);
                                 MDLRXNV2000Reader mdlrxnV2000Reader = new MDLRXNV2000Reader(fileReader)) {
                                Reaction reaction = mdlrxnV2000Reader.read(new Reaction());
                                if (reactionMap.containsKey(ecNumber)) {
                                    reactionMap.get(ecNumber).addReaction(reaction, reactionName);
                                } else {
                                    MultiReactionContainer r = new MultiReactionContainer(ecNumber);
                                    r.addReaction(reaction, reactionName);
                                    reactionMap.put(ecNumber, r);
                                }
                            } catch (FileNotFoundException ex) {
                                LOGGER.error(SEVERE, null, ex);
                            } catch (CDKException | IOException ex) {
                                LOGGER.error(SEVERE, null, ex);
                            }
                        }
                    }
                }
            }
        }

        LOGGER.debug("Number of EC parsed " + reactionMap.size());

        int ec1Counter = 0;
        int ec2Counter = 0;
        int ec3Counter = 0;
        int ec4Counter = 0;
        int ec5Counter = 0;
        int ec6Counter = 0;

        int ec1_RGroupReactionCounter = 0;
        int ec2_RGroupReactionCounter = 0;
        int ec3_RGroupReactionCounter = 0;
        int ec4_RGroupReactionCounter = 0;
        int ec5_RGroupReactionCounter = 0;
        int ec6_RGroupReactionCounter = 0;

        int no_common_fragment_in_non_r_group = 0;
        int no_common_fragment_in_r_group = 0;
        int no_common_fragment_in_either = 0;

        Map<String, Set<String>> commonCommonMap = new TreeMap<>();
        Map<String, Set<String>> commonDifferenceMap = new TreeMap<>();
        Map<String, Set<String>> commonUnionMap = new TreeMap<>();

        Set<String> r_group_ec = new TreeSet<>();

        for (String ec : reactionMap.keySet()) {
            if (reactionMap.get(ec).isRGroup()) {
                r_group_ec.add(ec);
            }

            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()
                    && !reactionMap.get(ec).isRGroup()) {
                no_common_fragment_in_non_r_group++;
                LOGGER.debug("Processing EC: " + ec
                        + ", R-found: " + reactionMap.get(ec).isRGroup()
                        + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                        + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                        + ", reaction count: " + reactionMap.get(ec).getReactionCount());
            }
            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()
                    && reactionMap.get(ec).isRGroup()) {
                no_common_fragment_in_r_group++;

                LOGGER.debug("Processing EC: " + ec
                        + ", R-found: " + reactionMap.get(ec).isRGroup()
                        + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                        + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                        + ", reaction count: " + reactionMap.get(ec).getReactionCount());
            }

            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()) {
                no_common_fragment_in_either++;
                LOGGER.debug("Processing EC: " + ec
                        + ", R-found: " + reactionMap.get(ec).isRGroup()
                        + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                        + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                        + ", reaction count: " + reactionMap.get(ec).getReactionCount());
            }

            String ec3 = reactionMap.get(ec).getEnzyme1Level() + "." + reactionMap.get(ec).getEnzyme2Level() + "." + reactionMap.get(ec).getEnzyme3Level();
            if (!commonCommonMap.containsKey(ec3)) {
                commonCommonMap.put(ec3, reactionMap.get(ec).getCommonCommonFP());
            } else if (commonCommonMap.containsKey(ec3)) {
                commonCommonMap.get(ec3).retainAll(reactionMap.get(ec).getCommonCommonFP());
            }
            if (!commonDifferenceMap.containsKey(ec3)) {
                commonDifferenceMap.put(ec3, reactionMap.get(ec).getCommonDifferenceFP());
            } else if (commonDifferenceMap.containsKey(ec3)) {
                commonDifferenceMap.get(ec3).retainAll(reactionMap.get(ec).getCommonDifferenceFP());
            }

            if (!commonUnionMap.containsKey(ec3)) {
                commonUnionMap.put(ec3, reactionMap.get(ec).getAllFP());
            } else if (commonUnionMap.containsKey(ec3)) {
                commonUnionMap.get(ec3).retainAll(reactionMap.get(ec).getAllFP());
            }

            switch (reactionMap.get(ec).getEnzyme1Level()) {
                case "1":
                    ec1Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec1_RGroupReactionCounter++;
                    }
                    break;
                case "2":
                    ec2Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec2_RGroupReactionCounter++;
                    }
                    break;
                case "3":
                    ec3Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec3_RGroupReactionCounter++;
                    }
                    break;
                case "4":
                    ec4Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec4_RGroupReactionCounter++;
                    }
                    break;
                case "5":
                    ec5Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec5_RGroupReactionCounter++;
                    }
                    break;
                case "6":
                    ec6Counter++;
                    if (reactionMap.get(ec).isRGroup()) {
                        ec6_RGroupReactionCounter++;
                    }
                    break;
                default:
                    LOGGER.debug("UNKNOW EC CLASS");
                    break;

            }
        }

        int total_EC = ec1Counter + ec2Counter + ec3Counter + ec4Counter + ec5Counter + ec6Counter;
        int total_r_group_counter
                = ec1_RGroupReactionCounter + ec2_RGroupReactionCounter
                + ec3_RGroupReactionCounter + ec4_RGroupReactionCounter
                + ec5_RGroupReactionCounter + ec6_RGroupReactionCounter;

        LOGGER.debug("------------------------------------------------------");
        LOGGER.debug("R-group EC numbers");
        LOGGER.debug(r_group_ec);

        LOGGER.debug("------------------------------------------------------");
        LOGGER.debug("# EC 1 Numbers " + ec1Counter);
        LOGGER.debug("# R-Group EC 1 Numbers " + ec1_RGroupReactionCounter);

        LOGGER.debug("# EC 2 Numbers " + ec2Counter);
        LOGGER.debug("# R-Group EC 2 Numbers " + ec2_RGroupReactionCounter);

        LOGGER.debug("# EC 3 Numbers " + ec3Counter);
        LOGGER.debug("# R-Group EC 3 Numbers " + ec3_RGroupReactionCounter);

        LOGGER.debug("# EC 4 Numbers " + ec4Counter);
        LOGGER.debug("# R-Group EC 4 Numbers " + ec4_RGroupReactionCounter);

        LOGGER.debug("# EC 5 Numbers " + ec5Counter);
        LOGGER.debug("# R-Group EC 5 Numbers " + ec5_RGroupReactionCounter);

        LOGGER.debug("# EC 6 Numbers " + ec6Counter);
        LOGGER.debug("# R-Group EC 6 Numbers " + ec6_RGroupReactionCounter);

        LOGGER.debug("# EC Numbers " + total_EC);
        LOGGER.debug("# R-Group EC Numbers " + total_r_group_counter);
        LOGGER.debug("# EC with No Common Fragment in Non R-Groups " + no_common_fragment_in_non_r_group);
        LOGGER.debug("# EC with No Common Fragment in R-Groups " + no_common_fragment_in_r_group);
        LOGGER.debug("# EC with No Common Fragment in Either R-Groups " + no_common_fragment_in_either);

        LOGGER.debug("------------------------------------------------------");
        /*
        EC Third Level Signature
         */
        int empty_signature = 0;
        int one_signature_common = 0;
        int more_than_one_signature_common = 0;
        int empty_union_common = 0;
        int one_union_common = 0;
        int more_than_one_union_common = 0;

        for (String ec : commonCommonMap.keySet()) {
            if (commonCommonMap.get(ec).isEmpty() && commonDifferenceMap.get(ec).isEmpty()) {
                empty_signature++;
            } else if (!commonCommonMap.get(ec).isEmpty() || !commonDifferenceMap.get(ec).isEmpty()) {
                int size = commonCommonMap.get(ec).size() + commonDifferenceMap.get(ec).size();
                if (size == 1) {
                    one_signature_common++;
                } else if (size > 1) {
                    more_than_one_signature_common++;
                }
            }

            if (commonUnionMap.get(ec).isEmpty()) {
                empty_union_common++;
            } else if (!commonUnionMap.get(ec).isEmpty()) {
                int size = commonUnionMap.get(ec).size();
                if (size == 1) {
                    one_union_common++;
                } else if (size > 1) {
                    more_than_one_union_common++;
                }
            }

        }
        LOGGER.debug("------------------------------------------------------");
        LOGGER.debug("# 3rd level EC Count " + commonCommonMap.size());
        LOGGER.debug("# 3rd level EC with No Common Fragment " + empty_signature);
        LOGGER.debug("# 3rd level EC with One Common Fragment " + one_signature_common);
        LOGGER.debug("# 3rd level EC with more than One Common Fragment " + more_than_one_signature_common);

        LOGGER.debug("# 3rd level EC Count " + commonUnionMap.size());
        LOGGER.debug("# 3rd level EC with No Common Union Fragment " + empty_union_common);
        LOGGER.debug("# 3rd level EC with One Common Union Fragment " + one_union_common);
        LOGGER.debug("# 3rd level EC with more than One Common Union Fragment " + more_than_one_union_common);
        LOGGER.debug("------------------------------------------------------");
    }

}
