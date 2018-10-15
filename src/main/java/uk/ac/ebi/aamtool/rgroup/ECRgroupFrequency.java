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
package uk.ac.ebi.aamtool.rgroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.System.out;
import java.util.ArrayList;
import static java.util.Arrays.asList;
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
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class ECRgroupFrequency {

    /**
     *
     */
    protected final static boolean DEBUG = false;

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
        out.println("------------------------------------------------------");
        Map<String, MultiReactionContainer> reactionMap = new TreeMap<>();
        for (String dir : args) {
            File f = new File(dir);
            if (f.isDirectory()) {
                File[] files = f.listFiles();
                //
                if (DEBUG) {
                    List<File> l = new ArrayList<>();
                    l.addAll(asList(files));
                    List<File> subList = l.subList(1, 100);
                    files = subList.toArray(new File[subList.size()]);
                }
                for (File ec : files) {
                    String ecNumber = ec.getName();
                    if (ec.isDirectory()) {
                        File[] listReactionFiles = ec.listFiles();
                        for (File reactionFile : listReactionFiles) {
                            String reactionName = reactionFile.getName().split("\\.")[0];
                            MDLRXNV2000Reader mdlrxnV2000Reader;
                            try {
                                mdlrxnV2000Reader = new MDLRXNV2000Reader(new FileReader(reactionFile));
                                Reaction reaction = mdlrxnV2000Reader.read(new Reaction());
                                mdlrxnV2000Reader.close();
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

        if (DEBUG) {
            out.println("Number of EC parsed " + reactionMap.size());
        }

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
            if (DEBUG) {
                if (reactionMap.get(ec).isRGroup()) {
                    out.println("Processing EC: " + ec
                            + ", R-found: " + reactionMap.get(ec).isRGroup()
                            + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP().size()
                            + ", reaction count: " + reactionMap.get(ec).getReactionCount());
                }
            }

            if (reactionMap.get(ec).isRGroup()) {
                r_group_ec.add(ec);
            }

            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()
                    && !reactionMap.get(ec).isRGroup()) {
                no_common_fragment_in_non_r_group++;
                if (DEBUG) {
                    out.println("Processing EC: " + ec
                            + ", R-found: " + reactionMap.get(ec).isRGroup()
                            + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                            + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                            + ", reaction count: " + reactionMap.get(ec).getReactionCount());
                }
            }
            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()
                    && reactionMap.get(ec).isRGroup()) {
                no_common_fragment_in_r_group++;
                if (DEBUG) {

                    out.println("Processing EC: " + ec
                            + ", R-found: " + reactionMap.get(ec).isRGroup()
                            + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                            + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                            + ", reaction count: " + reactionMap.get(ec).getReactionCount());
                }
            }

            if (reactionMap.get(ec).getCommonCommonFP().isEmpty()
                    && reactionMap.get(ec).getCommonDifferenceFP().isEmpty()) {
                no_common_fragment_in_either++;
                if (DEBUG) {
                    out.println("Processing EC: " + ec
                            + ", R-found: " + reactionMap.get(ec).isRGroup()
                            + ", common fragment: " + reactionMap.get(ec).getCommonCommonFP()
                            + ", difference fragment: " + reactionMap.get(ec).getCommonDifferenceFP()
                            + ", reaction count: " + reactionMap.get(ec).getReactionCount());
                }
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
                    out.println("UNKNOW EC CLASS");
                    break;

            }
        }

        int total_EC = ec1Counter + ec2Counter + ec3Counter + ec4Counter + ec5Counter + ec6Counter;
        int total_r_group_counter
                = ec1_RGroupReactionCounter + ec2_RGroupReactionCounter
                + ec3_RGroupReactionCounter + ec4_RGroupReactionCounter
                + ec5_RGroupReactionCounter + ec6_RGroupReactionCounter;

        out.println("------------------------------------------------------");
        out.println("R-group EC numbers");
        out.println(r_group_ec);

        out.println("------------------------------------------------------");
        out.println("# EC 1 Numbers " + ec1Counter);
        out.println("# R-Group EC 1 Numbers " + ec1_RGroupReactionCounter);

        out.println("# EC 2 Numbers " + ec2Counter);
        out.println("# R-Group EC 2 Numbers " + ec2_RGroupReactionCounter);

        out.println("# EC 3 Numbers " + ec3Counter);
        out.println("# R-Group EC 3 Numbers " + ec3_RGroupReactionCounter);

        out.println("# EC 4 Numbers " + ec4Counter);
        out.println("# R-Group EC 4 Numbers " + ec4_RGroupReactionCounter);

        out.println("# EC 5 Numbers " + ec5Counter);
        out.println("# R-Group EC 5 Numbers " + ec5_RGroupReactionCounter);

        out.println("# EC 6 Numbers " + ec6Counter);
        out.println("# R-Group EC 6 Numbers " + ec6_RGroupReactionCounter);

        out.println("# EC Numbers " + total_EC);
        out.println("# R-Group EC Numbers " + total_r_group_counter);
        out.println("# EC with No Common Fragment in Non R-Groups " + no_common_fragment_in_non_r_group);
        out.println("# EC with No Common Fragment in R-Groups " + no_common_fragment_in_r_group);
        out.println("# EC with No Common Fragment in Either R-Groups " + no_common_fragment_in_either);

        out.println("------------------------------------------------------");
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
            if (DEBUG) {
                if (commonUnionMap.get(ec).isEmpty()) {
                    out.println("EC: " + ec
                            + ", Common Signature: " + commonCommonMap.get(ec)
                            + ", Difference Signature " + commonDifferenceMap.get(ec)
                            + ", Union Common Signature " + commonUnionMap.get(ec));
                }
            }

            //if (commonUnionMap.get(ec).isEmpty()) {
//            System.out.println("EC: " + ec
//                    + ", Common Signature: " + commonCommonMap.get(ec)
//                    + ", Difference Signature " + commonDifferenceMap.get(ec)
//                    + ", Union Common Signature " + commonUnionMap.get(ec));
            //}
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
        out.println("------------------------------------------------------");
        out.println("# 3rd level EC Count " + commonCommonMap.size());
        out.println("# 3rd level EC with No Common Fragment " + empty_signature);
        out.println("# 3rd level EC with One Common Fragment " + one_signature_common);
        out.println("# 3rd level EC with more than One Common Fragment " + more_than_one_signature_common);

        out.println("# 3rd level EC Count " + commonUnionMap.size());
        out.println("# 3rd level EC with No Common Union Fragment " + empty_union_common);
        out.println("# 3rd level EC with One Common Union Fragment " + one_union_common);
        out.println("# 3rd level EC with more than One Common Union Fragment " + more_than_one_union_common);
        out.println("------------------------------------------------------");
    }

}
