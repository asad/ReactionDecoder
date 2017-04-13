/*
 * Copyright (C) 2007-2017 Syed Asad Rahman <asad@ebi.ac.uk>.
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
package uk.ac.ebi.generic;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import static java.lang.System.err;
import static java.lang.System.out;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.TestUtility;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class UnblancedReactionChecker extends TestUtility {

    private static final boolean DEBUG = false;
    private static final File DIR = new File(RHEA_RXN_DIR);

    private static final Logger LOG = getLogger(UnblancedReactionChecker.class.getName());

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (DIR.isDirectory()) {
            String[] list = DIR.list();
            for (String f : list) {
//                System.out.println("F " + f);
                IReaction rxnReactions;
                try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new FileReader(new File(DIR, f)));) {
                    try {
                        rxnReactions = reader.read(new Reaction());
                        reader.close();
                        rxnReactions.setID(f.split(".rxn")[0]);
                        if (!isBalanced(rxnReactions)) {
                            out.println("Unbalanced Reaction " + f);
                        }
                    } catch (IOException | CDKException ex) {
                        err.println("ERROR in Reading Reaction file " + f + "\n" + ex);
                    }
                } catch (Exception ex) {
                    getLogger(UnblancedReactionChecker.class.getName()).log(SEVERE, null, ex);
                }
            }
        }
    }

    private static boolean isBalanced(IReaction r) {
        Map<String, Integer> atomUniqueCounter1 = new TreeMap<>();
        Map<String, Integer> atomUniqueCounter2 = new TreeMap<>();

        int leftHandAtomCount = 0;
        for (IAtomContainer q : r.getReactants().atomContainers()) {
            for (IAtom a : q.atoms()) {
                if (a.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter1.containsKey(a.getSymbol())) {
                    atomUniqueCounter1.put(a.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter1.get(a.getSymbol()) + 1;
                    atomUniqueCounter1.put(a.getSymbol(), counter);
                }
                leftHandAtomCount++;
            }
            if (DEBUG) {
                try {
                    out.println("Q=mol " + generic().create(q));
                } catch (CDKException ex) {
                    getLogger(ReactionMechanismTool.class.getName()).log(SEVERE, null, ex);
                }
            }
        }

        int rightHandAtomCount = 0;
        for (IAtomContainer t : r.getProducts().atomContainers()) {
            for (IAtom b : t.atoms()) {
                if (b.getSymbol().equals("H")) {
                    continue;
                }
                if (!atomUniqueCounter2.containsKey(b.getSymbol())) {
                    atomUniqueCounter2.put(b.getSymbol(), 1);
                } else {
                    int counter = atomUniqueCounter2.get(b.getSymbol()) + 1;
                    atomUniqueCounter2.put(b.getSymbol(), counter);
                }
                rightHandAtomCount++;
            }
            if (DEBUG) {
                try {
                    out.println("T=mol " + generic().create(t));
                } catch (CDKException ex) {
                    getLogger(ReactionMechanismTool.class.getName()).log(SEVERE, null, ex);
                }
            }
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + leftHandAtomCount);
            out.println("atomUniqueCounter2 " + rightHandAtomCount);
        }

        if (leftHandAtomCount != rightHandAtomCount) {
            err.println();
            err.println("Number of atom(s) on the Left side " + leftHandAtomCount
                    + " =/= Number of atom(s) on the Right side " + rightHandAtomCount);
            atomUniqueCounter1.keySet().stream().map((s) -> {
                if (atomUniqueCounter2.containsKey(s)) {
                    if (atomUniqueCounter1.get(s) != atomUniqueCounter2.get(s).intValue()) {
                        err.println(s + "(" + atomUniqueCounter1.get(s) + ")" + " =/= " + s + "(" + atomUniqueCounter2.get(s) + ")");
                    }
                }
                return s;
            }).filter((s) -> (!atomUniqueCounter2.containsKey(s))).forEach((s) -> {
                err.println(s + "(" + atomUniqueCounter1.get(s) + ")" + " =/= " + s + "(" + 0 + ")");
            });
            atomUniqueCounter2.keySet().stream().filter((s) -> (!atomUniqueCounter1.containsKey(s))).forEach((s) -> {
                err.println(s + "(" + 0 + ")" + " =/= " + s + "(" + atomUniqueCounter2.get(s) + ")");
            });
            return false;
        } else if (!atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet())) {
            err.println();
            err.println("Number of unique atom types(s) on the Left side " + atomUniqueCounter1.size()
                    + " =/= Number of unique atom types(s)on the Right side " + atomUniqueCounter2.size());
            atomUniqueCounter1.keySet().stream().map((s) -> {
                if (atomUniqueCounter2.containsKey(s)) {
                    if (atomUniqueCounter1.get(s) != atomUniqueCounter2.get(s).intValue()) {
                        err.println("Number of reactant Atom: " + s + "(" + atomUniqueCounter1.get(s) + ")" + " =/= Number of product atom: " + s + "(" + atomUniqueCounter2.get(s) + ")");
                    }
                }
                return s;
            }).filter((s) -> (!atomUniqueCounter2.containsKey(s))).forEach((s) -> {
                err.println("Number of reactant Atom: " + s + "(" + atomUniqueCounter1.get(s) + ")" + " =/= Number of product atom: " + s + "(" + 0 + ")");
            });
            atomUniqueCounter2.keySet().stream().filter((s) -> (!atomUniqueCounter1.containsKey(s))).forEach((s) -> {
                err.println("Number of reactant Atom: " + s + "(" + 0 + ")" + " =/= Number of product atom: " + s + "(" + atomUniqueCounter2.get(s) + ")");
            });
            return false;
        }

        if (DEBUG) {
            out.println("atomUniqueCounter1 " + atomUniqueCounter1);
            out.println("atomUniqueCounter2 " + atomUniqueCounter2);
        }
        return atomUniqueCounter1.keySet().equals(atomUniqueCounter2.keySet());
    }

}
