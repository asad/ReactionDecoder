/* Copyright (C) 2009-2018  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus1;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class MoleculeHandler {

    private final IAtomContainer atomContainer;
    int bondNumber = 0;
    int atomNumber = 0;
    int startHatom_num = 0;
    private List<IAtom> atomString;
    public List<Integer> intTable = new LinkedList<>();
    public List<String> charTable = new LinkedList<>();
    protected List<Integer> specified_int_tab = new LinkedList<>();

    protected List<Integer> int_tab = new LinkedList<>();
    protected List<String> char_tab = new LinkedList<>();
    private final boolean matchBonds;

    /**
     * Creates a new instance of MolFileReader
     *
     * @param atomContainer
     * @param matchBonds
     */
    public MoleculeHandler(IAtomContainer atomContainer, boolean matchBonds) {
        this.atomContainer = atomContainer;
        this.bondNumber = atomContainer.getBondCount();
        this.atomNumber = atomContainer.getAtomCount();
        this.matchBonds = matchBonds;

        setAtomString();
        setIntConnectionTable();
        setCharConnectionTable();

        enumerate_startHatom_num();

        boolean no_H2 = true;
        if ((atomNumber == 2) && (bondNumber == 1)
                && (atomString.get(0).getSymbol().equals("H"))
                && (atomString.get(1).getSymbol().equals("H"))) {
            no_H2 = false;
        }
        discard_H_bonds(no_H2);

        int_tab_specifier(this.atomContainer.getAtomCount() * this.atomContainer.getBondCount());

    }

    public int getBondNumber() {
        return bondNumber;
    }

    public int indexOf() {
        return atomNumber;
    }

    public List<IAtom> getAtomString() {
        return atomString;
    }

    public int getStartHatom_num() {
        setHydrogenNumber();
        return startHatom_num;
    }

    private void setHydrogenNumber() {
        startHatom_num = getAtomContainer().getAtomCount();

        for (int atom = 0; atom < atomNumber; atom++) {
            if ((getAtomContainer().getAtom(atom).getSymbol()).equals("H")) {
                startHatom_num--;
            }

        }

    }

    private void setAtomString() {
        ArrayList<IAtom> temp = new ArrayList<>();
        for (int atom = 0; atom < atomNumber; atom++) {
            IAtom atomType = getAtomContainer().getAtom(atom);
            temp.add(atomType);
        }

        //System.err.println("In atomContainer: getString(temp) " +temp.size()+ " "+ getString(temp));
        this.atomString = temp;

    }

    public void setIntConnectionTable() {

        IAtomContainer ac = (IAtomContainer) getAtomContainer();

        for (int i = 0; i < bondNumber; i++) {
            IBond bond = ac.getBond(i);
            /*This will fetch the connected ATOM as integer and its Bond order ex: 2 as double, 1 as single */
            // System.out.println(ac.indexOf(bond.getAtom(0))+" "+ac.indexOf(bond.getAtom(1))+" "+(int)bond.getOrder());
            intTable.add((ac.indexOf(bond.getAtom(0)) + 1));//Plus one because Java Indexing is one less
            intTable.add((ac.indexOf(bond.getAtom(1)) + 1));//Plus one because Java indexing is one less
            if (matchBonds) {
                intTable.add((int) bond.getOrder().numeric());
            } else {
                intTable.add(1);
            }

            /*This will fetch the Connected ATOM Symbol*/
//            System.out.println(bond.getAtom(0).getSymbol() + " " + bond.getAtom(1).getSymbol()
//                    + " , bond: " + (int) bond.getOrder().numeric() + " Stored: " + intTable.get(i * 3 + 2));
        }
    }

    public void setCharConnectionTable() {
        IAtomContainer ac = (IAtomContainer) getAtomContainer();
        for (int i = 0; i < bondNumber; i++) {
            IBond bond = ac.getBond(i);
            /*This will fetch the Connected ATOM Symbol*/
            String atom1 = bond.getAtom(0).getSymbol();
            String atom2 = bond.getAtom(1).getSymbol();
            charTable.add(atom1);
            charTable.add(atom2);
        }
    }

    void enumerate_startHatom_num() {

        startHatom_num = atomNumber;
        int a = atomNumber - 1;
        while ((a >= 0) && (atomString.get(a).getSymbol().equals("H"))) {
            startHatom_num--;
            a--;
        }
        //System.out.println("startHatom_num: " + startHatom_num);
    }

    void discard_H_bonds(boolean is_no_H2) {

        int count_bonds = 0;

        if (is_no_H2) {
            for (int x = 0; x < bondNumber; x++) {
                if ((charTable.get(x * 2 + 0).equals("H")) || (charTable.get(x * 2 + 1).equals("H"))) {
                    atomNumber--; // mit jeder gestrichenen H-Bindung veringert sich Atomzahl um 1
                }
                if (!(charTable.get(x * 2 + 0).equals("H")) && !(charTable.get(x * 2 + 1).equals("H"))) {
                    char_tab.add(charTable.get(x * 2 + 0));
                    char_tab.add(charTable.get(x * 2 + 1));
                    int_tab.add(intTable.get(x * 3 + 0));
                    int_tab.add(intTable.get(x * 3 + 1));
                    int_tab.add(intTable.get(x * 3 + 2));
                    count_bonds++;
                }
            }
        } else { //falls es H2 ist:
            for (int x = 0; x < bondNumber; x++) {
                char_tab.add(charTable.get(x * 2 + 0));
                char_tab.add(charTable.get(x * 2 + 1));
                int_tab.add(intTable.get(x * 3 + 0));
                int_tab.add(intTable.get(x * 3 + 1));
                int_tab.add(intTable.get(x * 3 + 2));
                count_bonds++;
            }
        }

        bondNumber = count_bonds;
    }

    /*
     * needed to generate reaction mappings - the vector specified_int_tab is given to
     *
     */
    void int_tab_specifier(int specifier_value) {

        for (int a = 0; a < bondNumber; a++) {
            specified_int_tab.add(int_tab.get(a * 3 + 0) + specifier_value);
            specified_int_tab.add(int_tab.get(a * 3 + 1) + specifier_value);
            specified_int_tab.add(int_tab.get(a * 3 + 2));
        }
    }

    /**
     * @return the atomContainer
     */
    public IAtomContainer getAtomContainer() {
        return atomContainer;
    }
}
