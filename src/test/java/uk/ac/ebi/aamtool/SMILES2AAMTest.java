/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad@ebi.ac.uk>.
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

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import uk.ac.ebi.reactionblast.tools.TestUtility;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class SMILES2AAMTest extends TestUtility {

    public SMILES2AAMTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

//    /**
//     * Test of main method, of class ReactionDecoder.
//     */
//    @Test
//    public void testMain() {
//        System.out.println("main");
//        // TODO code application logic here
//        String RXN_DIR = INFORCHEM_RXN + "rxn" + separator + "balanced" + separator;
//        out.println("RXN File Directory: " + RXN_DIR);
//        /*
//        Instance of SMILES with AAM
//         */
//        SmilesGenerator smilesAAM = generic().withAtomClasses();
//        File dir = new File(RXN_DIR);
//        File[] files = dir.listFiles();
//        for (File file : files) {
//            //System.out.println("Iterating: " + file);
//            IReaction readReaction;
//            try {
//                readReaction = readReactionFile(file.getName().split("\\.")[0], RXN_DIR, true, false);
//                if (readReaction != null) {
//                    IReaction mapReaction = mapReaction(readReaction, true);
//                    out.println(" Mapped Reaction SMILES for Reaction " + file.getName() + ": "
//                            + smilesAAM.createReactionSMILES(mapReaction));
//                    assertEquals(true, mapReaction.getMappingCount() > 0);
//                }
//            } catch (CDKException ex) {
//                getLogger(SMILES2AAMTest.class.getName()).log(SEVERE, null, ex);
//            } catch (Exception ex) {
//                getLogger(SMILES2AAMTest.class.getName()).log(SEVERE, null, ex);
//            }
//        }
//    }

}
