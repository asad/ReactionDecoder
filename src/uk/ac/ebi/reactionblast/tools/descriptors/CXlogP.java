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
package uk.ac.ebi.reactionblast.tools.descriptors;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * This XlogP is calculated by Abdullah Kharamann. It uses summation on XlogP
 * vlaues of Atoms Presently Nitrogen group is not handled
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class CXlogP {

    private static final IAtomContainer molecule = null;
    private static final IAtomContainer eminMol = null;
    private static final Logger LOG = getLogger(CXlogP.class.getName());

    private CXlogP() {
    }

    /**
     *
     * @param Mol
     * @param verbose
     * @throws java.lang.Exception
     */
//    public CXlogP(IAtomContainer Mol, boolean verbose) throws Exception {
//        CXlogP.molecule = Mol;
//        // read in force field parameter. Required for energy calculation
//        MyForceFieldConfigurator ffc = new MyForceFieldConfigurator();
//        String mmfile = "cX" + File.separator + "my_mmff94.prm";
//        ffc.setForceFieldConfigurator(mmfile);
//        new CdkMolStandardize(true).standardize4MMFF94(molecule, verbose);
//
//        // do energy calculation
//        MyMMFF94EnergyFunction en = new MyMMFF94EnergyFunction(molecule,
//                ffc.getParameterSet(), verbose);
//        totalEnergy = en.energyFunctionOfAMolecule(molecule, verbose);
//        // minimize
//        ForceField f = new ForceField((AtomContainer) molecule);
//        GVector moleculeCoords = new GVector(3);
//        moleculeCoords.setSize(molecule.getAtomCount() * 3);
//        moleculeCoords.set(ForceFieldTools.getCoordinates3xNVector(molecule));
////        long start = System.currentTimeMillis();
////        use either both or one of the minimization methods. The first (steepestDescentsMinimization)
////        gives good approximation to the energy minimum. The second  (conjugateGradientMinimization)
////        is however superior to the first if the conformation is not too far to the minimum.
////        Thus best to do is, first to perform a steepestDescentsMinimization followed by a
////        conjugateGradientMinimization.
//
//        f.steepestDescentsMinimization(moleculeCoords, en);
//        f.conjugateGradientMinimization(moleculeCoords, en);
//        eminMol = f.getAtomContainer();
////      long stop = System.currentTimeMillis();
////      System.out.print(en.energyFunctionOfAMolecule(f.getAtomContainer(),verbose) + "\n");
////      write out minimized structure in file min2.sdf
////        MDLV2000Writer writer = new MDLV2000Writer(new FileWriter("min2.sdf"));
////        writer.write(f.getAtomContainer());
////        writer.close();
//
//    // to perform conformational analysis by rotating torsion angles in the molecule
//
////        ConformationalAnalysis confAnal = new ConformationalAnalysis(ffc.getParameterSet());
////        IAtomContainerSet molSet = confAnal.rotateTorsionsByTestingForStericClashesOnly(molecule, 10, verbose);
//
//
//    }
//
//    public double getTotalEnergy() {
//        return totalEnergy;
//    }
//
//    public IAtomContainer getEnergyMinimsedMolecule() {
//        return eminMol;
//    }
//    public double getXlogP() throws CDKException {
//
//        MyXlogP xlog = new MyXlogP();
//        double xlogP = Double.valueOf(xlog.calculate(molecule, false).getValue().toString());
//        return xlogP;
//    }
//
//    public Map<IAtom, Double> getAtomXlogPMap() {
//
//        Map<IAtom, Double> map = new HashMap<IAtom, Double>();
//        for (int i = 0; i < molecule.getAtomCount(); i++) {
//            if(molecule.getAtom(i).getProperty(cX.logP.MyXlogP.xlogpPropName)==null){
//                this.getXlogP();
//            }
//            map.put(molecule.getAtom(i), (Double) molecule.getAtom(i).getProperty(cX.logP.MyXlogP.xlogpPropName));
//
//        }
//    }
}
