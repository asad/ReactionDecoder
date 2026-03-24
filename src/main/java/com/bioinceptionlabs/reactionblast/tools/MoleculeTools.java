/*
 * MoleculeTools - consolidated small utility classes.
 * Merged from: Suffix, EBIDoubleUtility, BasicDebugger, AtomContainerSetComparator, EBIMolSplitter, ExtReactionManipulatorTool, ValencyCalculator
 */
package com.bioinceptionlabs.reactionblast.tools;

import com.bioinceptionlabs.reactionblast.fingerprints.IPatternFingerprinter;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.Feature;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter.IFeature;
import com.bioinceptionlabs.reactionblast.fingerprints.PatternFingerprinter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Comparator;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IElectronContainer;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
import org.openscience.smsd.ExtAtomContainerManipulator;
import static com.bioinceptionlabs.reactionblast.mechanism.BondChange.convertBondOrder;
import static java.lang.String.valueOf;
import static java.lang.System.getProperty;
import static java.util.Calendar.DATE;
import static java.util.Calendar.HOUR;
import static java.util.Calendar.MILLISECOND;
import static java.util.Calendar.MINUTE;
import static java.util.Calendar.MONTH;
import static java.util.Calendar.YEAR;
import static java.util.Collections.unmodifiableMap;
import static java.util.logging.Level.WARNING;
import static org.openscience.cdk.CDKConstants.UNSET;
import static org.openscience.cdk.CDKConstants.VISITED;
import static org.openscience.cdk.config.Isotopes.getInstance;
import static org.openscience.cdk.graph.PathTools.breadthFirstSearch;
import static org.openscience.cdk.math.RandomNumbersTool.randomInt;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getSingleBondEquivalentSum;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getElementCount;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getGroup;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol;

/**
 * Container class for miscellaneous molecule utility operations.
 */
public final class MoleculeTools {


    //~--- classes ----------------------------------------------------------------
    /**
     *
     * @contact Syed Asad Rahman, BioInception
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class Suffix {

        private static String suffix = "";
        private static String timeSuffix = "";
        private static String randonNumberSuffix = "";
        private static Suffix ref = null;
        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(Suffix.class);

        //~--- get methods --------------------------------------------------------
        /**
         * Creates a new instance of Suffix
         *
         * @return
         * @throws IOException
         */
        public static Suffix getInstance() throws IOException {
            if (ref == null) {

                // it's ok, we can call this constructor
                ref = new Suffix();
            }

            return ref;
        }

        //~--- constructors -------------------------------------------------------

        /**
         *
         * @throws IOException
         */

        protected Suffix() throws IOException {
            Calendar cal = new GregorianCalendar();
            int ms = cal.get(YEAR);
            timeSuffix = valueOf(ms);
            ms = cal.get(MONTH);
            timeSuffix = timeSuffix.concat(valueOf(ms));
            ms = cal.get(DATE);
            timeSuffix = timeSuffix.concat(valueOf(ms));
            ms = cal.get(HOUR);
            timeSuffix = timeSuffix.concat(valueOf(ms));
            ms = cal.get(MINUTE);
            timeSuffix = timeSuffix.concat(valueOf(ms));
            ms = cal.get(MILLISECOND);
            timeSuffix = timeSuffix.concat(valueOf(ms));

            randonNumberSuffix = valueOf(randomInt(1, 1000));
            suffix = timeSuffix + randonNumberSuffix;
        }

        /**
         *
         * @return time + randomnumber
         */
        public String getSuffix() {
            return suffix;
        }

        /**
         *
         * @return
         */
        public String getTimeSuffix() {
            return timeSuffix;
        }

        /**
         *
         * @return
         */
        public String getRandonNumberSuffix() {
            return randonNumberSuffix;
        }
    }



    /**
     *
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class EBIDoubleUtility {

        private static final long serialVersionUID = 7683452581122892189L;
        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(EBIDoubleUtility.class);

        /**
         *
         * @param val1
         * @param val2
         * @return fused array val1+val2
         * @throws CDKException
         */
        public static double[] append(double[] val1, double[] val2) throws CDKException {

            double[] feature = null;

            if (val1.length > 0 && val2.length > 0) {
                feature = new double[val1.length + val2.length];

                int index = 0;
                for (int i = 0; i < val1.length; i++) {
                    feature[index++] = val1[i];
                }

                for (int j = 0; j < val2.length; j++) {
                    feature[index++] = val2[j];
                }

            } else {
                throw new CDKException("Index < 0: ");
            }

            return feature;

        }

        /**
         *
         * @param val1
         * @param val2
         * @return fused array val1+val2
         * @throws CDKException
         */
        public static IPatternFingerprinter Union(IPatternFingerprinter val1, IPatternFingerprinter val2) throws CDKException {
            PatternFingerprinter patternFingerprinter = new PatternFingerprinter(val1.getFingerprintSize() + val2.getFingerprintSize());

            if (val1.getFingerprintSize() > 0 && val2.getFingerprintSize() > 0) {

                for (int i = 0; i < val1.getFeatureCount(); i++) {
                    IFeature feature = val1.getFeature(i);
                    patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
                }

                for (int j = 0; j < val2.getFeatureCount(); j++) {
                    IFeature feature = val2.getFeature(j);
                    patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));

                }

            } else {
                throw new CDKException("Index < 0: ");
            }
            return patternFingerprinter;
        }

        /**
         *
         * @param val1
         * @param val2
         * @param val3
         * @return fused array val1+val2+val3
         * @throws CDKException
         */
        public static double[] append(double[] val1, double[] val2, double[] val3) throws CDKException {

            double[] feature = null;

            if (val1.length > 0 && val2.length > 0 && val3.length > 0) {
                feature = new double[val1.length + val2.length + val3.length];

                int index = 0;
                for (int i = 0; i < val1.length; i++) {
                    feature[index++] = val1[i];
                }

                for (int j = 0; j < val2.length; j++) {
                    feature[index++] = val2[j];
                }

                for (int k = 0; k < val3.length; k++) {
                    feature[index++] = val3[k];
                }

            } else {
                throw new CDKException("Index < 0: ");
            }

            return feature;

        }

        /**
         *
         * @param val1
         * @param val2
         * @param val3
         * @return fused array val1+val2+val3
         * @throws CDKException
         */
        public static IPatternFingerprinter Union(IPatternFingerprinter val1, IPatternFingerprinter val2, IPatternFingerprinter val3) throws CDKException {

            PatternFingerprinter patternFingerprinter = new PatternFingerprinter(val1.getFingerprintSize() + val2.getFingerprintSize() + val3.getFingerprintSize());

            if (val1.getFingerprintSize() > 0 && val2.getFingerprintSize() > 0) {

                for (int i = 0; i < val1.getFeatureCount(); i++) {
                    IFeature feature = val1.getFeature(i);
                    patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
                }

                for (int j = 0; j < val2.getFeatureCount(); j++) {
                    IFeature feature = val2.getFeature(j);
                    patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
                }

                for (int k = 0; k < val3.getFeatureCount(); k++) {
                    IFeature feature = val3.getFeature(k);
                    patternFingerprinter.add(new Feature(feature.getPattern(), feature.getWeight()));
                }

            } else {
                throw new CDKException("Index < 0: ");
            }
            return patternFingerprinter;

        }

        /**
         *
         * @param val1
         * @param val2
         * @return is val1 contained in val2
         * @throws CDKException
         */
        public static boolean isSubset(double[] val1, double[] val2) throws CDKException {

            boolean flag = true;

            if (val1.length > 0 && val2.length > 0) {

                for (int i = 0; i < val1.length; i++) {
                    if (val1[i] > val2[i]) {
                        flag = false;
                        break;
                    }
                }

            } else {
                throw new CDKException("Index <0: ");
            }

            return flag;

        }

        /**
         *
         * @param val1
         * @param val2
         * @return is val2 contained in val1
         * @throws CDKException
         */
        public static boolean isSuperset(double[] val1, double[] val2) throws CDKException {

            boolean flag = true;

            if (val1.length > 0 && val2.length > 0) {

                for (int i = 0; i < val1.length; i++) {
                    if (val1[i] < val2[i]) {
                        flag = false;
                        break;
                    }
                }

            } else {
                throw new CDKException("Index <0: ");
            }

            return flag;
        }

        private EBIDoubleUtility() {
        }
    }



    /**
     *
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     *
     */
    public static abstract class BasicDebugger {

        private static final ILoggingTool LOGGER
                = LoggingToolFactory.createLoggingTool(BasicDebugger.class);

        protected static final String NEW_LINE = getProperty("line.separator");

        /**
         *
         * @param mappings
         */
        public static void printAtomAtomMapping(Map<IAtom, IAtom> mappings) {
            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            mappings.entrySet().stream().map((m) -> {
                sb.append("e:").append(m.getKey().getID()).append(NEW_LINE);
                return m;
            }).forEach((m) -> {
                sb.append("p:").append(m.getValue().getID()).append(NEW_LINE);
            });
            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param reaction
         */
        protected static void printReaction(IReaction reaction) {
            IAtomContainerSet Educt = reaction.getReactants();
            IAtomContainerSet Product = reaction.getProducts();

            StringBuilder sb = new StringBuilder();
            sb.append("*******************************").append(NEW_LINE);
            sb.append("Educt Mol Count: ").append(Educt.getAtomContainerCount()).append(NEW_LINE);
            sb.append("*******************************").append(NEW_LINE);

            for (int j = 0; j < Educt.getAtomContainerCount(); j++) {

                IAtomContainer M = Educt.getAtomContainer(j);
                sb.append("Mol ID: ").append(M.getID()).append(NEW_LINE);
                sb.append("SingleElectron: ").append(M.getSingleElectronCount()).append(NEW_LINE);
                sb.append("Stoic: ").append(reaction.getReactantCoefficient(M)).append(NEW_LINE);
                sb.append("Split Mol Atom Count: ").append(M.getAtomCount()).append(NEW_LINE);
                appendAtoms(sb, M);
            }

            sb.append("*******************************").append(NEW_LINE);
            sb.append("Product Mol Count: ").append(Product.getAtomContainerCount()).append(NEW_LINE);
            sb.append("*******************************").append(NEW_LINE);

            for (int j = 0; j < Product.getAtomContainerCount(); j++) {

                IAtomContainer M = Product.getAtomContainer(j);
                sb.append("Mol ID: ").append(M.getID()).append(NEW_LINE);
                sb.append("SingleElectron: ").append(M.getSingleElectronCount()).append(NEW_LINE);
                sb.append("Stoic: ").append(reaction.getProductCoefficient(M)).append(NEW_LINE);
                sb.append("Split Mol Atom Count: ").append(M.getAtomCount()).append(NEW_LINE);
                appendAtoms(sb, M);

            }

            sb.append(NEW_LINE).append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%").append(NEW_LINE).append(NEW_LINE);
            LOGGER.debug(sb.toString());
        }

        private static void appendAtoms(StringBuilder sb, IAtomContainer mol) {
            sb.append("Atom: ");
            for (IAtom a : mol.atoms()) {
                sb.append(a.getSymbol());
                if (a.getID() != null) {
                    sb.append("[").append(a.getID()).append("]");
                }
            }
            sb.append(NEW_LINE).append(NEW_LINE);
        }

        /**
         * Print Atoms in molecules
         *
         * @param mol
         */
        protected static void printAtoms(IAtomContainer mol) {
            StringBuilder sb = new StringBuilder();
            sb.append("Atom: ");
            for (IAtom a : mol.atoms()) {

                sb.append(a.getSymbol());
                if (a.getID() != null) {
                    sb.append("[").append(a.getID()).append("]");
                }

            }
            LOGGER.debug(sb.toString());
        }

        /**
         * Prints atoms in molecules
         *
         * @param molecule
         */
        protected static void printMolecule(IAtomContainer molecule) {

            StringBuilder sb = new StringBuilder();
            sb.append("AtomContainer ").append(molecule.getID()).append(": ").append(molecule.getAtomCount()).append(NEW_LINE);

            for (int i = 0; i < molecule.getAtomCount(); i++) {

                sb.append(molecule.getAtom(i).getSymbol()).append(" : ").append(molecule.getAtom(i).getID()).append(",  ");
            }

            LOGGER.debug(sb.toString());

        }
    }



    /**
     *
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     *
     * GraphAtomContainer Comparator
     */
    public static class AtomContainerSetComparator implements Comparator<IAtomContainer> {

        /**
         * Configure LoggingTool
         */
        private final ILoggingTool LOGGER
                = createLoggingTool(AtomContainerSetComparator.class);

        /**
         * Creates a new instance of AtomContainerComparator
         */
        public AtomContainerSetComparator() {
        }

        /*
         * <p>Compares two IAtomContainers for order with the following criteria with decreasing priority:</p>
         * <ul>
         *   <li>Compare atom count
         *   <li>Compare molecular weight (heavy atoms only)
         *   <li>Compare bond count
         *   <li>Compare sum of bond orders (heavy atoms only)
         * </ul>
         * <p>If no difference can be found with the above criteria, the IAtomContainers are
         * considered equal.</p>
         * <p>Returns a negative integer, zero, or a positive integer as the first argument is less than,
         * equal to, or greater than the second.</p>
         * <p>This method is null safe.</p>
         *
         * @param o1 the first IAtomContainer
         * @param o2 the second IAtomContainer
         * @return a negative integer, zero, or a positive integer as the first argument is less than, equal
         *         to, or greater than the second.
         */
        /**
         *
         * @param o1
         * @param o2
         * @return
         */
        @Override
        public int compare(IAtomContainer o1, IAtomContainer o2) {
            // Check for nulls
            if (o1 == null && o2 == null) {
                return 0;
            }
            if (o1 == null) {
                return -1;
            }
            if (o2 == null) {
                return 1;
            }

            // Check for correct instances
            if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
                return 0;
            }
            if (!(o1 instanceof IAtomContainer)) {
                return -1;
            }
            if (!(o2 instanceof IAtomContainer)) {
                return 1;
            }

            // Check for correct instances
            if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
                return 0;
            }
            if (!(o1 instanceof IAtomContainer)) {
                return -1;
            }
            if (!(o2 instanceof IAtomContainer)) {
                return 1;
            }

            IAtomContainer atomContainer1 = o1;
            IAtomContainer atomContainer2 = o2;

            // 1. Compare atom count
            if (atomContainer1.getAtomCount() > atomContainer2.getAtomCount()) {
                return -1;
            } else if (atomContainer1.getAtomCount() < atomContainer2.getAtomCount()) {
                return 1;
            } else {
                // 2. Atom count equal, compare molecular weight (heavy atoms only)
                double mw1;
                double mw2;
                try {
                    mw1 = getMolecularWeight(atomContainer1);
                    mw2 = getMolecularWeight(atomContainer2);
                } catch (CDKException e) {
                    LOGGER.warn("Exception in molecular mass calculation.");
                    return 0;
                }
                if (mw1 > mw2) {
                    return -1;
                } else if (mw1 < mw2) {
                    return 1;
                } else {
                    // 3. Molecular weight equal, compare bond count
                    if (atomContainer1.getBondCount() > atomContainer2.getBondCount()) {
                        return -1;
                    } else if (atomContainer1.getBondCount() < atomContainer2.getBondCount()) {
                        return 1;
                    } else {
                        // 4. Bond count equal, compare sum of bond orders (heavy atoms only)
                        double bondOrderSum1 = getSingleBondEquivalentSum(atomContainer1);
                        double bondOrderSum2 = getSingleBondEquivalentSum(atomContainer2);
                        if (bondOrderSum1 > bondOrderSum2) {
                            return -1;
                        } else if (bondOrderSum1 < bondOrderSum2) {
                            return 1;
                        }
                    }

                }
            }
            // AtomContainers are equal in terms of this comparator
            return 0;
        }

        /**
         * Returns the molecular weight (exact mass) of the major isotopes of all
         * heavy atoms of the given IAtomContainer.
         *
         * @param atomContainer an IAtomContainer to calculate the molecular weight
         * for
         * @throws org.openscience.cdk.exception.CDKException if an error occurs
         * with the IsotopeFactory
         * @return the molecular weight (exact mass) of the major isotopes of all
         * heavy atoms of the given IAtomContainer
         */
        private double getMolecularWeight(IAtomContainer atomContainer) throws CDKException {
            double mw = 0.0;
            try {
                for (IAtom atom : atomContainer.atoms()) {
                    if (!atom.getSymbol().equals("H") && !atom.getSymbol().equals("R")) {
                        try {
                            try {
                                IIsotope majorIsotope = getInstance().getMajorIsotope(atom.getSymbol());
                                mw += majorIsotope.getExactMass();
                            } catch (NullPointerException e) {
                                mw += getInstance().getMajorIsotope("Ra").getExactMass();
                                LOGGER.warn("Isotopes not defined in the CDK " + atom.getSymbol());
                            }
                        } catch (IOException e) {
                            LOGGER.warn("Molecular weight calculation failed for atom " + atom.getSymbol());
                        }
                    } else if (atom.getSymbol().equals("R")) {
                        mw += getInstance().getMajorIsotope("C").getExactMass();
                    }
                }
            } catch (IOException e) {
                LOGGER.warn("Molecular weight calculation failed for atleast one atom ");
            }
            return mw;
        }
    }



    /**
     *
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class EBIMolSplitter {

        /**
         * Check whether a set of atoms in an atomcontainer is connected
         *
         * @param atomContainer The GraphAtomContainer to be check for connectedness
         * @return true if the GraphAtomContainer is connected
         */
        public static boolean isConnected(IAtomContainer atomContainer) {
            boolean flag = false;

            IAtomContainer ac = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            IAtom atom = null;
            IAtomContainer molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            List<IAtom> sphere = new ArrayList<>();
            for (int f = 0; f < atomContainer.getAtomCount(); f++) {
                atom = atomContainer.getAtom(f);
                atom.setFlag(VISITED, false);
                ac.addAtom(atomContainer.getAtom(f));
            }

            Iterator<IBond> bonds = atomContainer.bonds().iterator();
            while (bonds.hasNext()) {
                IBond bond = bonds.next();
                bond.setFlag(VISITED, false);
                ac.addBond(bond);
            }
            atom = ac.getAtom(0);
            sphere.add(atom);
            atom.setFlag(VISITED, true);
            breadthFirstSearch(ac, sphere, molecule);
            if (molecule.getAtomCount() == atomContainer.getAtomCount()) {
                flag = true;
            }
            return flag;
        }

        /**
         * Partitions the atoms in an GraphAtomContainer into covalently connected
         * components.
         *
         * @param atomContainer The GraphAtomContainer to be partitioned into
         * connected components, i.e. molecules
         * @return A MoleculeSet.
         *
         *
         */
        public static IAtomContainerSet splitMolecules(IAtomContainer atomContainer) {
            IAtomContainer ac = atomContainer.getBuilder().newInstance(IAtomContainer.class);
            IAtom atom;
            IElectronContainer eContainer;
            IAtomContainer molecule;
            IAtomContainerSet molecules = atomContainer.getBuilder().newInstance(IAtomContainerSet.class);
            List<IAtom> sphere = new ArrayList<>();
            for (int f = 0; f < atomContainer.getAtomCount(); f++) {
                atom = atomContainer.getAtom(f);
                atom.setFlag(VISITED, false);
                ac.addAtom(atom);
            }
            Iterator<IElectronContainer> eContainers = atomContainer.electronContainers().iterator();
            while (eContainers.hasNext()) {
                eContainer = eContainers.next();
                eContainer.setFlag(VISITED, false);
                ac.addElectronContainer(eContainer);
            }
            while (ac.getAtomCount() > 0) {
                atom = ac.getAtom(0);
                molecule = atomContainer.getBuilder().newInstance(IAtomContainer.class);
                sphere.clear();
                sphere.add(atom);
                atom.setFlag(VISITED, true);
                breadthFirstSearch(ac, sphere, molecule);
                molecules.addAtomContainer(molecule);
                ac.remove(molecule);
            }
            return molecules;
        }

        private EBIMolSplitter() {
        }
    }



    /**
     *
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class ExtReactionManipulatorTool extends ReactionManipulator {


        /**
         *
         * @param reaction
         * @return deep clone of the reactions with mol IDs set and reaction ids set
         * plus flags copied
         * @throws CloneNotSupportedException
         */
        public static IReaction deepClone(IReaction reaction) throws CloneNotSupportedException {
            IReaction clone = new Reaction();
            // clone the reactants, products and agents

            for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                IAtomContainer acClone = new AtomContainer(ac).clone();
                /*Set IDs as CDK clone doesn't*/
                for (int i = 0; i < ac.getAtomCount(); i++) {
                    acClone.getAtom(i).setID(ac.getAtom(i).getID());
                }
                acClone.setID(ac.getID());
                acClone.addProperties(ac.getProperties());
                clone.getReactants().addAtomContainer(acClone);
            }

            for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                IAtomContainer acClone = new AtomContainer(ac).clone();
                /*Set IDs as CDK clone doesn't*/
                for (int i = 0; i < ac.getAtomCount(); i++) {
                    acClone.getAtom(i).setID(ac.getAtom(i).getID());
                }
                acClone.setID(ac.getID());
                acClone.addProperties(ac.getProperties());
                clone.getProducts().addAtomContainer(acClone);
            }

            for (IAtomContainer ac : reaction.getAgents().atomContainers()) {
                IAtomContainer acClone = new AtomContainer(ac).clone();
                acClone.setID(ac.getID());
                acClone.addProperties(ac.getProperties());
                clone.getAgents().addAtomContainer(acClone);
            }

            // create a Map of corresponding atoms for molecules (key: original Atom, 
            // value: clone Atom)
            Map<IChemObject, IChemObject> atomatom = new HashMap<>();
            for (int i = 0; i < reaction.getReactants().getAtomContainerCount(); ++i) {
                IAtomContainer mol = reaction.getReactants().getAtomContainer(i);
                IAtomContainer mol2 = clone.getReactants().getAtomContainer(i);
                for (int j = 0; j < mol.getAtomCount(); ++j) {
                    atomatom.put(mol.getAtom(j), mol2.getAtom(j));
                }
            }
            for (int i = 0; i < reaction.getProducts().getAtomContainerCount(); ++i) {
                IAtomContainer mol = reaction.getProducts().getAtomContainer(i);
                IAtomContainer mol2 = clone.getProducts().getAtomContainer(i);
                for (int j = 0; j < mol.getAtomCount(); ++j) {
                    atomatom.put(mol.getAtom(j), mol2.getAtom(j));
                }
            }
            //Add mapping to the clone
            for (IMapping mapping : reaction.mappings()) {
                clone.addMapping(new Mapping(atomatom.get(mapping.getChemObject(0)), atomatom.get(mapping.getChemObject(1))));

            }
            clone.setID(reaction.getID());
            return clone;
        }

        /**
         *
         * @param reaction
         * @return a new mol with explicit Hydrogens
         * @throws CloneNotSupportedException
         */
        public static IReaction addExplicitH(IReaction reaction) throws CloneNotSupportedException {
            IReaction r = reaction.getBuilder().newInstance(IReaction.class);
            for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
                IAtomContainer addExplicitH = ExtAtomContainerManipulator.addExplicitH(ac);
                r.addReactant(addExplicitH, reaction.getReactantCoefficient(ac));
            }

            for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
                IAtomContainer addExplicitH = ExtAtomContainerManipulator.addExplicitH(ac);
                r.addProduct(addExplicitH, reaction.getProductCoefficient(ac));
            }

            r.setDirection(reaction.getDirection());
            r.setID(reaction.getID() == null ? "" : reaction.getID());

            return r;
        }
    }



    /**
     * @refer for valency http://en.wikipedia.org/wiki/Periodic_table_(valence)
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static class ValencyCalculator {

        private static Map<String, Integer> valencElectronMap = null;
        private static boolean isInitialized = false;
        private final static ILoggingTool LOGGER
                = createLoggingTool(ValencyCalculator.class);

        private static void initialize() {
            if (isInitialized) {
                return;
            }
            valencElectronMap = new TreeMap<>();
            for (int i = 1; i < getElementCount(); i++) {
                String symbol = getSymbol(i);
                if (getGroup(symbol) != null
                        && (getGroup(symbol) < 3 || getGroup(symbol) > 12)) {

                    switch (getGroup(symbol)) {
                        case 1:
                            valencElectronMap.put(symbol, 1);
                            break;
                        case 2:
                            valencElectronMap.put(symbol, 2);
                            break;
                        case 13:
                            valencElectronMap.put(symbol, 3);
                            break;
                        case 14:
                            valencElectronMap.put(symbol, 4);
                            break;
                        case 15:
                            valencElectronMap.put(symbol, 5);
                            break;
                        case 16:
                            valencElectronMap.put(symbol, 6);
                            break;
                        case 17:
                            valencElectronMap.put(symbol, 7);
                            break;
                        case 18:
                            valencElectronMap.put(symbol, 8);
                            break;
                        default:
                            valencElectronMap.put(symbol, 0);
                            break;
                    }
                } else {
                    valencElectronMap.put(symbol, 99);
                }

            }
            /* 
             * Metal
             */
            valencElectronMap.put("Sc", 3);
            valencElectronMap.put("Ti", 4);
            valencElectronMap.put("V", 5);
            valencElectronMap.put("Cr", 6);
            valencElectronMap.put("Mn", 4);
            valencElectronMap.put("Ni", 2);
            valencElectronMap.put("Cu", 2);
            valencElectronMap.put("Zn", 2);
            valencElectronMap.put("Fe", 3);
            valencElectronMap.put("Co", 3);
            /* 
             * Generics
             */
            valencElectronMap.put("*", 1);
            valencElectronMap.put("R", 1);
            valencElectronMap.put("A", 1);
            valencElectronMap.put("X", 8);
            valencElectronMap.put("PsH", 1);

            isInitialized = true;
        }

        /**
         * This method calculates the valence of an atom.
         *
         * @param atom The IAtom for which the DescriptorValue is requested
         * @return atomValence The valency for the given atom
         * @throws CDKException
         */
        public static Integer getValenceElectron(IAtom atom) throws CDKException {
            initialize();
            Integer atomValence;
            String symbol = atom.getSymbol();
            if (valencElectronMap.containsKey(symbol)) {
                atomValence = valencElectronMap.get(symbol);
            } else {
                LOGGER.warn(WARNING, "Element {0} not found. Valence assigned 99.", symbol);
                atomValence = 99;
            }
            return atomValence;
        }

        /**
         *
         * @param m
         * @param atom
         * @param skipHydrogen
         * @return
         * @throws CDKException
         */
        public static Integer getFreeValenceElectrons(IAtomContainer m, IAtom atom, boolean skipHydrogen) throws CDKException {
            initialize();
            Integer totalConnectedBondOrder = 0;
            List<IAtom> connectedAtoms = m.getConnectedAtomsList(atom);
            int counterH = 0;
            for (IAtom connAtom : connectedAtoms) {
                if (skipHydrogen && connAtom.getSymbol().equalsIgnoreCase("H")) {
                    counterH++;
                }
                IBond bond = m.getBond(atom, connAtom);
                totalConnectedBondOrder += convertBondOrder(bond);
            }
            Integer charge = Objects.equals(atom.getFormalCharge(), UNSET) ? 0 : atom.getFormalCharge();
            return skipHydrogen ? (getValenceElectron(atom) - totalConnectedBondOrder + counterH - charge)
                    : (getValenceElectron(atom) - totalConnectedBondOrder - charge);
        }

        /**
         * @return Elements
         */
        public static String[] getElements() {
            initialize();
            String[] st = new String[valencElectronMap.size()];
            int i = 0;
            for (String s : valencElectronMap.keySet()) {
                st[i++] = s;
            }
            return st;
        }

        /**
         * @return Element Map
         */
        public static Map<String, Integer> getElementMap() {
            initialize();
            return unmodifiableMap(valencElectronMap);
        }

        /**
         *
         * @return
         */
        public static int getSize() {
            initialize();
            return valencElectronMap.size();
        }

        /**
         *
         * @return
         */
        public static Iterable<String> getKeySet() {
            initialize();
            return valencElectronMap.keySet();
        }

        /**
         *
         * @param key
         * @return the valence
         */
        public static int getValue(String key) {
            initialize();
            return valencElectronMap.get(key);
        }

        private ValencyCalculator() {
        }
    }


}
