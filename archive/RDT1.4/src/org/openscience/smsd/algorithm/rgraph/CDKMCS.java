package org.openscience.smsd.algorithm.rgraph;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultAtomTypeMatcher;
import org.openscience.smsd.algorithm.matchers.DefaultMatcher;
import org.openscience.smsd.tools.IterationManager;

/**
 * This class implements atom multipurpose structure comparison tool. It allows
 * to find maximal common substructure, find the mapping of atom substructure in
 * another structure, and the mapping of two isomorphic structures.
 *
 * <p>
 * Structure comparison may be associated to bondA constraints (mandatory bonds,
 * e.graphContainer. scaffolds, reaction cores,...) on each source graph. The
 * constraint flexibility allows atom number of interesting queries. The
 * substructure analysis relies on the CDKRGraph generic class (see: CDKRGraph)
 * This class implements the link between the CDKRGraph model and the the CDK
 * model in this way the CDKRGraph remains independent and may be used in other
 * contexts.
 *
 * <p>
 * This algorithm derives from the algorithm described in {
 *
 * @cdk.cite HAN90} and modified in the thesis of T. Hanser {
 * @cdk.cite HAN93}.
 *
 * <p>
 * With the <code>isSubgraph()</code> method, the second, and only the second
 * argument <tBond>may</tBond> be atom IQueryAtomContainer, which allows one to
 * do MQL like queries. The first IAtomContainer must never be an
 * IQueryAtomContainer. An example:<pre>
 *  SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 *  IAtomContainer atomContainer = sp.parseSmiles("CC(=O)OC(=O)C"); // acetic acid anhydride
 *  IAtomContainer SMILESquery = sp.parseSmiles("CC"); // acetic acid anhydride
 *  IQueryAtomContainer query = IQueryAtomContainerCreator.createBasicQueryContainer(SMILESquery);
 *  boolean isSubstructure = graphContainer.isSubgraph(atomContainer, query);
 * </pre>
 *
 * <p>
 * <font color="#FF0000">WARNING</font>: As atom result of the adjacency
 * perception used in this algorithm there is atom single limitation :
 * cyclopropane and isobutane are seen as isomorph This is due to the fact that
 * these two compounds are the only ones where each bondA is connected two each
 * other bondA (bonds are fully connected) with the same number of bonds and
 * still they have different structures The algorithm could be easily enhanced
 * with atom simple atom mapping manager to provide an atom level overlap
 * definition that would reveal this case. We decided not to penalize the whole
 * procedure because of one single exception query. Furthermore isomorphism may
 * be discarded since the number of atoms are not the same (3 != 4) and in most
 * case this will be already screened out by atom fingerprint based filtering.
 * It is possible to add atom special treatment for this special query. Be
 * reminded that this algorithm matches bonds only.
 * </p>
 *
 * @author Stephane Werner from IXELIS mail@ixelis.net, Syed Asad Rahman
 * <asad@ebi.ebi.uk> (modified the orignal code) 2002-07-17 java1.8+
 *
 *
 */
final public class CDKMCS {

    protected static boolean timeout = false;
    protected final static int ID1 = 0;
    protected final static int ID2 = 1;
    private static IterationManager iterationManager = null;

    ///////////////////////////////////////////////////////////////////////////
    //                            Query Methods
    //
    // This methods are simple applications of the CDKRGraph model on atom containers
    // using different constrains and search options. They give an exemple of the
    // most common queries but of course it is possible to define other type of
    // queries exploiting the constrain and option combinations
    //
    ////
    // Isomorphism search
    /**
     * Tests if g1 and g2 are isomorph.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return true if the 2 molecule are isomorph
     * @throws CDKException if the first molecule is an instance of
     * IQueryAtomContainer
     */
    public static boolean isIsomorph(IAtomContainer g1, IAtomContainer g2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        if (g1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        if (g2.getAtomCount() != g1.getAtomCount()) {
            return false;
        }
        // check single atom case
        if (g2.getAtomCount() == 1) {
            IAtom atom = g1.getAtom(0);
            IAtom atom2 = g2.getAtom(0);
            if (atom instanceof IQueryAtom) {
                IQueryAtom qAtom = (IQueryAtom) atom;
                return qAtom.matches(g2.getAtom(0));
            } else if (atom2 instanceof IQueryAtom) {
                IQueryAtom qAtom = (IQueryAtom) atom2;
                return qAtom.matches(g1.getAtom(0));
            } else {
                String atomSymbol = atom2.getSymbol();
                return g1.getAtom(0).getSymbol().equals(atomSymbol);
            }
        }
        return (getIsomorphMap(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType) != null);
    }

    /**
     * Returns the first isomorph mapping found or null.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the first isomorph mapping found projected of g1. This is a List
     * of CDKRMap objects containing Ids of matching bonds.
     * @throws CDKException
     */
    private static List<CDKRMap> getIsomorphMap(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        if (g1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        List<CDKRMap> result = null;

        List<List<CDKRMap>> rMapsList = search(g1, g2, getBitSet(g1), getBitSet(g2), false, false, shouldMatchBonds, shouldMatchRings, matchAtomType);

        if (!rMapsList.isEmpty()) {
            result = rMapsList.get(0);
        }

        return result;
    }

    /**
     * Returns the first isomorph 'atom mapping' found for g2 in g1.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the first isomorph atom mapping found projected on g1. This is a
     * List of CDKRMap objects containing Ids of matching atoms.
     * @throws CDKException if the first molecules is not an instance of
     * {@link IQueryAtomContainer}
     */
    private static List<CDKRMap> getIsomorphAtomsMap(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        if (g1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        List<CDKRMap> list = checkSingleAtomCases(g1, g2);
        if (list == null) {
            return makeAtomsMapOfBondsMap(getIsomorphMap(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType), g1, g2);
        } else if (list.isEmpty()) {
            return null;
        } else {
            return list;
        }
    }

    /**
     * Returns all the isomorph 'mappings' found between two atom containers.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the list of all the 'mappings'
     * @throws CDKException
     */
    public static List<List<CDKRMap>> getIsomorphMaps(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        return search(g1, g2, getBitSet(g1), getBitSet(g2), true, true, shouldMatchBonds, shouldMatchRings, matchAtomType);
    }

    /////
    // Subgraph search
    /**
     * Returns all the subgraph 'bondA mappings' found for g2 in g1. This is an
     * {@link List} of {@link List}s of {@link CDKRMap} objects.
     *
     * Note that if the query molecule is a single atom, then bondA mappings
     * cannot be defined. In such a case, the {@link CDKRMap} object refers
     * directly to atom - atom mappings. Thus CDKRMap.id1 is the index of the
     * target atom and CDKRMap.id2 is the index of the matching query atom (in
     * this case, it will always be 0). Note that in such a case, there is no
     * need to call
     * {@link #makeAtomsMapsOfBondsMaps(List, IAtomContainer, IAtomContainer)},
     * though if it is called, then the return value is simply the same as the
     * return value of this method.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the list of all the 'mappings' found projected of g1
     *
     * @throws CDKException
     * @see #makeAtomsMapsOfBondsMaps(List, IAtomContainer, IAtomContainer)
     */
    public static List<List<CDKRMap>> getSubgraphMaps(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        return search(g1, g2, new BitSet(), getBitSet(g2), true, true, shouldMatchBonds, shouldMatchRings, matchAtomType);
    }

    /**
     * Returns the first subgraph 'bondA mapping' found for g2 in g1.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the first subgraph bondA mapping found projected on g1. This is a
     * {@link List} of {@link CDKRMap} objects containing Ids of matching bonds.
     * @throws CDKException
     */
    public static List<CDKRMap> getSubgraphMap(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        List<CDKRMap> result = null;
        List<List<CDKRMap>> rMapsList = search(g1, g2, new BitSet(), getBitSet(g2), false, false,
                shouldMatchBonds, shouldMatchRings, matchAtomType);

        if (!rMapsList.isEmpty()) {
            result = rMapsList.get(0);
        }

        return result;
    }

    /**
     * Returns all subgraph 'atom mappings' found for g2 in g1, where g2 must be
     * a substructure of g1. If it is not a substructure, null will be returned.
     * This is an {@link List} of {@link List}s of {@link CDKRMap} objects.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 substructure to be mapped. May be an
     * {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return all subgraph atom mappings found projected on g1. This is a
     * {@link List} of {@link CDKRMap} objects containing Ids of matching atoms.
     * @throws CDKException
     */
    public static List<List<CDKRMap>> getSubgraphAtomsMaps(IAtomContainer g1,
            IAtomContainer g2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType)
            throws CDKException {
        List<CDKRMap> list = checkSingleAtomCases(g1, g2);
        if (list == null) {
            return makeAtomsMapsOfBondsMaps(
                    getSubgraphMaps(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType), g1, g2);
        } else {
            List<List<CDKRMap>> atomsMap = new ArrayList<>();
            atomsMap.add(list);
            return atomsMap;
        }
    }

    /**
     * Returns the first subgraph 'atom mapping' found for g2 in g1, where g2
     * must be a substructure of g1. If it is not a substructure, null will be
     * returned.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 substructure to be mapped. May be an
     * {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the first subgraph atom mapping found projected on g1. This is a
     * {@link List} of {@link CDKRMap} objects containing Ids of matching atoms.
     * @throws CDKException
     */
    private static List<CDKRMap> getSubgraphAtomsMap(IAtomContainer g1,
            IAtomContainer g2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType)
            throws CDKException {
        List<CDKRMap> list = checkSingleAtomCases(g1, g2);
        if (list == null) {
            return makeAtomsMapOfBondsMap(getSubgraphMap(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType), g1, g2);
        } else if (list.isEmpty()) {
            return null;
        } else {
            return list;
        }
    }

    /**
     * Tests if g2 a subgraph of g1.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return true if g2 a subgraph on g1
     * @throws CDKException
     */
    public static boolean isSubgraph(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        if (g1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        if (g2.getAtomCount() > g1.getAtomCount()) {
            return false;
        }
        // test for single atom case
        if (g2.getAtomCount() == 1) {
            IAtom atom = g2.getAtom(0);
            for (int i = 0; i < g1.getAtomCount(); i++) {
                IAtom atom2 = g1.getAtom(i);
                if (atom instanceof IQueryAtom) {
                    IQueryAtom qAtom = (IQueryAtom) atom;
                    if (qAtom.matches(atom2)) {
                        return true;
                    }
                } else if (atom2 instanceof IQueryAtom) {
                    IQueryAtom qAtom = (IQueryAtom) atom2;
                    if (qAtom.matches(atom)) {
                        return true;
                    }
                } else {
                    return atom2.getSymbol().equals(atom.getSymbol());
                }
            }
            return false;
        }
        if (!testSubgraphHeuristics(g1, g2)) {
            return false;
        }
        return (getSubgraphMap(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType) != null);
    }

    ////
    // Maximum common substructure search
    /**
     * Returns all the maximal common substructure between twp atom containers.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the list of all the maximal common substructure found projected
     * of g1 (list of GraphAtomContainer )
     * @throws CDKException
     */
    public static List<IAtomContainer> getOverlaps(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        List<List<CDKRMap>> rMapsList = search(g1, g2, new BitSet(), new BitSet(), true, false, shouldMatchBonds, shouldMatchRings, matchAtomType);

        // projection on G1
        List<IAtomContainer> graphList = projectList(rMapsList, g1, ID1);

        // reduction of set of solution (isomorphism and substructure
        // with different 'mappings'
        return getMaximum(graphList, shouldMatchBonds, shouldMatchRings, matchAtomType);
    }

    /**
     * Transforms an GraphAtomContainer into a {@link BitSet} (which's size =
     * number of bondA in the atomContainer, all the bit are set to true).
     *
     * @param ac {@link IAtomContainer} to transform
     * @return The bitSet
     */
    public static BitSet getBitSet(IAtomContainer ac) {
        BitSet bs;
        int n = ac.getBondCount();

        if (n != 0) {
            bs = new BitSet(n);
            for (int i = 0; i < n; i++) {
                bs.set(i);
            }
        } else {
            bs = new BitSet();
        }

        return bs;
    }

    //////////////////////////////////////////////////
    //          Internal methods
    /**
     * Builds the {@link CDKRGraph} ( resolution graph ), from two atomContainer
     * (description of the two molecules to compare) This is the interface point
     * between the CDK model and the generic MCSS algorithm based on the RGRaph.
     *
     * @param g1 Description of the first molecule
     * @param g2 Description of the second molecule
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return the rGraph
     * @throws CDKException
     */
    public static CDKRGraph buildRGraph(IAtomContainer g1, IAtomContainer g2,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        CDKRGraph rGraph = new CDKRGraph();
        nodeConstructor(rGraph, g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        arcConstructor(rGraph, g1, g2);
        return rGraph;
    }

    /**
     * General {@link CDKRGraph} parsing method (usually not used directly) This
     * method is the entry point for the recursive search adapted to the atom
     * container input.
     *
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @param c1 initial condition ( bonds from g1 that must be contains in the
     * solution )
     * @param c2 initial condition ( bonds from g2 that must be contains in the
     * solution )
     * @param findAllStructure if false stop at the first structure found
     * @param findAllMap if true search all the 'mappings' for one same
     * structure
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return a List of Lists of {@link CDKRMap} objects that represent the
     * search solutions
     * @throws CDKException
     */
    public static List<List<CDKRMap>> search(IAtomContainer g1, IAtomContainer g2, BitSet c1,
            BitSet c2, boolean findAllStructure, boolean findAllMap,
            boolean shouldMatchBonds, boolean shouldMatchRings,
            boolean matchAtomType) throws CDKException {
        // handle single query atom case separately
        if (g2.getAtomCount() == 1) {
            List<List<CDKRMap>> matches = new ArrayList<>();
            IAtom queryAtom = g2.getAtom(0);

            //use atom matcher from SMSD
            AtomMatcher defaultRGraphAtomMatcher;
            if (matchAtomType) {
                defaultRGraphAtomMatcher = new DefaultAtomTypeMatcher(queryAtom, shouldMatchRings);
            } else {
                defaultRGraphAtomMatcher = new DefaultAtomMatcher(queryAtom, shouldMatchRings);
            }
            for (IAtom atom : g1.atoms()) {
                if (defaultRGraphAtomMatcher.matches(atom)) {
                    List<CDKRMap> lmap = new ArrayList<>();
                    lmap.add(new CDKRMap(g1.indexOf(atom), 0));
                    matches.add(lmap);
                }
            }
            return matches;
        }

        // reset result
        List<List<CDKRMap>> rMapsList = new ArrayList<List<CDKRMap>>();

        // build the CDKRGraph corresponding to this problem
        CDKRGraph rGraph = buildRGraph(g1, g2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        // Set time data
        setIterationManager(new IterationManager((g1.getAtomCount() + g2.getAtomCount())));
        // parse the CDKRGraph with the given constrains and options
        rGraph.parse(c1, c2, findAllStructure, findAllMap);
        List<BitSet> solutionList = rGraph.getSolutions();

        // conversions of CDKRGraph's internal solutions to G1/G2 mappings
        for (BitSet set : solutionList) {
            rMapsList.add(rGraph.bitSetToRMap(set));
        }

        return rMapsList;
    }

    //////////////////////////////////////
    //    Manipulation tools
    /**
     * Projects a list of {@link CDKRMap} on a molecule.
     *
     * @param rMapList the list to project
     * @param g the molecule on which project
     * @param id the id in the {@link CDKRMap} of the molecule g
     * @return an GraphAtomContainer
     */
    private static IAtomContainer project(List<CDKRMap> rMapList, IAtomContainer g, int id) {
        IAtomContainer ac = g.getBuilder().newInstance(IAtomContainer.class);

        Map<IAtom, IAtom> table = new HashMap<>();
        IAtom a1;
        IAtom a2;
        IAtom a;
        IBond bond;
        for (CDKRMap rMap : rMapList) {
            if (id == ID1) {
                bond = g.getBond(rMap.getId1());
            } else {
                bond = g.getBond(rMap.getId2());
            }

            a = bond.getAtom(0);
            a1 = table.get(a);

            if (a1 == null) {
                try {
                    a1 = (IAtom) a.clone();
                } catch (CloneNotSupportedException e) {
                    e.printStackTrace();
                }
                ac.addAtom(a1);
                table.put(a, a1);
            }

            a = bond.getAtom(1);
            a2 = table.get(a);

            if (a2 == null) {
                try {
                    a2 = (IAtom) a.clone();
                } catch (CloneNotSupportedException e) {
                    e.printStackTrace();
                }
                ac.addAtom(a2);
                table.put(a, a2);
            }
            IBond newBond = g.getBuilder().newInstance(IBond.class, a1, a2, bond.getOrder());
            newBond.setFlag(
                    CDKConstants.ISAROMATIC,
                    bond.getFlag(CDKConstants.ISAROMATIC));
            ac.addBond(newBond);
        }
        return ac;
    }

    /**
     * Projects a list of RMapsList on a molecule.
     *
     * @param rMapsList list of RMapsList to project
     * @param g the molecule on which project
     * @param id the id in the CDKRMap of the molecule g
     * @return a list of GraphAtomContainer
     */
    private static List<IAtomContainer> projectList(List<List<CDKRMap>> rMapsList, IAtomContainer g, int id) {
        List<IAtomContainer> graphList = new ArrayList<>();

        for (List<CDKRMap> rMapList : rMapsList) {
            IAtomContainer ac = project(rMapList, g, id);
            graphList.add(ac);
        }
        return graphList;
    }

    /**
     * Removes all redundant solution.
     *
     * @param graphList the list of structure to clean
     * @return the list cleaned
     * @throws CDKException if there is a problem in obtaining subgraphs
     */
    private static List<IAtomContainer> getMaximum(List<IAtomContainer> graphList,
            boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) throws CDKException {
        List<IAtomContainer> reducedGraphList = new ArrayList<>();
        reducedGraphList.addAll(graphList);

        for (int i = 0; i < graphList.size(); i++) {
            IAtomContainer gi = graphList.get(i);

            for (int j = i + 1; j < graphList.size(); j++) {
                IAtomContainer gj = graphList.get(j);

                // Gi included in Gj or Gj included in Gi then
                // reduce the irrelevant solution
                if (isSubgraph(gj, gi, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                    reducedGraphList.remove(gi);
                } else if (isSubgraph(gi, gj, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                    reducedGraphList.remove(gj);
                }
            }
        }
        return reducedGraphList;
    }

    /**
     * Checks for single atom cases before doing subgraph/isomorphism search.
     *
     * @param g1 GraphAtomContainer to match on. Must not be an
     * {@link IQueryAtomContainer}.
     * @param g2 GraphAtomContainer as query. May be an
     * {@link IQueryAtomContainer}.
     * @return {@link List} of {@link List} of {@link CDKRMap} objects for the
     * Atoms (not Bonds!), null if no single atom case
     * @throws CDKException if the first molecule is an instance of
     * IQueryAtomContainer
     */
    static List<CDKRMap> checkSingleAtomCases(IAtomContainer g1, IAtomContainer g2) throws CDKException {
        if (g1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        if (g2.getAtomCount() == 1) {
            List<CDKRMap> arrayList = new ArrayList<>();
            IAtom atom = g2.getAtom(0);
            if (atom instanceof IQueryAtom) {
                IQueryAtom qAtom = (IQueryAtom) atom;
                for (int i = 0; i < g1.getAtomCount(); i++) {
                    if (qAtom.matches(g1.getAtom(i))) {
                        arrayList.add(new CDKRMap(i, 0));
                    }
                }
            } else {
                String atomSymbol = atom.getSymbol();
                for (int i = 0; i < g1.getAtomCount(); i++) {
                    if (g1.getAtom(i).getSymbol().equals(atomSymbol)) {
                        arrayList.add(new CDKRMap(i, 0));
                    }
                }
            }
            return arrayList;
        } else if (g1.getAtomCount() == 1) {
            List<CDKRMap> arrayList = new ArrayList<>();
            IAtom atom = g1.getAtom(0);
            for (int i = 0; i < g2.getAtomCount(); i++) {
                IAtom atom2 = g2.getAtom(i);
                if (atom2 instanceof IQueryAtom) {
                    IQueryAtom qAtom = (IQueryAtom) atom2;
                    if (qAtom.matches(atom)) {
                        arrayList.add(new CDKRMap(0, i));
                    }
                } else if (atom2.getSymbol().equals(atom.getSymbol())) {
                    arrayList.add(new CDKRMap(0, i));
                }
            }
            return arrayList;
        } else {
            return null;
        }
    }

    /**
     * This makes maps of matching atoms out of a maps of matching bonds as
     * produced by the get(Subgraph|Ismorphism)Maps methods.
     *
     * @param l The list produced by the getMap method.
     * @param g1 The first atom container. Must not be a
     * {@link IQueryAtomContainer}.
     * @param g2 The second one (first and second as in getMap). May be an
     * {@link IQueryAtomContainer}.
     * @return A List of {@link List}s of {@link CDKRMap} objects of matching
     * Atoms.
     */
    static List<List<CDKRMap>> makeAtomsMapsOfBondsMaps(List<List<CDKRMap>> l, IAtomContainer g1, IAtomContainer g2) {
        if (l == null) {
            return l;
        }
        if (g2.getAtomCount() == 1) {
            return l; // since the CDKRMap is already an atom-atom mapping
        }
        List<List<CDKRMap>> result = new ArrayList<>();
        for (List<CDKRMap> l2 : l) {
            result.add(makeAtomsMapOfBondsMap(l2, g1, g2));
        }
        return result;
    }

    /**
     * This makes a map of matching atoms out of a map of matching bonds as
     * produced by the get(Subgraph|Ismorphism)Map methods.
     *
     * @param l The list produced by the getMap method.
     * @param g1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param g2 second molecule. May be an {@link IQueryAtomContainer}.
     * @return The mapping found projected on g1. This is a {@link List} of
     * {@link CDKRMap} objects containing Ids of matching atoms.
     */
    private static List<CDKRMap> makeAtomsMapOfBondsMap(List<CDKRMap> l, IAtomContainer g1, IAtomContainer g2) {
        if (l == null) {
            return (l);
        }
        List<CDKRMap> result = new ArrayList<>();
        for (int i = 0; i < l.size(); i++) {
            IBond bond1 = g1.getBond(l.get(i).getId1());
            IBond bond2 = g2.getBond(l.get(i).getId2());
            IAtom[] atom1 = BondManipulator.getAtomArray(bond1);
            IAtom[] atom2 = BondManipulator.getAtomArray(bond2);
            for (int j = 0; j < 2; j++) {
                List<IBond> bondsConnectedToAtom1j = g1.getConnectedBondsList(atom1[j]);
                for (int k = 0; k < bondsConnectedToAtom1j.size(); k++) {
                    if (bondsConnectedToAtom1j.get(k) != bond1) {
                        IBond testBond = bondsConnectedToAtom1j.get(k);
                        for (int m = 0; m < l.size(); m++) {
                            IBond testBond2;
                            if (l.get(m).getId1() == g1.getBondNumber(testBond)) {
                                testBond2 = g2.getBond(l.get(m).getId2());
                                for (int n = 0; n < 2; n++) {
                                    List<IBond> bondsToTest = g2.getConnectedBondsList(atom2[n]);
                                    if (bondsToTest.contains(testBond2)) {
                                        CDKRMap map;
                                        if (j == n) {
                                            map = new CDKRMap(g1.indexOf(atom1[0]), g2.indexOf(atom2[0]));
                                        } else {
                                            map = new CDKRMap(g1.indexOf(atom1[1]), g2.indexOf(atom2[0]));
                                        }
                                        if (!result.contains(map)) {
                                            result.add(map);
                                        }
                                        CDKRMap map2;
                                        if (j == n) {
                                            map2 = new CDKRMap(g1.indexOf(atom1[1]), g2.indexOf(atom2[1]));
                                        } else {
                                            map2 = new CDKRMap(g1.indexOf(atom1[0]), g2.indexOf(atom2[1]));
                                        }
                                        if (!result.contains(map2)) {
                                            result.add(map2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }

    /**
     * Builds the nodes of the {@link CDKRGraph} ( resolution graph ), from two
     * atom containers (description of the two molecules to compare)
     *
     * @param gr the target CDKRGraph
     * @param ac1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param ac2 second molecule. May be an {@link IQueryAtomContainer}.
     * @throws CDKException if it takes too long to identify overlaps
     */
    private static void nodeConstructor(
            CDKRGraph gr,
            IAtomContainer ac1,
            IAtomContainer ac2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) throws CDKException {
        if (ac1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        // resets the target graph.
        gr.clear();

        // compares each bondA of G1 to each bondA of G2
        for (int i = 0; i < ac1.getBondCount(); i++) {
            for (int j = 0; j < ac2.getBondCount(); j++) {
//                // if both bonds are compatible then create an association node
//                // in the resolution graph

                if (!isMatchFeasible(ac1.getBond(i), ac2.getBond(j), shouldMatchBonds, shouldMatchRings, matchAtomType)) {
                } else {
                    gr.addNode(new CDKRNode(i, j));
                }
            }
        }
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     * @return
     */
    protected static boolean isMatchFeasible(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {

        if (bondA1 instanceof IQueryBond) {
            if (((IQueryBond) bondA1).matches(bondA2)) {
                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
                return atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0));
            }
            return false;
        } else {
            return DefaultMatcher.matches(bondA1, bondA2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        }
    }

    /**
     * Build edges of the {@link CDKRGraph}s. This method create the edge of the
     * CDKRGraph and calculates the incompatibility and neighborhood
     * relationships between CDKRGraph nodes.
     *
     * @param gr the rGraph
     * @param ac1 first molecule. Must not be an {@link IQueryAtomContainer}.
     * @param ac2 second molecule. May be an {@link IQueryAtomContainer}.
     * @throws CDKException if it takes too long to get the overlaps
     */
    private static void arcConstructor(CDKRGraph gr, IAtomContainer ac1, IAtomContainer ac2) throws CDKException {
        // each node is incompatible with himself
        for (int i = 0; i < gr.getGraph().size(); i++) {
            CDKRNode x = gr.getGraph().get(i);
            x.getForbidden().set(i);
        }

        IBond a1;
        IBond a2;
        IBond b1;
        IBond b2;

        gr.setFirstGraphSize(ac1.getBondCount());
        gr.setSecondGraphSize(ac2.getBondCount());

        for (int i = 0; i < gr.getGraph().size(); i++) {
            CDKRNode x = gr.getGraph().get(i);

            // two nodes are neighbors if their adjacency
            // relationship in are equivalent in G1 and G2
            // else they are incompatible.
            for (int j = i + 1; j < gr.getGraph().size(); j++) {
                CDKRNode y = gr.getGraph().get(j);

                a1 = ac1.getBond(gr.getGraph().get(i).getRMap().getId1());
                a2 = ac2.getBond(gr.getGraph().get(i).getRMap().getId2());
                b1 = ac1.getBond(gr.getGraph().get(j).getRMap().getId1());
                b2 = ac2.getBond(gr.getGraph().get(j).getRMap().getId2());

                if (a2 instanceof IQueryBond) {
                    if (a1.equals(b1) || a2.equals(b2)
                            || !queryAdjacencyAndOrder(a1, b1, a2, b2)) {
                        x.getForbidden().set(j);
                        y.getForbidden().set(i);
                    } else if (hasCommonAtom(a1, b1)) {
                        x.getExtension().set(j);
                        y.getExtension().set(i);
                    }
                } else if (a1.equals(b1) || a2.equals(b2)
                        || (!getCommonSymbol(a1, b1).equals(getCommonSymbol(a2, b2)))) {
                    x.getForbidden().set(j);
                    y.getForbidden().set(i);
                } else if (hasCommonAtom(a1, b1)) {
                    x.getExtension().set(j);
                    y.getExtension().set(i);
                }
            }
        }
    }

    /**
     * Determines if two bonds have at least one atom in common.
     *
     * @param a first bondA
     * @param b second bondA
     * @return the symbol of the common atom or "" if the 2 bonds have no common
     * atom
     */
    private static boolean hasCommonAtom(IBond a, IBond b) {
        return a.contains(b.getAtom(0)) || a.contains(b.getAtom(1));
    }

    /**
     * Determines if 2 bondA have 1 atom in common and returns the common
     * symbol.
     *
     * @param a first bondA
     * @param b second bondA
     * @return the symbol of the common atom or "" if the 2 bonds have no common
     * atom
     */
    private static String getCommonSymbol(IBond a, IBond b) {
        String symbol = "";

        if (a.contains(b.getAtom(0))) {
            symbol = b.getAtom(0).getSymbol();
        } else if (a.contains(b.getAtom(1))) {
            symbol = b.getAtom(1).getSymbol();
        }

        return symbol;
    }

    /**
     * Determines if 2 bondA have 1 atom in common if second is a query
     * GraphAtomContainer.
     *
     * @param a1 first bondA
     * @param b1 second bondA
     * @return the symbol of the common atom or "" if the 2 bonds have no common
     * atom
     */
    private static boolean queryAdjacency(IBond a1, IBond b1, IBond a2, IBond b2) {

        IAtom atom1 = null;
        IAtom atom2 = null;

        if (a1.contains(b1.getAtom(0))) {
            atom1 = b1.getAtom(0);
        } else if (a1.contains(b1.getAtom(1))) {
            atom1 = b1.getAtom(1);
        }

        if (a2.contains(b2.getAtom(0))) {
            atom2 = b2.getAtom(0);
        } else if (a2.contains(b2.getAtom(1))) {
            atom2 = b2.getAtom(1);
        }

        if (atom1 != null && atom2 != null) {
            // well, this looks fishy: the atom2 is not always a IQueryAtom !
            return ((IQueryAtom) atom2).matches(atom1);
        } else {
            return atom1 == null && atom2 == null;
        }

    }

    /**
     * Determines if 2 bondA have 1 atom in common if second is a query
     * GraphAtomContainer and whether the order of the atoms is correct (atoms
     * match).
     *
     * @param bond1 first bondA
     * @param bond2 second bondA
     * @param queryBond1 first query bondA
     * @param queryBond2 second query bondA
     * @return the symbol of the common atom or "" if the 2 bonds have no common
     * atom
     */
    private static boolean queryAdjacencyAndOrder(IBond bond1, IBond bond2, IBond queryBond1, IBond queryBond2) {

        IAtom centralAtom = null;
        IAtom centralQueryAtom = null;

        if (bond1.contains(bond2.getAtom(0))) {
            centralAtom = bond2.getAtom(0);
        } else if (bond1.contains(bond2.getAtom(1))) {
            centralAtom = bond2.getAtom(1);
        }

        if (queryBond1.contains(queryBond2.getAtom(0))) {
            centralQueryAtom = queryBond2.getAtom(0);
        } else if (queryBond1.contains(queryBond2.getAtom(1))) {
            centralQueryAtom = queryBond2.getAtom(1);
        }

        if (centralAtom != null && centralQueryAtom != null
                && ((IQueryAtom) centralQueryAtom).matches(centralAtom)) {
            IQueryAtom queryAtom1 = (IQueryAtom) queryBond1.getConnectedAtom(centralQueryAtom);
            IQueryAtom queryAtom2 = (IQueryAtom) queryBond2.getConnectedAtom(centralQueryAtom);
            IAtom atom1 = bond1.getConnectedAtom(centralAtom);
            IAtom atom2 = bond2.getConnectedAtom(centralAtom);
            if (queryAtom1.matches(atom1) && queryAtom2.matches(atom2)
                    || queryAtom1.matches(atom2) && queryAtom2.matches(atom1)) {
                return true;
            } else {
                return false;
            }
        } else {
            return centralAtom == null && centralQueryAtom == null;
        }

    }

    /**
     * Checks some simple heuristics for whether the subgraph query can
     * realistically be a subgraph of the supergraph. If, for example, the
     * number of nitrogen atoms in the query is larger than that of the
     * supergraph it cannot be part of it.
     *
     * @param ac1 the supergraph to be checked. Must not be an
     * {@link IQueryAtomContainer}.
     * @param ac2 the subgraph to be tested for. May be an
     * {@link IQueryAtomContainer}.
     * @return true if the subgraph ac2 has a chance to be a subgraph of ac1
     * @throws CDKException if the first molecule is an instance of
     * {@link IQueryAtomContainer}
     */
    private static boolean testSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2)
            throws CDKException {
        if (ac1 instanceof IQueryAtomContainer) {
            throw new CDKException(
                    "The first IAtomContainer must not be an IQueryAtomContainer");
        }

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;
        int ac1SCount = 0;
        int ac1OCount = 0;
        int ac1NCount = 0;
        int ac1FCount = 0;
        int ac1ClCount = 0;
        int ac1BrCount = 0;
        int ac1ICount = 0;
        int ac1CCount = 0;

        int ac2SCount = 0;
        int ac2OCount = 0;
        int ac2NCount = 0;
        int ac2FCount = 0;
        int ac2ClCount = 0;
        int ac2BrCount = 0;
        int ac2ICount = 0;
        int ac2CCount = 0;

        IBond bond;
        IAtom atom;
        for (int i = 0; i < ac1.getBondCount(); i++) {
            bond = ac1.getBond(i);
            if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                ac1AromaticBondCount++;
            } else if (bond.getOrder() == IBond.Order.SINGLE) {
                ac1SingleBondCount++;
            } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                ac1DoubleBondCount++;
            } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                ac1TripleBondCount++;
            }
        }
        for (int i = 0; i < ac2.getBondCount(); i++) {
            bond = ac2.getBond(i);
            if (bond instanceof IQueryBond) {
                continue;
            }
            if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                ac2AromaticBondCount++;
            } else if (bond.getOrder() == IBond.Order.SINGLE) {
                ac2SingleBondCount++;
            } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                ac2DoubleBondCount++;
            } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                ac2TripleBondCount++;
            }
        }

        if (ac2SingleBondCount > ac1SingleBondCount) {
            return false;
        }
        if (ac2AromaticBondCount > ac1AromaticBondCount) {
            return false;
        }
        if (ac2DoubleBondCount > ac1DoubleBondCount) {
            return false;
        }
        if (ac2TripleBondCount > ac1TripleBondCount) {
            return false;
        }

        for (int i = 0; i < ac1.getAtomCount(); i++) {
            atom = ac1.getAtom(i);
            if (atom.getSymbol().equals("S")) {
                ac1SCount++;
            } else if (atom.getSymbol().equals("N")) {
                ac1NCount++;
            } else if (atom.getSymbol().equals("O")) {
                ac1OCount++;
            } else if (atom.getSymbol().equals("F")) {
                ac1FCount++;
            } else if (atom.getSymbol().equals("Cl")) {
                ac1ClCount++;
            } else if (atom.getSymbol().equals("Br")) {
                ac1BrCount++;
            } else if (atom.getSymbol().equals("I")) {
                ac1ICount++;
            } else if (atom.getSymbol().equals("C")) {
                ac1CCount++;
            }
        }
        for (int i = 0; i < ac2.getAtomCount(); i++) {
            atom = ac2.getAtom(i);
            if (atom instanceof IQueryAtom) {
                continue;
            }
            if (atom.getSymbol().equals("S")) {
                ac2SCount++;
            } else if (atom.getSymbol().equals("N")) {
                ac2NCount++;
            } else if (atom.getSymbol().equals("O")) {
                ac2OCount++;
            } else if (atom.getSymbol().equals("F")) {
                ac2FCount++;
            } else if (atom.getSymbol().equals("Cl")) {
                ac2ClCount++;
            } else if (atom.getSymbol().equals("Br")) {
                ac2BrCount++;
            } else if (atom.getSymbol().equals("I")) {
                ac2ICount++;
            } else if (atom.getSymbol().equals("C")) {
                ac2CCount++;
            }
        }

        if (ac1SCount < ac2SCount) {
            return false;
        }
        if (ac1NCount < ac2NCount) {
            return false;
        }
        if (ac1OCount < ac2OCount) {
            return false;
        }
        if (ac1FCount < ac2FCount) {
            return false;
        }
        if (ac1ClCount < ac2ClCount) {
            return false;
        }
        if (ac1BrCount < ac2BrCount) {
            return false;
        }
        if (ac1ICount < ac2ICount) {
            return false;
        }
        return ac1CCount >= ac2CCount;

    }

    /**
     * @return the timeout
     */
    public static boolean isTimeout() {
        return timeout;
    }

    /**
     * @return the iterationManager
     */
    protected static IterationManager getIterationManager() {
        return iterationManager;
    }

    /**
     * @param aIterationManager the iterationManager to set
     */
    private static void setIterationManager(IterationManager aIterationManager) {
        iterationManager = aIterationManager;
    }
}
