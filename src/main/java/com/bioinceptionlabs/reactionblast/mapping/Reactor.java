/*
 * Reactor - consolidated reactor and helper classes.
 * Merged: AbstractReactor (inlined), Debugger, MappingHandler into Reactor
 */
package com.bioinceptionlabs.reactionblast.mapping;

import com.bioinceptionlabs.reactionblast.mapping.IMappingAlgorithm;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer.MoleculeMoleculeMapping;
import com.bioinceptionlabs.reactionblast.mapping.ReactionContainer;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.CalculationProcess;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.MappingChecks;
import com.bioinceptionlabs.reactionblast.mapping.algorithm.Holder;
import com.bioinceptionlabs.reactionblast.tools.CDKSMILES;
import com.bioinceptionlabs.reactionblast.legacy.EBIMatrix;
import com.bioinceptionlabs.reactionblast.legacy.ImageGenerator;
import com.bioinceptionlabs.reactionblast.tools.MoleculeTools.BasicDebugger;
import java.io.IOException;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.ChemicalFilters.IAtomMapping;
import org.openscience.smsd.ExtAtomContainerManipulator;
import org.openscience.smsd.BaseMapping;
import static com.bioinceptionlabs.reactionblast.tools.MoleculeTools.ExtReactionManipulatorTool.deepClone;
import static java.io.File.separator;
import static java.lang.Integer.parseInt;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.valueOf;
import static java.lang.System.getProperty;
import static java.text.NumberFormat.getInstance;
import static java.util.Arrays.sort;
import static java.util.Collections.unmodifiableList;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.MAPPED;
import static org.openscience.cdk.geometry.GeometryUtil.has2DCoordinates;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getTotalFormalCharge;
import static org.openscience.smsd.ExtAtomContainerManipulator.aromatizeMolecule;
import static org.openscience.smsd.ExtAtomContainerManipulator.cloneWithIDs;
import static org.openscience.smsd.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static org.openscience.smsd.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;


/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class Reactor extends BasicDebugger implements Serializable {

    private static final long serialVersionUID = 197816786981017L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(Reactor.class);
    private final Map<Integer, Integer> rLabelledAtoms;
    private final Map<Integer, Integer> pLabelledAtoms;
    private final Map<Integer, Integer> inputRankLabelledAtomsReactant;
    private final Map<Integer, Integer> inputRankLabelledAtomsProduct;
    private final Map<Integer, IAtomContainer> educts;
    private final Map<Integer, IAtomContainer> products;
    private final List<IBond> rBonds;
    private final List<IBond> pBonds;
    private final IReaction reactionWithSTOICHIOMETRY;
    private final boolean partialMapping;
    private final IMappingAlgorithm algorithm;
    private MoleculeMoleculeMapping reactionBlastMolMapping;
    private Integer substrateAtomCounter;
    private Integer productAtomCounter;
    private int delta;
    private boolean balanceFlag;
    private IReaction reactionWithUniqueSTOICHIOMETRY;
    private final SmilesGenerator smiles;

    //~--- constructors -------------------------------------------------------
    /**
     *
     * @param reaction Reaction
     * @param balanceReaction eg. balance hydrogens in the reaction if its not
     * balanced
     * @param partialMapping (without Hydrogens is set True, its faster)
     * @throws Exception
     */
    Reactor(IReaction reaction,
            boolean partialMapping,
            IMappingAlgorithm algorithm)
            throws Exception {
        if (partialMapping) {
            //else CDKToBeam throws an error "Aromatic bond connects non-aromatic atomic atoms"
            smiles = new SmilesGenerator(
                    SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo);
        } else {
            smiles = new SmilesGenerator(
                    SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo);
        }

        this.partialMapping = partialMapping;
        this.algorithm = algorithm;
        this.reactionWithSTOICHIOMETRY = reaction.getBuilder().newInstance(IReaction.class);
        this.reactionWithUniqueSTOICHIOMETRY = reaction.getBuilder().newInstance(IReaction.class);
        this.balanceFlag = true;

        this.inputRankLabelledAtomsReactant = new HashMap<>();
        this.inputRankLabelledAtomsProduct = new HashMap<>();
        this.rLabelledAtoms = new HashMap<>();
        this.pLabelledAtoms = new HashMap<>();
        this.rBonds = new ArrayList<>();
        this.pBonds = new ArrayList<>();

        this.educts = new TreeMap<>();
        this.products = new TreeMap<>();

        this.substrateAtomCounter = 1;
        this.productAtomCounter = 1;
        LOGGER.debug("|++++++++++++++++++++++++++++|");
        LOGGER.debug("|i. Reactor Initialized");
        LOGGER.debug("|++++++++++++++++++++++++++++|");
        printReaction(reaction);
        LOGGER.debug("|ii. Create Mapping Objects");
        copyReferenceReaction(reaction);
        MappingHandler.cleanMapping(reactionWithSTOICHIOMETRY);
        expandReaction();
        checkReactionBalance();
        LOGGER.debug("|iii. Compute atom-atom Mappings");
        calculateAtomAtomMapping();
        printReaction(reactionWithUniqueSTOICHIOMETRY);
        LOGGER.debug("|iv. Done|");
    }

    @Override
    public String toString() {

        String createReactionSMILES = "";
        try {
            createReactionSMILES = smiles.create(reactionWithUniqueSTOICHIOMETRY);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return "Reactor{" + "partialMapping=" + partialMapping + ", algorithm=" + algorithm
                + ", mapping=" + createReactionSMILES + '}';
    }

    /**
     * Copy reference reaction molecules into the stoichiometry reaction.
     * Uses direct clone + perceive instead of expensive SMILES round-trip
     * (serialize → parse → perceive was ~15% of total mapping time).
     */
    private void copyReferenceReaction(IReaction referenceReaction) throws CDKException, IOException, Exception {
        try {
            for (int i = 0; i < referenceReaction.getReactantCount(); i++) {
                IAtomContainer refMol = referenceReaction.getReactants().getAtomContainer(i);
                IAtomContainer cloneMolecule = cloneWithIDs(refMol);
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cloneMolecule);
                cloneMolecule = prepareMol(cloneMolecule);
                cloneMolecule.setID(refMol.getID());
                Double st = referenceReaction.getReactantCoefficient(refMol);
                aromatizeMolecule(cloneMolecule);
                reactionWithSTOICHIOMETRY.addReactant(cloneMolecule, st);
            }
        } catch (CloneNotSupportedException | CDKException e) {
            LOGGER.error(SEVERE, null, e);
        }
        try {
            for (int i = 0; i < referenceReaction.getProductCount(); i++) {
                IAtomContainer refMol = referenceReaction.getProducts().getAtomContainer(i);
                IAtomContainer cloneMolecule = cloneWithIDs(refMol);
                ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cloneMolecule);
                cloneMolecule = prepareMol(cloneMolecule);
                cloneMolecule.setID(refMol.getID());
                Double st = referenceReaction.getProductCoefficient(refMol);
                aromatizeMolecule(cloneMolecule);
                reactionWithSTOICHIOMETRY.addProduct(cloneMolecule, st);
            }
            reactionWithSTOICHIOMETRY.setID(referenceReaction.getID());
            reactionWithSTOICHIOMETRY.setDirection(referenceReaction.getDirection());
        } catch (CloneNotSupportedException | CDKException e) {
            LOGGER.error(SEVERE, "Error in Reactor class", e.getMessage());
        }
    }

    private void expandReaction() throws CloneNotSupportedException {

        for (int i = 0; i < reactionWithSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer _react = reactionWithSTOICHIOMETRY.getReactants().getAtomContainer(i);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getReactantCoefficient(_react);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer _reactDup = cloneWithIDs(_react);
                _reactDup.setID(_react.getID());
                _reactDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addReactant(_reactDup, 1.0);
            }
        }

        for (int j = 0; j < reactionWithSTOICHIOMETRY.getProductCount(); j++) {

            IAtomContainer _prod = reactionWithSTOICHIOMETRY.getProducts().getAtomContainer(j);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getProductCoefficient(_prod);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer prodDup = cloneWithIDs(_prod);
                prodDup.setID(_prod.getID());
                prodDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addProduct(prodDup, 1.0);
            }

        }

        reactionWithUniqueSTOICHIOMETRY.setID(
                reactionWithSTOICHIOMETRY.getID() == null
                ? "MappedReaction (ecBLAST)"
                : reactionWithSTOICHIOMETRY.getID());
        reactionWithUniqueSTOICHIOMETRY.setDirection(reactionWithSTOICHIOMETRY.getDirection() == null
                ? BIDIRECTIONAL
                : reactionWithSTOICHIOMETRY.getDirection());

        LabelAtoms();
        BondCollection();
    }

    private void LabelAtoms() {
        int new_atom_rank_index_reactant = 1;
        int new_atom_rank_index_product = 1;
        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (substrateAtomCounter).toString();
                substrateAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
                rLabelledAtoms.put(atom.hashCode(), i);
                if (atom.getProperty("OLD_RANK") != null) {
                    inputRankLabelledAtomsReactant.put((int) atom.getProperty("OLD_RANK"), (new_atom_rank_index_reactant++));
                }
            }

            educts.put(i, container);
        }

        for (int j = 0; j < reactionWithUniqueSTOICHIOMETRY.getProductCount(); j++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(j);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (productAtomCounter).toString();
                productAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
                pLabelledAtoms.put(atom.hashCode(), j);
                if (atom.getProperty("OLD_RANK") != null) {
                    inputRankLabelledAtomsProduct.put((int) atom.getProperty("OLD_RANK"), (new_atom_rank_index_product++));
                }
            }

            products.put(j, container);
        }

    }

    private void BondCollection() {

        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer mol = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
            for (int j = 0; j < mol.getBondCount(); j++) {
                IBond bond = mol.getBond(j);
                if (bond != null) {
                    rBonds.add(bond);
                }
            }
        }

        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getProductCount(); i++) {
            IAtomContainer mol = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(i);
            for (int j = 0; j < mol.getBondCount(); j++) {
                IBond bond = mol.getBond(j);
                if (bond != null) {
                    pBonds.add(bond);
                }
            }
        }

    }

    private void checkReactionBalance() throws IOException {
        IAtomContainerSet reactantSet = getExpandedReactants();
        IAtomContainerSet productSet = getExpandedProducts();
        HashMap<String, Integer> AtomMap = new HashMap<>();
        for (int i = 0; i < reactantSet.getAtomContainerCount(); i++) {
            IAtomContainer rMol = reactantSet.getAtomContainer(i);
            Iterator<IAtom> rAtomIterator = rMol.atoms().iterator();
            while (rAtomIterator.hasNext()) {
                IAtom rAtom = rAtomIterator.next();
                if (!rAtom.getSymbol().equals("H")) {
                    if (AtomMap.containsKey(rAtom.getSymbol())) {
                        int count = AtomMap.get(rAtom.getSymbol()) + 1;
                        AtomMap.put(rAtom.getSymbol(), count);
                    } else {
                        AtomMap.put(rAtom.getSymbol(), 1);
                    }
                }
            }
        }

        for (int j = 0; j < productSet.getAtomContainerCount(); j++) {
            IAtomContainer pMol = productSet.getAtomContainer(j);
            Iterator<IAtom> pAtomIterator = pMol.atoms().iterator();
            while (pAtomIterator.hasNext()) {
                IAtom pAtom = pAtomIterator.next();
                if (!pAtom.getSymbol().equals("H")) {
                    if (AtomMap.containsKey(pAtom.getSymbol())) {
                        int count = AtomMap.get(pAtom.getSymbol()) - 1;
                        AtomMap.put(pAtom.getSymbol(), count);
                    } else if (!AtomMap.containsKey(pAtom.getSymbol())) {
                        AtomMap.put(pAtom.getSymbol(), 1);
                        this.balanceFlag = false;
                        break;
                    }
                }
            }
        }

        for (Map.Entry<String, Integer> I : AtomMap.entrySet()) {
            if (I.getValue() != 0) {
                this.balanceFlag = false;
                break;
            }
        }
    }

    private void calculateAtomAtomMapping() throws IOException, Exception {

        try {
            IReaction reactionCopy = copyReaction(reactionWithUniqueSTOICHIOMETRY, partialMapping);
            CalculationProcess calP
                    = new CalculationProcess(partialMapping, reactionCopy, getAlgorithm());
            delta = calP.getDelta();
            IReaction mappedReaction = calP.getMappedReaction();
            reactionWithUniqueSTOICHIOMETRY = getMapping(mappedReaction);
            setReactionBlastMolMapping(calP.getReactionBlastMolMapping());
        } catch (Exception ex) {
            LOGGER.error(SEVERE, "Error in Reactor class", ex);
        }
    }

    private IReaction getMapping(IReaction coreMappedReaction) throws IOException, CDKException, CloneNotSupportedException {

        IReaction mappedReaction = deepClone(reactionWithUniqueSTOICHIOMETRY);
        MappingHandler.cleanMapping(mappedReaction);

        /*
        * This section set the mappingMap ID for the mapped atoms
         */
        int counter = 1;

        counter = MappingHandler.setMappingFlags(mappedReaction, reactionWithUniqueSTOICHIOMETRY, coreMappedReaction, counter);

        /*
        * This section set the mappingMap ID for the unmapped atoms
        *
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && "-1".equalsIgnoreCase(atom.getID())) {
                    String atomLabel = Integer.toString(counter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                }
                counter += 1;
            }
        }

        for (int pMol = 0; pMol < mappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = mappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = mappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && "-1".equalsIgnoreCase(atom.getID())) {
                    String atomLabel = Integer.toString(counter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                    counter += 1;
                }
            }
        }

        /*
        * This section will mark map common H atoms. example H-R + H <=> R-H + H Here R-H will be mapped to the R-H.
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && !"-1".equalsIgnoreCase(atom.getID())) {
                    List<IAtom> eductConnAtoms = eMolecule.getConnectedAtomsList(atom);
                    List<IAtom> productHAtoms = markHAroundCoreAtoms(atom.getID(), mappedReaction.getProducts());
                    for (IAtom eAtomH : eductConnAtoms) {
                        //Collect ummmarked H and map common ones
                        if ("-1".equalsIgnoreCase(eAtomH.getID()) && eAtomH.getSymbol().equalsIgnoreCase("H")) {
                            if (!productHAtoms.isEmpty()) {
                                String atomLabel = Integer.toString(counter);
                                eAtomH.setID(atomLabel);
                                eAtomH.setFlag(MAPPED, true);
                                IAtom pAtomH = productHAtoms.iterator().next();
                                pAtomH.setID(atomLabel);
                                pAtomH.setFlag(MAPPED, true);
                                productHAtoms.remove(pAtomH);
                                counter += 1;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        }

        /*
        *
        * This section will mark single unmapped H atoms on both sides. example R-H + H <=> H + R-H Here H and H will
        * be marked and matched
        *
         */
        List<IAtom> unMappedSingleHAtEduct = collectUnMappedSingleHAtoms(mappedReaction.getReactants());
        List<IAtom> unMappedSingleHAtProduct = collectUnMappedSingleHAtoms(mappedReaction.getProducts());

        /*
        * Mark single unmapped Hs on both sides
         */
        for (IAtom eAtomH : unMappedSingleHAtEduct) {
            if (!unMappedSingleHAtProduct.isEmpty()) {
                String atomLabel = Integer.toString(counter);
                eAtomH.setID(atomLabel);
                eAtomH.setFlag(MAPPED, true);
                IAtom pAtomH = unMappedSingleHAtProduct.iterator().next();
                pAtomH.setID(atomLabel);
                pAtomH.setFlag(MAPPED, true);
                unMappedSingleHAtProduct.remove(pAtomH);
                counter += 1;
            } else {
                break;
            }
        }

        /*
        *
        * This section will mark unmapped H atoms on both sides. example R-H + H <=> H-H + R Here H and H will be
        * marked and matched
        *
         */
        List<IAtom> unMappedHAtEduct = collectUnMappedHAtoms(mappedReaction.getReactants());
        List<IAtom> unMappedHAtProduct = collectUnMappedHAtoms(mappedReaction.getProducts());

        /*
        * Mark single unmapped Hs on both sides
         */
        for (IAtom eAtomH : unMappedHAtEduct) {
            if (!unMappedHAtProduct.isEmpty()) {
                String atomLabel = Integer.toString(counter);
                eAtomH.setID(atomLabel);
                eAtomH.setFlag(MAPPED, true);
                IAtom pAtomH = unMappedHAtProduct.iterator().next();
                pAtomH.setID(atomLabel);
                pAtomH.setFlag(MAPPED, true);
                unMappedHAtProduct.remove(pAtomH);
                counter += 1;
            } else {
                break;
            }
        }

        /*
         * Mark unmapped H atoms i.e. protonations
         */
        counter = markUnMappedHAtoms(mappedReaction, counter);

        /*
        * Create atom-atom mappingMap objects to be stored in a map
        *
         */
        Map<IAtom, IAtom> mappings = new HashMap<>();
        for (IAtomContainer ac1 : mappedReaction.getReactants().atomContainers()) {
            for (IAtom atom1 : ac1.atoms()) {
                IAtom atom2 = getContainerAtomByID(mappedReaction.getProducts(), atom1.getID());
                if (atom2 != null) {
                    mappings.put(atom1, atom2);
                }
            }
        }

        /*
         * Store atom-atom mappingMap objects in the reaction
         *
         */
        mappings.keySet().stream().filter((key)
                -> (key != null && mappings.get(key) != null)).map((key)
                -> mappedReaction.getBuilder().newInstance(IMapping.class, key, mappings.get(key))).forEachOrdered((mappingObject) -> {
            mappedReaction.addMapping(mappingObject);
        });

        /*
        * Canonical labelling of each molecule is done and mappingMap number corresponds to the lables
        *
         */
        counter = setCanonicalMappingLabels(mappedReaction);
        return mappedReaction;
    }

//~--- get methods --------------------------------------------------------
    /**
     *
     * @return reactantSet expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
        public IAtomContainerSet getExpandedReactants() throws IOException {
        return reactionWithUniqueSTOICHIOMETRY.getReactants();
    }

    /**
     *
     * @return productSet expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
        public IAtomContainerSet getExpandedProducts() throws IOException {
        return reactionWithUniqueSTOICHIOMETRY.getProducts();
    }

    /**
     *
     * @return IReaction object with unique atom labeling
     * @throws Exception
     */
        public IReaction getReactionWithAtomAtomMapping() throws Exception {
        return reactionWithUniqueSTOICHIOMETRY;
    }

    /**
     *
     * @param i Index Ith position
     * @return Stoichiometry weight of the reactant molecule at ith Position
     */
        public Double getExpandedReactantStoichiometry(
            int i) {
        IAtomContainer Mol = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
        return reactionWithUniqueSTOICHIOMETRY.getReactantCoefficient(Mol);

    }

    /**
     *
     * @param i Index at I th position
     * @return Stoichiometry weight of the product molecule at i th Position
     *
     */
        public Double getExpandedProductStoichiometry(int i) {
        IAtomContainer Mol = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(i);
        return reactionWithUniqueSTOICHIOMETRY.getProductCoefficient(Mol);
    }

    /**
     *
     * @return true if its a balanced reaction else false
     * @throws IOException
     *
     */
        public boolean getReactionBalanceFlag() throws IOException {

        boolean flag = true;
        if (!Objects.equals(this.getLabledReactantAtomsCount(), this.getLabledProductAtomsCount())) {
            flag = false;
        }
        if (!getReactionBalanceFlagWithoutHydrogen()) {
            flag = false;
        }
        return flag;
    }

    /**
     *
     * @return @throws IOException
     */
        public boolean getReactionBalanceFlagWithChargeBalance() throws IOException {
        boolean flag = true;
        if (!Objects.equals(this.getLabledReactantAtomsCount(), this.getLabledProductAtomsCount())) {
            flag = false;
        }
        if (getTotalFormalCharge(this.getExpandedReactants()) != getTotalFormalCharge(this.getExpandedProducts())) {
            flag = false;
        }
        if (!getReactionBalanceFlagWithoutHydrogen()) {
            flag = false;
        }
        return flag;
    }

    /**
     *
     * @return true if its a balanced reaction else false Note: This does not
     * consider whether Hydrogens are balanced or not
     *
     *
     */
        public boolean getReactionBalanceFlagWithoutHydrogen() {
        return this.balanceFlag;
    }

    /**
     *
     * @return this will return IAtom Vector of total Reactant atom count with
     * unique labeling
     *
     *
     */
    private List<IAtom> getLabledReactantAtoms() {
        List<IAtom> reactantAtoms = new ArrayList<>();
        IAtomContainerSet MSet = reactionWithUniqueSTOICHIOMETRY.getReactants();
        for (int j = 0; j
                < MSet.getAtomContainerCount(); j++) {
            IAtomContainer M = MSet.getAtomContainer(j);
            for (int k = 0; k
                    < M.getAtomCount(); k++) {
                reactantAtoms.add(M.getAtom(k));
            }
        }
        return unmodifiableList(reactantAtoms);
    }

    /**
     *
     * @return this will return IAtom Vector of total Product atom count with
     * unique labelling
     *
     */
    private List<IAtom> getLabledProductAtoms() {
        List<IAtom> productAtoms = new ArrayList<>();
        IAtomContainerSet MSet = reactionWithUniqueSTOICHIOMETRY.getProducts();
        for (int j = 0; j
                < MSet.getAtomContainerCount(); j++) {
            IAtomContainer M = MSet.getAtomContainer(j);
            for (int k = 0; k
                    < M.getAtomCount(); k++) {
                productAtoms.add(M.getAtom(k));
            }
        }
        return unmodifiableList(productAtoms);
    }

    /**
     *
     * @return this will return the total reactant
     *
     * atom count with unique labeling
     *
     *
     */
    private Integer getLabledReactantAtomsCount() {
        return getLabledReactantAtoms().size();
    }

    /**
     *
     * @return this will return the total product
     *
     * atom count with unique labeling
     *
     */
    private Integer getLabledProductAtomsCount() {
        return getLabledProductAtoms().size();

    }

    /**
     *
     * @return bonds of reactantSet
     */
        public List<IBond> getEductBonds() {
        return unmodifiableList(rBonds);
    }

    /**
     * @return bonds of productSet
     */
        public List<IBond> getProductBonds() {
        return unmodifiableList(pBonds);
    }

    /**
     *
     * @return
     */
        public int getMappingCount() {
        return reactionWithUniqueSTOICHIOMETRY.getMappingCount();
    }

    private IReaction copyReaction(IReaction orignalReaction, boolean removeHydrogen) throws Exception {
        IReaction copiedReaction = reactionWithUniqueSTOICHIOMETRY.getBuilder().newInstance(IReaction.class);

        for (int i = 0; i < orignalReaction.getReactantCount(); i++) {
            IAtomContainer mol = orignalReaction.getReactants().getAtomContainer(i);
            Double st = orignalReaction.getReactantCoefficient(mol);
            IAtomContainer newMol = cloneWithIDs(mol);
            for (int index = 0; index < mol.getAtomCount(); index++) {
                mol.getAtom(index).setProperty("index", index);
                IAtom atom = newMol.getAtom(index);
                atom.setProperty("index", index);
            }

            percieveAtomTypesAndConfigureAtoms(newMol);
            if (removeHydrogen) {
                newMol = removeHydrogensExceptSingleAndPreserveAtomID(newMol);
            }
            copiedReaction.addReactant(newMol, st);
        }
        for (int i = 0; i < orignalReaction.getProductCount(); i++) {
            IAtomContainer mol = orignalReaction.getProducts().getAtomContainer(i);
            Double st = orignalReaction.getProductCoefficient(mol);
            IAtomContainer newMol = cloneWithIDs(mol);
            for (int index = 0; index < mol.getAtomCount(); index++) {
                mol.getAtom(index).setProperty("index", index);
                IAtom atom = newMol.getAtom(index);
                atom.setProperty("index", index);
            }

            percieveAtomTypesAndConfigureAtoms(newMol);
            if (removeHydrogen) {
                newMol = removeHydrogensExceptSingleAndPreserveAtomID(newMol);
            }
            copiedReaction.addProduct(newMol, st);
        }
        copiedReaction.setFlags(orignalReaction.getFlags());
        copiedReaction.setID(orignalReaction.getID());
        copiedReaction.setDirection(orignalReaction.getDirection());
        copiedReaction.notifyChanged();
        return copiedReaction;
    }

    /**
     *
     * @param id
     * @param molSet
     * @return
     */
    private List<IAtom> markHAroundCoreAtoms(String id, IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int pMol = 0; pMol < molSet.getAtomContainerCount(); pMol++) {
            IAtomContainer pMolecule = molSet.getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = molSet.getAtomContainer(pMol).getAtom(pAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && !"-1".equalsIgnoreCase(atom.getID())) {
                    if (id != null && id.equalsIgnoreCase(atom.getID())) {
                        List<IAtom> conAtoms = pMolecule.getConnectedAtomsList(atom);
                        conAtoms.stream().filter((atomH) -> ("-1".equalsIgnoreCase(atomH.getID()) && atomH.getSymbol().equalsIgnoreCase("H"))).forEach((atomH) -> {
                            list.add(atomH);
                        });
                    }
                }
            }
        }
        return list;
    }

    /**
     * @param molSet
     * @return
     */
    private List<IAtom> collectUnMappedSingleHAtoms(IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int index = 0; index < molSet.getAtomContainerCount(); index++) {
            IAtomContainer mol = molSet.getAtomContainer(index);
            if (mol.getAtomCount() == 1) {
                for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                    IAtom atom = molSet.getAtomContainer(index).getAtom(atomIndex);
                    if (atom.getSymbol().equalsIgnoreCase("H")
                            && !atom.getFlag(MAPPED)
                            && "-1".equalsIgnoreCase(atom.getID())) {
                        list.add(atom);
                    }
                }
            }
        }
        return list;
    }

    /**
     * @param molSet
     * @return
     */
    private List<IAtom> collectUnMappedHAtoms(IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int index = 0; index < molSet.getAtomContainerCount(); index++) {
            IAtomContainer mol = molSet.getAtomContainer(index);
            for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                IAtom atom = molSet.getAtomContainer(index).getAtom(atomIndex);
                if (atom.getSymbol().equalsIgnoreCase("H")
                        && !atom.getFlag(MAPPED)
                        && "-1".equalsIgnoreCase(atom.getID())) {
                    list.add(atom);
                }
            }
        }
        return list;
    }

    /**
     *
     * @param mappedReaction
     * @param counter
     * @return updated Counter
     */
    private int markUnMappedHAtoms(IReaction mappedReaction, int counter) {

        int localCounter = counter;

        /*
         * Mark unmapped H atoms
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (atom.getSymbol().equalsIgnoreCase("H") && !atom.getFlag(MAPPED)
                        && "-1".equalsIgnoreCase(atom.getID())) {
                    String atomLabel = Integer.toString(localCounter);
                    atom.setFlag(MAPPED, false);
                    atom.setID(atomLabel);
                    localCounter += 1;
                }
            }
        }

        for (int pMol = 0; pMol < mappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = mappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = mappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                if (atom.getSymbol().equalsIgnoreCase("H") && !atom.getFlag(MAPPED)
                        && "-1".equalsIgnoreCase(atom.getID())) {
                    String atomLabel = Integer.toString(localCounter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                    localCounter += 1;
                }
            }
        }
        return localCounter;
    }

    /**
     * @return the delta
     */
    public int getDelta() {
        return delta;
    }

    /**
     * @return the reactionBlastMolMapping
     */
    public MoleculeMoleculeMapping getReactionBlastMolMapping() {
        return reactionBlastMolMapping;
    }

    /**
     * @param reactionBlastMolMapping the reactionBlastMolMapping to set
     */
    private void setReactionBlastMolMapping(MoleculeMoleculeMapping reactionBlastMolMapping) {
        this.reactionBlastMolMapping = reactionBlastMolMapping;
    }

    private IAtom getContainerAtomByID(IAtomContainerSet products, String mappingID) {
        if (mappingID == null) {
            return null;
        }
        for (IAtomContainer ac : products.atomContainers()) {
            for (IAtom atom : ac.atoms()) {
                if (mappingID.equals(atom.getID())) {
                    return atom;
                }
            }
        }
        return null;
    }

    private int setCanonicalMappingLabels(IReaction mappedReaction) throws CDKException {
        IAtomContainerSet rMolSet = mappedReaction.getReactants();
        IAtomContainerSet pMolSet = mappedReaction.getProducts();

        Map<IAtom, IAtom> mappingMap = new HashMap<>();

        for (IMapping aaMapping : mappedReaction.mappings()) {
            aaMapping.getChemObject(0).removeProperty(ATOM_ATOM_MAPPING);
            aaMapping.getChemObject(1).removeProperty(ATOM_ATOM_MAPPING);
            mappingMap.put((IAtom) aaMapping.getChemObject(0), (IAtom) aaMapping.getChemObject(1));
        }

        /*
         * Re-arrange the molecule index for mappings
         */
        for (IAtomContainer mol : rMolSet.atomContainers()) {
            List<Integer> atom_index = new ArrayList<>();

            for (IAtom a : mol.atoms()) {
                if (!a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            for (IAtom a : mol.atoms()) {
                if (a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            int[] array = new int[atom_index.size()];

            int index = 0;
            for (int c : atom_index) {
                array[index] = c;
                index++;
            }
            permuteWithoutClone(array, mol);
        }

        for (IAtomContainer mol : pMolSet.atomContainers()) {

            List<Integer> atom_index = new ArrayList<>();
            for (IAtom a : mol.atoms()) {
                if (!a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            for (IAtom a : mol.atoms()) {
                if (a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            int[] array = new int[atom_index.size()];

            int index = 0;
            for (int c : atom_index) {
                array[index] = c;
                index++;
            }
            permuteWithoutClone(array, mol);
        }

        int counter = 1;
        for (IAtom qAtom : collectMappedAtomsByOriginalRank(rMolSet, mappingMap, false)) {
            assignMappedLabel(qAtom, mappingMap.get(qAtom), counter++);
        }

        for (IAtom qAtom : collectMappedAtomsByOriginalRank(rMolSet, mappingMap, true)) {
            assignMappedLabel(qAtom, mappingMap.get(qAtom), counter++);
        }

        for (IAtom qAtom : collectUnmappedAtomsByOriginalRank(rMolSet, mappingMap.keySet(), false)) {
            assignUnmappedLabel(qAtom, counter++);
        }

        for (IAtom qAtom : collectUnmappedAtomsByOriginalRank(rMolSet, mappingMap.keySet(), true)) {
            assignUnmappedLabel(qAtom, counter++);
        }

        for (IAtom tAtom : collectUnmappedAtomsByOriginalRank(pMolSet, new HashSet<>(mappingMap.values()), false)) {
            assignUnmappedLabel(tAtom, counter++);
        }

        for (IAtom tAtom : collectUnmappedAtomsByOriginalRank(pMolSet, new HashSet<>(mappingMap.values()), true)) {
            assignUnmappedLabel(tAtom, counter++);
        }
        /*
        Finally permute molecules based on the atom mapping rank
         */
        for (IAtomContainer mol : pMolSet.atomContainers()) {
            TreeMap<Integer, Integer> mapping_rank = new TreeMap<>();
            for (IAtom a : mol.atoms()) {
                mapping_rank.put((Integer) a.getProperty(ATOM_ATOM_MAPPING), mol.indexOf(a));
            }
            int[] mappingIndexPermutation = new int[mapping_rank.size()];
            int index = 0;
            for (int i : mapping_rank.values()) {
                mappingIndexPermutation[index] = i;
                index++;
            }
            permuteWithoutClone(mappingIndexPermutation, mol);
        }

        mappingMap.clear();
        return counter;
    }

    private void assignMappedLabel(IAtom reactantAtom, IAtom productAtom, int counter) {
        String id = valueOf(counter);
        reactantAtom.setID(id);
        reactantAtom.setProperty(ATOM_ATOM_MAPPING, counter);
        reactantAtom.setMapIdx(counter);

        if (productAtom != null) {
            productAtom.setID(id);
            productAtom.setProperty(ATOM_ATOM_MAPPING, counter);
            productAtom.setMapIdx(counter);
        }
    }

    private void assignUnmappedLabel(IAtom atom, int counter) {
        String id = valueOf(counter);
        atom.setID(id);
        atom.setProperty(ATOM_ATOM_MAPPING, counter);
        atom.setMapIdx(counter);
    }

    private List<IAtom> collectMappedAtomsByOriginalRank(IAtomContainerSet molSet,
            Map<IAtom, IAtom> mappingMap, boolean hydrogens) {
        List<IAtom> atoms = new ArrayList<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                if (mappingMap.containsKey(atom)
                        && atom.getSymbol().equalsIgnoreCase("H") == hydrogens) {
                    atoms.add(atom);
                }
            }
        }
        sortByOriginalRank(atoms);
        return atoms;
    }

    private List<IAtom> collectUnmappedAtomsByOriginalRank(IAtomContainerSet molSet,
            Set<IAtom> mappedAtoms, boolean hydrogens) {
        List<IAtom> atoms = new ArrayList<>();
        for (IAtomContainer mol : molSet.atomContainers()) {
            for (IAtom atom : mol.atoms()) {
                if (!mappedAtoms.contains(atom)
                        && atom.getSymbol().equalsIgnoreCase("H") == hydrogens) {
                    atoms.add(atom);
                }
            }
        }
        sortByOriginalRank(atoms);
        return atoms;
    }

    private void sortByOriginalRank(List<IAtom> atoms) {
        atoms.sort((left, right) -> {
            int rankComparison = Integer.compare(getOriginalRank(left), getOriginalRank(right));
            if (rankComparison != 0) {
                return rankComparison;
            }

            int labelComparison = Integer.compare(getStableAtomPosition(left), getStableAtomPosition(right));
            if (labelComparison != 0) {
                return labelComparison;
            }
            return left.getSymbol().compareTo(right.getSymbol());
        });
    }

    private int getOriginalRank(IAtom atom) {
        Object oldRank = atom.getProperty("OLD_RANK");
        if (oldRank instanceof Integer) {
            return (Integer) oldRank;
        }
        if (oldRank != null) {
            try {
                return parseInt(oldRank.toString());
            } catch (NumberFormatException ignore) {
                // Fall through to the stable fallback below.
            }
        }
        return getStableAtomPosition(atom);
    }

    private int getStableAtomPosition(IAtom atom) {
        Object label = atom.getProperty("label");
        if (label instanceof Integer) {
            return (Integer) label;
        }
        if (label != null) {
            try {
                return parseInt(label.toString());
            } catch (NumberFormatException ignore) {
                // Fall through to the generic fallback below.
            }
        }
        Object index = atom.getProperty("index");
        if (index instanceof Integer) {
            return (Integer) index;
        }
        if (index != null) {
            try {
                return parseInt(index.toString());
            } catch (NumberFormatException ignore) {
                return Integer.MAX_VALUE;
            }
        }
        return Integer.MAX_VALUE;
    }

    /**
     * @return the algorithm
     */
    public IMappingAlgorithm getAlgorithm() {
        return algorithm;
    }

    private IAtomContainer prepareMol(IAtomContainer cloneMolecule)
            throws CloneNotSupportedException, CDKException {

        LOGGER.debug("Original");
        printAtoms(cloneMolecule);
        /*
        Use the Canonical labelling from the SMILES
        IMP: Suggested by John May
         */
        int[] p = new int[cloneMolecule.getAtomCount()];

        try {
            //this helps to avoid concurrent modification error, reason unknown
            String sm = unique().create(cloneMolecule, p);
            LOGGER.debug("smiles " + sm);
        } catch (CDKException e) {
            LOGGER.error(SEVERE, null, e);
        }
        permuteWithoutClone(p, cloneMolecule);

        LOGGER.debug("mol after: ");
        printAtoms(cloneMolecule);

        /*
        Generate 2D Diagram without cloning
         */
        if (!has2DCoordinates(cloneMolecule)) {
            try {
                /*
                Clone it else it will loose mol ID
                 */
                StructureDiagramGenerator sdg = new StructureDiagramGenerator();
                sdg.setMolecule(cloneMolecule, false);
                sdg.generateCoordinates();
            } catch (CDKException e) {
                LOGGER.error(SEVERE, "Error in 2D Generation ", e.getMessage());
            }
        }

        /*
        Set the IDs to -1 very IMP
         */
        for (IAtom atom : cloneMolecule.atoms()) {
            atom.setID("-1");
        }

        LOGGER.debug("Processed");
        printAtoms(cloneMolecule);

        return cloneMolecule;
    }

    /*
     * This is a very imp code modified by John May
     * The idea is to canonicalise the atoms and bonds
     */
    private void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
        int n = atomContainer.getAtomCount();
        LOGGER.debug("permuting " + java.util.Arrays.toString(p));
        IAtom[] permutedAtoms = new IAtom[n];

        for (int i = 0; i < n; i++) {
            IAtom atom = atomContainer.getAtom(i);
            permutedAtoms[p[i]] = atom;
            atom.setProperty("label", p[i]);
        }
        atomContainer.setAtoms(permutedAtoms);

        IBond[] bonds = java.util.Arrays.stream(getBondArray(atomContainer))
                .filter(Objects::nonNull)
                .toArray(IBond[]::new);
        sort(bonds, (IBond o1, IBond o2) -> {
            int u = o1.getAtom(0).getProperty("label");
            int v = o1.getAtom(1).getProperty("label");
            int x = o2.getAtom(0).getProperty("label");
            int y = o2.getAtom(1).getProperty("label");
            int min1 = min(u, v);
            int min2 = min(x, y);
            int max1 = max(u, v);
            int max2 = max(x, y);

            int minCmp = Integer.compare(min1, min2);
            if (minCmp != 0) {
                return minCmp;
            }
            int maxCmp = Integer.compare(max1, max2);
            if (maxCmp != 0) {
                return maxCmp;
            }
            LOGGER.debug("pokemon!");
            throw new InternalError();
        });
        atomContainer.setBonds(bonds);
    }

    /**
     * Old Atom Rank in the reactant mapped to new Rank
     *
     * @return the inputRankLabelledAtomsReactant
     */
    public Map<Integer, Integer> getInputRankLabelledAtomsReactant() {
        return inputRankLabelledAtomsReactant;
    }

    /**
     * Old Atom Rank in the product mapped to new Rank
     *
     * @return the inputRankLabelledAtomsProduct
     */
    public Map<Integer, Integer> getInputRankLabelledAtomsProduct() {
        return inputRankLabelledAtomsProduct;
    }


    /**
     * @contact Syed Asad Rahman, BioInception.
     * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     */
    public static abstract class Debugger extends BasicDebugger {

        private final static ILoggingTool LOGGER
                = createLoggingTool(Debugger.class);

        /**
         * Prints reactant and product atom container in the matrix
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printMatrixAtomContainer(Holder mh, List<String> EdMap, List<String> PdMap) {
            try {
                ReactionContainer _rSTMap = mh.getReactionContainer();
                StringBuilder sb = new StringBuilder();
                sb.append("<--------Atom Size in the Container-------->").append(NEW_LINE);
                for (int i = 0; i < EdMap.size(); i++) {
                    sb.append("Educt ").append(EdMap.get(i)).append(" : ").append(_rSTMap.getEduct(i).getAtomCount()).append(NEW_LINE);
                    if (!_rSTMap.getEduct(i).isEmpty()) {
                        CDKSMILES sm = new CDKSMILES(_rSTMap.getEduct(i), true, false);
                        sb.append("SMILES: ").append(sm.getCanonicalSMILES()).append(NEW_LINE);
                    }
                    printAtoms(_rSTMap.getEduct(i));
                }
                sb.append(NEW_LINE);
                for (int i = 0; i < PdMap.size(); i++) {
                    sb.append("Product ").append(PdMap.get(i)).append(" : ").append(_rSTMap.getProduct(i).getAtomCount()).append(NEW_LINE);
                    if (!_rSTMap.getProduct(i).isEmpty()) {
                        CDKSMILES sm = new CDKSMILES(_rSTMap.getProduct(i), true, false);
                        sb.append("SMILES: ").append(sm.getCanonicalSMILES()).append(NEW_LINE);
                    }
                    printAtoms(_rSTMap.getProduct(i));
                }
                LOGGER.debug(sb.toString());
            } catch (IOException | CDKException | CloneNotSupportedException ex) {
                LOGGER.error(SEVERE, null, ex);
            }
        }

        /**
         * Prints Clique Matrix
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printCliqueMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {

            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);
            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;
                sb.append("Clique Matrix").append(NEW_LINE);
                sb.append("\t\t");
                for (int j = 0; j < PdMap.size(); j++) {
                    sb.append(" ").append(PdMap.get(j)).append(":(").append(reactionStructureInformationContainer.getProduct(j).getAtomCount()).append(")");
                }
                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i < EdMap.size(); i++) {
                    sb.append(" ").append(EdMap.get(i)).append(":(").append(reactionStructureInformationContainer.getEduct(i).getAtomCount()).append(")");
                    for (int j = 0; j < PdMap.size(); j++) {
                        val = mh.getCliqueMatrix().getValue(i, j);
                        result = format.format(val);
                        sb.append("   ").append(result);
                    }
                    sb.append(NEW_LINE);
                }
            } catch (IOException | CDKException e) {
                LOGGER.debug("Parser Error" + e);
            }
            LOGGER.debug(sb.toString());
        }

        /**
         * Prints Similarity Matrix
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printSimMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
            ReactionContainer reactionStructureInformationContainer = mh.getReactionContainer();
            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);
            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;
                sb.append("Similarity Matrix").append(NEW_LINE);
                sb.append("\t\t");
                for (int j = 0; j < PdMap.size(); j++) {
                    sb.append(" ").append(PdMap.get(j)).append(":(").append(reactionStructureInformationContainer.getProduct(j).getAtomCount()).append(")");
                }
                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i < EdMap.size(); i++) {
                    sb.append(" ").append(EdMap.get(i)).append(":(").append(reactionStructureInformationContainer.getEduct(i).getAtomCount()).append(")");
                    for (int j = 0; j < PdMap.size(); j++) {
                        val = mh.getGraphSimilarityMatrix().getValue(i, j);
                        result = format.format(val);
                        sb.append("   ").append(result);
                    }
                    sb.append(NEW_LINE);
                }
            } catch (IOException | CDKException e) {
                LOGGER.debug("Parser Error" + e);
            }
            LOGGER.debug(sb.toString());

        }

        /**
         *
         * @param winner
         * @param EdMap
         * @param PdMap
         */
        protected void printFlagMatrix(MappingChecks.ChooseWinner winner, List<String> EdMap, List<String> PdMap) {

            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);
            boolean[][] FlagMatrix = winner.getFlagMatrix();
            sb.append("Flag Matrix").append(NEW_LINE);
            sb.append("\t\t");
            PdMap.forEach((PdMap1) -> {
                sb.append("  ").append(PdMap1).append(" ");
            });

            sb.append(NEW_LINE);
            for (int i = 0; i < EdMap.size(); i++) {
                sb.append(" ").append(EdMap.get(i));
                for (int j = 0; j < PdMap.size(); j++) {
                    if (FlagMatrix[i][j]) {
                        sb.append("      ").append(1).append("  ");
                    } else {
                        sb.append("      ").append(0).append("  ");
                    }

                }
                sb.append(NEW_LINE);
            }
            LOGGER.debug(sb.toString());

        }

        /**
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printStereoMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
            EBIMatrix StereoMatrix = mh.getStereoMatrix();

            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);

            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;

                sb.append("Stereo Matrix").append(NEW_LINE);
                sb.append("\t\t");

                PdMap.forEach((PdMap1) -> {
                    sb.append(" ").append(PdMap1);
                });

                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i
                        < EdMap.size(); i++) {
                    sb.append(" ").append(EdMap.get(i));
                    for (int j = 0; j
                            < PdMap.size(); j++) {
                        val = StereoMatrix.getValue(i, j);
                        result
                                = format.format(val);
                        sb.append("   ").append(result);
                    }

                    sb.append(NEW_LINE);
                }

            } catch (Exception e) {
                LOGGER.debug("Parser Error" + e);
            }

            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printFragmentMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
            EBIMatrix fragmentMatrix = mh.getFragmentMatrix();

            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);

            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;

                sb.append("Fragment Matrix").append(NEW_LINE);
                sb.append("\t\t");

                PdMap.forEach((PdMap1) -> {
                    sb.append(" ").append(PdMap1);
                });

                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i
                        < EdMap.size(); i++) {
                    sb.append(" ").append(EdMap.get(i));
                    for (int j = 0; j
                            < PdMap.size(); j++) {
                        val = fragmentMatrix.getValue(i, j);
                        result
                                = format.format(val);
                        sb.append("   ").append(result);
                    }

                    sb.append(NEW_LINE);
                }

            } catch (Exception e) {
                LOGGER.debug("Parser Error" + e);
            }

            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printCarbonMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
            EBIMatrix carbonMatrix = mh.getCarbonOverlapMatrix();

            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);

            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;

                sb.append("Fragment Matrix").append(NEW_LINE);
                sb.append("\t\t");

                PdMap.forEach((PdMap1) -> {
                    sb.append(" ").append(PdMap1);
                });

                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i < EdMap.size(); i++) {
                    sb.append(" ").append(EdMap.get(i));
                    for (int j = 0; j
                            < PdMap.size(); j++) {
                        val = carbonMatrix.getValue(i, j);
                        result = format.format(val);
                        sb.append("   ").append(result);
                    }

                    sb.append(NEW_LINE);
                }

            } catch (Exception e) {
                LOGGER.debug("Parser Error" + e);
            }

            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param mh
         * @param EdMap
         * @param PdMap
         */
        protected void printEnergyMatrix(Holder mh, List<String> EdMap, List<String> PdMap) {
            EBIMatrix energyMatrixProfile = mh.getEnergyMatrix();

            StringBuilder sb = new StringBuilder();
            sb.append(NEW_LINE);
            sb.append("********* MATRIX **********").append(NEW_LINE);

            try {
                NumberFormat format = new DecimalFormat("0.00");
                String result;

                sb.append("Energy Matrix").append(NEW_LINE);
                sb.append("\t\t");

                PdMap.forEach((PdMap1) -> {
                    sb.append("\t").append(PdMap1);
                });

                sb.append(NEW_LINE);
                double val;
                for (int i = 0; i
                        < EdMap.size(); i++) {
                    sb.append("\t").append(EdMap.get(i));
                    for (int j = 0; j
                            < PdMap.size(); j++) {
                        val = energyMatrixProfile.getValue(i, j);
                        result = format.format(val);
                        sb.append("\t").append(result);
                    }

                    sb.append(NEW_LINE);
                }

            } catch (Exception e) {
                LOGGER.debug("Parser Error" + e);
            }

            LOGGER.debug(sb.toString());
        }

        /**
         * Print Graph matching solutions
         *
         * @param comparison
         * @param mol1
         * @param mol2
         */
        protected void printGraphMatching(IAtomMapping comparison, IAtomContainer mol1, IAtomContainer mol2) {
            int count_final_sol = 0;
            StringBuilder sb = new StringBuilder();
            sb.append("Output of the final Mappings: ").append(NEW_LINE);
            sb.append("Mol1: ").append(mol1.getID()).append(NEW_LINE);
            sb.append("Mol2: ").append(mol2.getID()).append(NEW_LINE);
            try {
                if (comparison.getMappingCount() > 0) {

                    for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                        int final_solution_size = final_solution.getCount();
                        sb.append("Final mapping Nr. ").append(++count_final_sol)
                                .append(" Size:").append(final_solution_size).append(NEW_LINE);

                        final int solIndex = count_final_sol;
                        final_solution.getMappingsByAtoms().entrySet().forEach((mapping) -> {
                            IAtom eAtom = mapping.getKey();
                            IAtom pAtom = mapping.getValue();

                            sb.append(mol1.indexOf(eAtom) + 1).append(" ").append(mol2.indexOf(pAtom) + 1).append(NEW_LINE);

                            sb.append(eAtom.getSymbol()).append(" ")
                                    .append(pAtom.getSymbol()).append(NEW_LINE);
                        });
                        sb.append("").append(NEW_LINE);

                        sb.append("Stereo Match: ").append(comparison.getStereoScore(count_final_sol - 1)).append(NEW_LINE);
                        sb.append("Stereo different: ").append(comparison.isStereoMisMatch()).append(NEW_LINE);
                        sb.append("Fragment Size: ").append(comparison.getFragmentSize(count_final_sol - 1)).append(NEW_LINE);
                    }

                    sb.append("").append(NEW_LINE);
                }
            } catch (Exception ex) {
                LOGGER.debug("Parser Error" + ex);
            }
            LOGGER.debug(sb.toString());
        }

        /**
         *
         * @param outPutFileName
         * @param query
         * @param target
         * @param smsd
         */
        protected void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, BaseMapping smsd) {

            ImageGenerator imageGenerator = new ImageGenerator();

            ////set the format right for the Tanimoto score (only two digits printed)
            NumberFormat nf = getInstance();
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);
            LOGGER.debug("Output of the final Mappings: ");
            int counter = 1;
            for (AtomAtomMapping mapping : smsd.getAllAtomMapping()) {

                String tanimoto = nf.format(smsd.getTanimotoSimilarity());
                String stereo = "NA";
                if (smsd.getStereoScore(counter - 1) != null) {
                    stereo = nf.format(smsd.getStereoScore(counter - 1));
                }
                String label = "Scores [" + "Tanimoto: " + tanimoto + ", Stereo: " + stereo + "]";
                try {
                    imageGenerator.addImages(query, target, label, mapping);
                } catch (Exception ex) {
                    LOGGER.error(SEVERE, null, ex);
                }
                counter++;
            }
            String filePNG = getProperty("user.dir") + separator + outPutFileName;
            imageGenerator.createImage(filePNG, "Query", "Target");
        }

        /**
         *
         * @param mh
         * @param substrateIndex
         * @param productIndex
         * @throws CloneNotSupportedException
         * @throws IOException
         * @throws CDKException
         */
        protected void printSMILES(Holder mh, int substrateIndex, int productIndex)
                throws CloneNotSupportedException, IOException, CDKException {
            ReactionContainer reactionStructureInformation = mh.getReactionContainer();
            String cdkSmilesE = new CDKSMILES(reactionStructureInformation.getEduct(substrateIndex), false, false).getCanonicalSMILES();
            String cdkSmilesP = new CDKSMILES(reactionStructureInformation.getProduct(productIndex), false, false).getCanonicalSMILES();

            StringBuilder sb = new StringBuilder();
            sb.append("A: ").append(reactionStructureInformation.getEduct(substrateIndex).getID()).append(" ").append(cdkSmilesE)
                    .append(" B: ").append(reactionStructureInformation.getProduct(productIndex).getID()).append(" ").append(cdkSmilesP).append(NEW_LINE);

            sb.append("A: ").append(reactionStructureInformation.getEduct(substrateIndex).getAtomCount())
                    .append(" B: ").append(reactionStructureInformation.getProduct(productIndex).getAtomCount()).append(NEW_LINE);

            sb.append(" GetValue: ").append(mh.getGraphSimilarityMatrix().getValue(substrateIndex, productIndex))
                    .append(", ").append(mh.getStereoMatrix().getValue(substrateIndex, productIndex));

            LOGGER.debug(sb.toString());
        }
    }



    /**
     * @Author: Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
     * @Date: 2009/06/3
     * @Revision: 1.10
     */
    public static class MappingHandler extends BasicDebugger {

        /**
         *
         * @param MappedReaction
         */
        public static void cleanMapping(IReaction MappedReaction) {
            int count = MappedReaction.getMappingCount();
            for (int i = count - 1; i >= 0; i--) {
                MappedReaction.removeMapping(i);
            }

            for (int eMol = 0; eMol < MappedReaction.getReactantCount(); eMol++) {
                IAtomContainer eMolecule = MappedReaction.getReactants().getAtomContainer(eMol);
                for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {

                    IAtom atomEMap = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                    String atomLabel = Integer.toString(-1);
                    atomEMap.setFlag(MAPPED, false);
                    atomEMap.setID(atomLabel);
                }
            }

            for (int pMol = 0; pMol < MappedReaction.getProductCount(); pMol++) {
                IAtomContainer pMolecule = MappedReaction.getProducts().getAtomContainer(pMol);
                for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                    IAtom atomPMap = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                    String atomLabel = Integer.toString(-1);
                    atomPMap.setFlag(MAPPED, false);
                    atomPMap.setID(atomLabel);
                }
            }
        }

        /**
         *
         * @param expLabReaction
         * @param MappedReaction
         * @param counter
         * @return
         */
        protected static int setMappingFlags(IReaction expLabReaction, IReaction MappedReaction, int counter) {
            IAtomContainerSet expEductSet = expLabReaction.getReactants();
            IAtomContainerSet expProductSet = expLabReaction.getProducts();

            for (IMapping map : expLabReaction.mappings()) {

                IAtom I_Atom = (IAtom) map.getChemObject(0);
                IAtom J_Atom = (IAtom) map.getChemObject(1);

                if (I_Atom != null && J_Atom != null) {

                    /*
                    *******************************
                    * Mapping the Reactants ******************************
                     */
                    boolean eFlag = false;
                    IAtom firstAtom = null;
                    IAtom secondAtom = null;
                    for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                        IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                        for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                            if (I_Atom.getID().trim().equalsIgnoreCase(eMolecule.getAtom(eAtom).getID().trim())) {

                                String atomLabel = Integer.toString(counter);

                                firstAtom = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                                firstAtom.setID(atomLabel);
                                firstAtom.setFlag(MAPPED, true);
                                eFlag = true;
                                break;
                            }
                        }

                        if (eFlag) {
                            break;
                        }

                    }
                    /*
                    *******************************
                    * Mapping the Products ******************************
                     */
                    boolean pFlag = false;
                    for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                        IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                        for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                            if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {

                                String atomLabel = Integer.toString(counter);

                                secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                                secondAtom.setID(atomLabel);
                                secondAtom.setFlag(MAPPED, true);
                                IMapping mappingObject = MappedReaction.getBuilder().newInstance(IMapping.class, firstAtom, secondAtom);
                                MappedReaction.addMapping(mappingObject);
                                counter++;
                                pFlag = true;
                                break;
                            }
                        }

                        if (pFlag) {
                            break;
                        }
                    }
                }
            }
            return counter;
        }

        /**
         *
         * @param MappedReaction
         * @param ReactionWithUniqueSTOICHIOMETRY
         * @param coreMappedReaction
         * @param counter
         * @return
         */
        protected static int setMappingFlags(IReaction MappedReaction, IReaction ReactionWithUniqueSTOICHIOMETRY, IReaction coreMappedReaction, int counter) {
            IAtomContainerSet expEductSet = ReactionWithUniqueSTOICHIOMETRY.getReactants();
            IAtomContainerSet expProductSet = ReactionWithUniqueSTOICHIOMETRY.getProducts();

            for (IMapping map : coreMappedReaction.mappings()) {

                IAtom I_Atom = (IAtom) map.getChemObject(0);
                IAtom J_Atom = (IAtom) map.getChemObject(1);


                if (I_Atom != null && J_Atom != null) {

                    /*
                     * Mapping the Reactants 
                     */
                    boolean eFlag = false;
                    IAtom firstAtom = null;
                    IAtom secondAtom = null;
                    for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                        IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                        for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                            if (I_Atom.getID().trim().equalsIgnoreCase(eMolecule.getAtom(eAtom).getID().trim())) {

                                String atomLabel = Integer.toString(counter);

                                firstAtom = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                                firstAtom.setID(atomLabel);
                                firstAtom.setFlag(MAPPED, true);
                                eFlag = true;
                                break;
                            }
                        }

                        if (eFlag) {
                            break;
                        }

                    }
                    /*
                     * Mapping the Products 
                     */
                    boolean pFlag = false;
                    for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                        IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                        for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                            if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {

                                String atomLabel = Integer.toString(counter);


                                secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                                secondAtom.setID(atomLabel);
                                secondAtom.setFlag(MAPPED, true);
                                counter++;
                                pFlag = true;
                                break;
                            }
                        }

                        if (pFlag) {
                            break;
                        }

                    }
                }

            }

            return counter;
        }
    }


}
