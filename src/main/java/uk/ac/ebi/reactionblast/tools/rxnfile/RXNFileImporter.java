/*
 * RXNFileImporter.java
 *
 * Created on 12 January 2007, 10:35
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package uk.ac.ebi.reactionblast.tools.rxnfile;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import static java.lang.System.out;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.io.IChemObjectReader.Mode.STRICT;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.reactionblast.tools.EBIMolSplitter.splitMolecules;

/**
 * @author Syed Asad Rahman, EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
public class RXNFileImporter {

    /**
     *
     */
    public static final String RXN = "$RXN";
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(RXNFileImporter.class);

    private IReaction reaction;

    /**
     *
     */
    protected IAtomContainer atomContainer;
    int nProducts;
    int nReactants;

    /**
     * Creates a new instance of RXNFileImporter
     */
    public RXNFileImporter() {
    }

    /**
     *
     * @param File
     * @throws java.io.IOException
     */
    public void readFile(String File) throws IOException {

        try {

            InputStream RXNFile = new BufferedInputStream(new FileInputStream(File));
            reaction = getInstance().newInstance(IReaction.class);

            try (MDLRXNV2000Reader reader = new MDLRXNV2000Reader(RXNFile, STRICT)) {
                reaction = reader.read(reaction);
            }
        } catch (CDKException cdkerr) {
            out.println("Error: only RXN V2000 file format is "
                    + "supported by this Software");
            LOGGER.debug("Error: " + cdkerr);
        } catch (FileNotFoundException e) {
            LOGGER.debug("Error: RXN File not found" + e);
        }
    }

    /**
     *
     * @return all the reactants in the read reaction
     */
    public IAtomContainerSet getReactants() {
        return getMolecules(reaction.getReactants());
    }

    /**
     *
     * @return all the products in the read reaction
     */
    public IAtomContainerSet getProducts() {
        return getMolecules(reaction.getProducts());

    }

    /**
     *
     * @return the Reaction
     */
    public IReaction getReaction() {
        return reaction;
    }

    private IAtomContainerSet getMolecules(IAtomContainerSet BigMol) {

        IAtomContainerSet SplitMoleculeList = getInstance().newInstance(IAtomContainerSet.class);

        for (int j = 0; j < BigMol.getAtomContainerCount(); j++) {
            IAtomContainer molecules = BigMol.getAtomContainer(j);
//            System.out.println("Big Mol: " + BigMol.getAtomContainerCount() + " Big Mol Atom Cont:" + molecules.getAtomCount());
            if (!isConnected(molecules)) {
                IAtomContainerSet splitMol = splitMolecules(molecules);
                for (int i = 0; i < splitMol.getAtomContainerCount(); i++) {
                    IAtomContainer mol = splitMol.getAtomContainer(i);
//                System.out.println("Split Mol Atom Count: " + mol.getAtomCount());
                    SplitMoleculeList.addAtomContainer(mol);
//                printAtoms(mol);
//                printRadicalAtoms(mol);
                }
            } else {
                SplitMoleculeList.addAtomContainer(molecules);
            }
        }
        return SplitMoleculeList;
    }

    private void printAtoms(IAtomContainer mol) {
        out.print("Atom: ");
        for (IAtom a : mol.atoms()) {
            out.print(a.getSymbol());
        }
        out.println();
        out.println();
    }

    private void printRadicalAtoms(IAtomContainer mol) {
        out.print("Atom: ");
        for (Integer a = 0; a < mol.getSingleElectronCount(); a++) {
            out.print(mol.getSingleElectron(a).getAtom().getSymbol());
            out.print(mol.getSingleElectron(a).getElectronCount());
        }
        out.println();
        out.println();
    }
}
