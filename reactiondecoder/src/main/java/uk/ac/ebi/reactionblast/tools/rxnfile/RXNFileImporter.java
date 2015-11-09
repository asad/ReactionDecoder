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

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import uk.ac.ebi.reactionblast.tools.EBIMolSplitter;

/**
 * @author Syed Asad Rahman, EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
public class RXNFileImporter {

    private IReaction reaction;
    protected IAtomContainer atomContainer;
    int nProducts;
    int nReactants;
    public final static String RXN = "$RXN";

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
            reaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);

            MDLRXNV2000Reader reader = new MDLRXNV2000Reader(RXNFile, Mode.STRICT);
            reaction = reader.read(reaction);
            reader.close();
        } catch (CDKException cdkerr) {
            System.out.println("Error: only RXN V2000 file format is "
                    + "supported by this Software");
            System.err.println("Error: " + cdkerr);
        } catch (FileNotFoundException e) {
            System.err.println("Error: RXN File not found" + e);
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

        IAtomContainerSet SplitMoleculeList = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

        for (int j = 0; j < BigMol.getAtomContainerCount(); j++) {
            IAtomContainer molecules = BigMol.getAtomContainer(j);
//            System.out.println("Big Mol: " + BigMol.getAtomContainerCount() + " Big Mol Atom Cont:" + molecules.getAtomCount());
            if (!ConnectivityChecker.isConnected(molecules)) {
                IAtomContainerSet splitMol = EBIMolSplitter.splitMolecules(molecules);
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
        System.out.print("Atom: ");
        for (IAtom a : mol.atoms()) {
            System.out.print(a.getSymbol());
        }
        System.out.println();
        System.out.println();
    }

    private void printRadicalAtoms(IAtomContainer mol) {
        System.out.print("Atom: ");
        for (Integer a = 0; a < mol.getSingleElectronCount(); a++) {
            System.out.print(mol.getSingleElectron(a).getAtom().getSymbol());
            System.out.print(mol.getSingleElectron(a).getElectronCount());
        }
        System.out.println();
        System.out.println();
    }
}
