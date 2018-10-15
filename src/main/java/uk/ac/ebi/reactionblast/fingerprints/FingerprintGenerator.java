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
package uk.ac.ebi.reactionblast.fingerprints;

import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import static org.openscience.cdk.fingerprint.CircularFingerprinter.CLASS_ECFP4;
import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import static org.openscience.cdk.graph.ConnectivityChecker.isConnected;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFingerprintGenerator;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FingerprintGenerator implements IFingerprintGenerator {

    private final static ILoggingTool LOGGER
            = createLoggingTool(FingerprintGenerator.class);

    /**
     * Size of the fingerprint
     *
     * @return
     */
    public static int getFingerprinterSize() {
        return new CircularFingerprinter(CLASS_ECFP4).getSize();
    }

    //define the FINGERPRINT_SIZE of the fingerprint
    //NOTE: this should be a multiple of 64 and preferably not 1024 or 2048
    //as for these values we often get the random numbers for one-atom or
    //two-atom paths the same!
    final CircularFingerprinter fingerprinter;

    /**
     *
     */
    public FingerprintGenerator() {
        fingerprinter = new CircularFingerprinter(CLASS_ECFP4);
    }

    /**
     *
     * @param mol
     * @return
     * @throws CDKException
     */
    @Override
    public synchronized BitSet getFingerprint(IAtomContainer mol) throws CDKException {
        if (!has2DCoordinates(mol)) {
            StructureDiagramGenerator structureDiagramGenerator = new StructureDiagramGenerator();
            structureDiagramGenerator.setMolecule(mol, true);
            if (isConnected(mol)) {
                structureDiagramGenerator.generateCoordinates();
                mol = structureDiagramGenerator.getMolecule();
            } else {
                LOGGER.debug("Disconnected components needs to be layout separately");
            }
        }
        return fingerprinter.getBitFingerprint(mol).asBitSet();
    }

}
