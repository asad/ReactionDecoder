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
package uk.ac.ebi.reactionblast.io.filesystem;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import uk.ac.ebi.reactionblast.interfaces.IDataStore;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000Writer;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FilesystemMoleculeDataStore implements IDataStore<IAtomContainer> {

    private File moleculeDir;

    private MDLV2000Writer molWriter;

    /**
     *
     * @param moleculePath
     */
    public FilesystemMoleculeDataStore(String moleculePath) {
        if (moleculePath != null) {
            this.moleculeDir = new File(moleculePath);
            molWriter = new MDLV2000Writer();
        }
    }

    @Override
    public void store(IAtomContainer molecule) {
        String id = molecule.getID();
        File file = new File(moleculeDir, id + ".mol");
        try {
            FileWriter writer = new FileWriter(file);
            molWriter = new MDLV2000Writer(writer);
            molWriter.write(molecule);
            molWriter.close();
        } catch (CDKException | IOException e) {
            e.printStackTrace();
        }

    }

}
