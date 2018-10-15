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
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.interfaces.IDataStore;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000RXNWriter;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FilesystemReactionDataStore implements IDataStore<IReaction> {

    private File reactionDir;

    private MDLV2000RXNWriter rxnWriter;

    /**
     *
     * @param reactionPath
     */
    public FilesystemReactionDataStore(String reactionPath) {
        if (reactionPath != null) {
            this.reactionDir = new File(reactionPath);
            if (!reactionDir.exists()) {
                reactionDir.mkdir();
            }
            rxnWriter = new MDLV2000RXNWriter();
        }
    }

    @Override
    public void store(IReaction reaction) {
        String id = reaction.getID();
        File file = new File(reactionDir, id + ".rxn");
        try {
            FileWriter writer = new FileWriter(file);
            rxnWriter = new MDLV2000RXNWriter(writer);
//            rxnWriter.setWriter(writer);
            rxnWriter.write(reaction);
//            writer.flush();
//            writer.close();
            rxnWriter.close();
        } catch (CDKException | IOException e) {
            e.printStackTrace();
        }

    }

}
