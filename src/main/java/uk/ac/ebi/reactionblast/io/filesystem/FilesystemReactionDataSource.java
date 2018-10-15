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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.interfaces.IDataSource;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FilesystemReactionDataSource implements IDataSource<IReaction> {

    private File reactionDir;

    private MDLRXNV2000Reader rxnReader;

    private ITransformation<IReaction> transformation;

    /**
     *
     * @param reactionPath
     */
    public FilesystemReactionDataSource(String reactionPath) {
        if (reactionPath != null) {
            this.reactionDir = new File(reactionPath);
            rxnReader = new MDLRXNV2000Reader();
        }
    }

    @Override
    public void setTransformation(ITransformation<IReaction> transformation) {
        this.transformation = transformation;
    }

    @Override
    public IReaction get(String id) {
        File reactionFile = new File(reactionDir, id + ".rxn");
        try {
            rxnReader.setReader(new FileReader(reactionFile));
            IReaction reaction = rxnReader.read(new Reaction());
            if (transformation == null) {
                return reaction;
            } else {
                return transformation.transform(reaction);
            }
        } catch (CDKException | FileNotFoundException c) {
            c.printStackTrace();
            return null;
        }
    }

    @Override
    public Iterable<IReaction> getAll() {
        final String[] fileNames = reactionDir.list();
        final int numberOfReactions = fileNames.length;
        return new Iterable<IReaction>() {

            @Override
            public Iterator<IReaction> iterator() {
                return new Iterator<IReaction>() {

                    private int currentIndex;

                    @Override
                    public boolean hasNext() {
                        return currentIndex < numberOfReactions;
                    }

                    @Override
                    public IReaction next() {
                        String filename = fileNames[currentIndex];
                        File rxnFile = new File(reactionDir, filename);
                        currentIndex++;
                        try {
                            IReaction reaction;
                            try (FileReader reader = new FileReader(rxnFile)) {
                                rxnReader.setReader(reader);
                                reaction = rxnReader.read(new Reaction());
                                reaction.setID(filename.substring(0, filename.indexOf('.')));
                            }
                            return reaction;
                        } catch (FileNotFoundException e) {
                            e.printStackTrace();
                            return null;
                        } catch (CDKException | IOException e) {
                            e.printStackTrace();
                            return null;
                        }
                    }

                    @Override
                    public void remove() {
                    }

                };
            }
        };
    }

    @Override
    public List<String> getIDList() {
        List<String> ids = new ArrayList<>();
        for (String filename : reactionDir.list()) {
            ids.add(filename.substring(0, filename.indexOf('.')));
        }
        return ids;
    }

    @Override
    public void close() {
        try {
            rxnReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
