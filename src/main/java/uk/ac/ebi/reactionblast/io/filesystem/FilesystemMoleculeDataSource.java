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


import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import uk.ac.ebi.reactionblast.interfaces.IDataSource;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class FilesystemMoleculeDataSource implements IDataSource<IAtomContainer> {
    
    private File moleculeDir;
    
    private MDLV2000Reader molReader;
    
    private ITransformation<IAtomContainer> transformation;
    
    /**
     *
     * @param moleculePath
     */
    public FilesystemMoleculeDataSource(String moleculePath) {
        if (moleculePath != null) {
            this.moleculeDir = new File(moleculePath);
            molReader = new MDLV2000Reader();
        }
    }

    @Override
    public void setTransformation(ITransformation<IAtomContainer> transformation) {
        this.transformation = transformation;
    }

    @Override
    public IAtomContainer get(String id) {
        File molFile = new File(moleculeDir, id + ".mol");
        try {
            molReader.setReader(new FileReader(molFile));
            IAtomContainer mol = molReader.read(new AtomContainer());
            if (transformation == null) {
                return mol;
            } else {
                return transformation.transform(mol);
            }
        } catch (CDKException | FileNotFoundException c) {
            return null;
        }
    }

    @Override
    public Iterable<IAtomContainer> getAll() {
        final String[] fileNames = moleculeDir.list();
        final int numberOfReactions = fileNames.length;
        return new Iterable<IAtomContainer>() {

            @Override
            public Iterator<IAtomContainer> iterator() {
                return new Iterator<IAtomContainer>() {
                    
                    private int currentIndex;

                    @Override
                    public boolean hasNext() {
                        return currentIndex < numberOfReactions;
                    }

                    @Override
                    public IAtomContainer next() {
                        String filename = fileNames[currentIndex];
                        File molFile = new File(moleculeDir, filename);
                        currentIndex++;
                        try {
                            IAtomContainer mol;
                            try (FileReader reader = new FileReader(molFile)) {
                                molReader.setReader(reader);
                                mol = molReader.read(new AtomContainer());
                                mol.setID(filename.substring(0, filename.indexOf('.')));
                            }
                            if (transformation == null) {
                                return mol;
                            } else {
                                return transformation.transform(mol); 
                            }
                        } catch (FileNotFoundException e) {
                            e.printStackTrace();
                            return null;
                        } catch (CDKException | IOException e) {
                            e.printStackTrace();
                            return null;                        
                        }
                    }

                    @Override
                    public void remove() {}
                    
                };
            }
        };
    }

    @Override
    public List<String> getIDList() {
        List<String> ids = new ArrayList<>();
        for (String filename : moleculeDir.list()) {
            ids.add(filename.substring(0, filename.indexOf('.')));
        }
        return ids;
    }

    @Override
    public void close() {
        try {
            molReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
