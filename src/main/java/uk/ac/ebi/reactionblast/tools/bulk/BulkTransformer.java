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
package uk.ac.ebi.reactionblast.tools.bulk;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.interfaces.IDataSource;
import uk.ac.ebi.reactionblast.interfaces.IDataStore;
import uk.ac.ebi.reactionblast.interfaces.ITransformation;

/**
 * Transforms a number of reactions (or molecules in reactions) in turn.
 *
 * @author maclean
 *
 */
public class BulkTransformer {

    /**
     * The transformation, or series of transformations, to apply.
     */
    private final ITransformation transformation;

    /**
     *
     * @param transformation
     */
    public BulkTransformer(ITransformation transformation) {
        this.transformation = transformation;
    }

    /**
     * Transform all the entries in the data source.
     *
     * @param dataSource
     * @param dataStore
     */
    public void transform(IDataSource dataSource, IDataStore dataStore) {
        switch (transformation.getTargetType()) {
            case REACTION:
                transformReactions(dataSource, dataStore);
                break;
            case MOLECULE:
                transformMolecules(dataSource, dataStore);
                break;
            default:
        }
    }

    private void transformReactions(IDataSource<IReaction> dataSource, IDataStore dataStore) {
        IDataSource<IReaction> rxnDataSource = dataSource;
        for (IReaction reaction : rxnDataSource.getAll()) {
            if (reaction != null) {
                IReaction transformedReaction
                        = (IReaction) transformation.transform(reaction);
                if (transformedReaction != null) {
                    dataStore.store(transformedReaction);
                }
            }
        }
    }

    private void transformMolecules(IDataSource<IAtomContainer> dataSource, IDataStore dataStore) {
        IDataSource<IAtomContainer> molDataSource = dataSource;
        for (IAtomContainer molecule : molDataSource.getAll()) {
            if (molecule != null) {
                IAtomContainer transformedMolecule
                        = (IAtomContainer) transformation.transform(molecule);
                if (transformedMolecule != null) {
                    dataStore.store(transformedMolecule);
                }
            }
        }
    }
}
