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
package uk.ac.ebi.reactionblast.interfaces;

import java.util.List;


/**
 * A source of reactions or molecules, for example a database or filesystem.
 * 
 * @author maclean
 * @param <T>
 *
 */
public interface IDataSource<T> {
    
    /**
     * Get the reaction with this ID.
     * 
     * @param id the identifier for this reaction.
     * @return a reaction
     */
    T get(String id);
    
    /**
     * Get all the reactions in the data source.
     * 
     * @return an iterable for all the reactions
     */
    Iterable<T> getAll();
    
    /**
     * Set the transformation to apply to the object before returning it.
     * 
     * @param transformation
     */
    void setTransformation(ITransformation<T> transformation);
    
    /**
     * Get a list of the IDs in this data source.
     * 
     * @return a list of ID strings suitable for passing to get(String id)
     */
    List<String> getIDList();
    
    /**
     * Close the data source.
     */
    void close();
  
}
