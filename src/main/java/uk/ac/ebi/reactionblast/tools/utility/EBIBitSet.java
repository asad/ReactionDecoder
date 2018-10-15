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
package uk.ac.ebi.reactionblast.tools.utility;

import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class EBIBitSet extends BitSet {

    private static final long serialVersionUID = 3997698588986878753L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(EBIBitSet.class);

    /**
     * 
     * @param fp1 fingerprint
     * @param fp2 fingerprint
     * @return
     * @throws CDKException 
     */
    public static BitSet append(BitSet fp1, BitSet fp2) throws CDKException {
        BitSet append=null;
        
        if(fp1.size()> 0 && fp2.size() > 0){
            append=new BitSet(fp1.size()+fp2.size());
            
            int index=0;
            for(int i=0;i<=fp1.size();i++){
                append.set(index++, fp1.get(i));
            }
            
            for(int j=0; j<=fp2.size();j++){
                append.set(index++, fp2.get(j));
            }
            
        }else{
            throw new CDKException("BitIndex <0: ");
        }
        
        return append;
    }
    
    /**
     *
     */
    public EBIBitSet() {
        super();
    }

    /**
     *
     * @param nbits
     */
    public EBIBitSet(int nbits) {
        super(nbits);
    }

    @Override
    public Object clone() {
        return super.clone(); //To change body of generated methods, choose Tools | Templates.
    }
}
