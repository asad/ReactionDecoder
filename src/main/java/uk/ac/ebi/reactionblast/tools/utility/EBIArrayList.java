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

import java.util.ArrayList;
import java.util.Collection;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 *
 * @param <E>
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class EBIArrayList<E> extends ArrayList<E> {

    private static final long serialVersionUID = 7683452581122892189L;
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(EBIArrayList.class);

    /**
     *
     * @param L1
     * @param L2
     * @return
     */
    @SuppressWarnings({"unchecked", "unchecked"})
    public static ArrayList append(ArrayList L1, ArrayList L2) {

        ArrayList temp1 = new ArrayList(L1);
        ArrayList temp2 = new ArrayList(L2);
        temp1.trimToSize();
        temp2.trimToSize();
        int initialCapacity = temp1.size() + temp2.size();

        ArrayList mergedList = new ArrayList(initialCapacity);

        int index = 0;
        for (int i = 0; i <= temp1.size(); i++) {
            mergedList.set(index++, temp1.get(i));
        }

        for (int j = 0; j <= temp2.size(); j++) {
            mergedList.set(index++, temp1.get(j));
        }
        return mergedList;
    }

    /**
     *
     * @param initialCapacity
     */
    public EBIArrayList(int initialCapacity) {
        super(initialCapacity);
    }

    /**
     *
     */
    public EBIArrayList() {
        super();
    }

    /**
     *
     * @param c
     */
    @SuppressWarnings("unchecked")
    public EBIArrayList(Collection<? extends E> c) {
        super(c);
    }

    @Override
    public Object clone() {
        return super.clone(); //To change body of generated methods, choose Tools | Templates.
    }
}
