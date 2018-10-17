/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.io;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import static java.lang.System.getProperty;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import uk.ac.ebi.centres.Digraph;
import uk.ac.ebi.centres.Ligand;

/**
 * Allows a digraph to be created
 *
 * @author John May
 * @param <A>
 */
public abstract class CytoscapeWriter<A> implements Closeable {

    static final String NEW_LINE = getProperty("line.separator");
    private final Digraph<A> digraph;
    private Writer sif;
    private File folder;
    private final Map<String, Map<String, String>> attributes = new HashMap<>();

    /**
     *
     * @param folder
     * @param digraph
     * @throws IOException
     */
    public CytoscapeWriter(File folder, Digraph<A> digraph) throws IOException {

        this.digraph = digraph;

        if (folder.exists() && !folder.isDirectory()) {
            throw new IllegalArgumentException("Folder should be a directory");
        }

        if (!folder.exists() && !folder.mkdirs()) {
            throw new IllegalArgumentException("Unable to create folder");
        }

        this.folder = folder;
        this.sif = new FileWriter(new File(folder, folder.getName().replace(" ", "-") + ".sif"));

    }

    /**
     *
     * @throws IOException
     */
    public void writeSif() throws IOException {
        write(digraph.getProximal(), "1");
    }

    /**
     *
     * @throws IOException
     */
    public void writeAttributes() throws IOException {
        // do nothing
        for (Map.Entry<String, Map<String, String>> entry : attributes.entrySet()) {
            try (FileWriter attributeWriter = new FileWriter(new File(folder, entry.getKey() + ".noa"))) {
                attributeWriter.write(entry.getKey().replaceAll(" ", ".") + " (class=String)" + NEW_LINE);
                for (Map.Entry<String, String> nodeEntry : entry.getValue().entrySet()) {
                    attributeWriter.write(nodeEntry.getKey() + " = " + nodeEntry.getValue() + NEW_LINE);
                }
            }
        }
    }

    private void write(List<Ligand<A>> ligands, String sourceId) throws IOException {

        for (int i = 0; i < ligands.size(); i++) {

            Ligand<A> ligand = ligands.get(i);

            String targetId = sourceId + (Integer.toString(i + 1));

            sif.write(sourceId + "\t->\t" + targetId + "\n");

            // invert map properties
            Map<String, String> map = new HashMap<>();
            mapAttributes(ligand.getAtom(), map);
            map.entrySet().stream().map((Map.Entry<String, String> e) -> {
                if (!attributes.containsKey(e.getKey())) {
                    attributes.put(e.getKey(), new HashMap<>());
                }
                return e;
            }).forEachOrdered((e) -> {
                attributes.get(e.getKey()).put(targetId, e.getValue());
            });

            write(ligands.get(i).getLigands(), targetId);

        }

    }

    /**
     *
     * @param atom
     * @param map
     */
    public abstract void mapAttributes(A atom, Map<String, String> map);

    @Override
    public void close() throws IOException {
        sif.close();
    }
}
