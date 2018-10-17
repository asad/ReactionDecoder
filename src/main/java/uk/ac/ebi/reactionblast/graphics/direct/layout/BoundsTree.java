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
package uk.ac.ebi.reactionblast.graphics.direct.layout;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import static java.lang.String.format;
import static java.lang.System.getProperty;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * A tree of bounding boxes for objects.
 *
 * @author maclean
 *
 */
public final class BoundsTree implements Iterable<Rectangle2D> {

    static final String NEW_LINE = getProperty("line.separator");
    private Rectangle2D root;

    private String rootLabel;

    private Map<String, Rectangle2D> childMap;

    /**
     * Make an empty instance, with the specified label.
     *
     * @param rootLabel
     */
    public BoundsTree(String rootLabel) {
        root = null;    // important : the first added bounds must replace root
        this.rootLabel = rootLabel;
        childMap = new HashMap<>();
    }

    /**
     * Make an instance that contains a single bounds.
     *
     * @param rootLabel
     * @param firstLabel
     * @param firstBox
     */
    public BoundsTree(String rootLabel, String firstLabel, Rectangle2D firstBox) {
        this(rootLabel);
        add(firstLabel, firstBox);
    }

    /**
     *
     * @param rootLabel
     * @param boundsTrees
     */
    public BoundsTree(String rootLabel, BoundsTree... boundsTrees) {
        this(rootLabel);
        for (BoundsTree tree : boundsTrees) {
            add(rootLabel, tree);
        }
    }

    /**
     *
     * @param prefix
     * @return
     */
    public BoundsTree getSubtree(String prefix) {
        BoundsTree subtree = new BoundsTree(rootLabel);
        childMap.keySet().stream().filter((label) -> (label.startsWith(prefix))).forEachOrdered((label) -> {
            subtree.add(label, childMap.get(label));
        });
        return subtree;
    }

    /**
     *
     * @return
     */
    public Rectangle2D getRoot() {
        return this.root;
    }

    /**
     * Adds a rectangular bounds to the tree and updates the root bounds.
     *
     * @param label
     * @param bounds the bounding box of the
     */
    public void add(String label, Rectangle2D bounds) {
        // don't add empty bounding boxes to the root
        boolean isEmpty = (bounds.getCenterX() == 0 && bounds.getCenterY() == 0
                && bounds.getWidth() == 0 && bounds.getHeight() == 0);

        childMap.put(label, bounds);
        if (root == null && !isEmpty) {
            root = new Rectangle2D.Double(
                    bounds.getMinX(), bounds.getMinY(),
                    bounds.getWidth(), bounds.getHeight());
            childMap.put(rootLabel, root);
        } else if (!isEmpty) {
            root.add(bounds);
        }
        if (root != null) {
//            System.out.println("root " + BoundsPrinter.toString(root) + " added " + BoundsPrinter.toString(bounds) + " " + label);
        }
    }

    /**
     *
     * @param label
     * @param point
     */
    public void add(String label, Point2D point) {
        Rectangle2D bounds = new Rectangle2D.Double(point.getX(), point.getY(), 0, 0);
        childMap.put(label, bounds);
        if (root == null) {
            root = new Rectangle2D.Double(
                    bounds.getMinX(), bounds.getMinY(),
                    bounds.getWidth(), bounds.getHeight());
        } else {
            root.add(point);
        }
    }

    /**
     *
     * @param labels
     * @return
     */
    public Rectangle2D getBounds(List<String> labels) {
        Rectangle2D totalBounds = null;
        for (String label : labels) {
            Rectangle2D bounds = get(label);
            if (bounds == null) {
                continue;
            }
            if (totalBounds == null) {
                totalBounds = new Rectangle2D.Double(
                        bounds.getMinX(), bounds.getMinY(),
                        bounds.getWidth(), bounds.getHeight());
            } else {
                totalBounds.add(bounds);
            }
        }
        if (totalBounds == null) {
            // hack - a dummy bounding box, in case of errors
            return new Rectangle2D.Double(0, 0, 100, 100);
        } else {
            return totalBounds;
        }
    }

    // XXX - this method is dangerous, consider removing!
    /**
     *
     * @param root
     */
    public void setRoot(Rectangle2D root) {
        this.root = root;
    }

    /**
     * Add all the members of another tree, prefixing their labels with the
     * supplied label, separated by an underscore. So if the prefix was 'mol1',
     * and the tree had labels {'atom1', 'atom2'}, the resulting bounds would be
     * labeled {'mol1_atom1', 'mol2_atom2'}.
     *
     * @param prefix
     * @param label
     * @param tree
     */
    public void add(String prefix, BoundsTree tree) {
        tree.getBoundLabels().forEach((label) -> {
            add(prefix + "_" + label, tree.get(label));
        });
    }

    /**
     *
     * @return
     */
    public List<String> getBoundLabels() {
        return new ArrayList<>(childMap.keySet());
    }

    /**
     *
     * @param dx
     * @param dy
     */
    public void shift(double dx, double dy) {
        childMap.keySet().stream().map((key) -> childMap.get(key)).forEachOrdered((bounds) -> {
            //            System.out.print(key + " Before : " + BoundsPrinter.toString(bounds));
            bounds.setRect(bounds.getMinX() + dx, bounds.getMinY() + dy,
                    bounds.getWidth(), bounds.getHeight());
//            System.out.println(" After: " + BoundsPrinter.toString(bounds) + " " + dx + " " + dy);
        });
    }

    /**
     * Get the bonding box of the element with this label.
     *
     * @param label
     * @return
     */
    public Rectangle2D get(String label) {
        return childMap.get(label);
    }

    /**
     *
     * @return
     */
    public double getWidth() {
        if (root == null) {
            return 0;
        } else {
            return root.getWidth();
        }
    }

    /**
     *
     * @return
     */
    public double getHeight() {
        if (root == null) {
            return 0;
        } else {
            return root.getHeight();
        }
    }

    @Override
    public Iterator<Rectangle2D> iterator() {
        return childMap.values().iterator();
    }

    /**
     *
     * @param transform
     * @return
     */
    public BoundsTree transform(AffineTransform transform) {
        BoundsTree transformedTree = new BoundsTree(rootLabel);
        childMap.keySet().forEach((key) -> {
            Rectangle2D shape = childMap.get(key);

            // annoyingly, createTransformedShape returns a Path2D! 
            // (so we can't just cast to R2D)...
            transformedTree.add(key, transform.createTransformedShape(shape).getBounds2D());
        });
        return transformedTree;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        childMap.keySet().stream().map((key) -> {
            Rectangle2D rect = get(key);
            sb.append(key).append("=").append(format("[(%2.0f, %2.0f), (%2.0f, %2.0f)]",
                    rect.getMinX(), rect.getMinY(), rect.getMaxX(), rect.getMaxY()));
            return key;
        }).forEachOrdered((_item) -> {
            sb.append(NEW_LINE);
        });
        return sb.toString();
    }

}
