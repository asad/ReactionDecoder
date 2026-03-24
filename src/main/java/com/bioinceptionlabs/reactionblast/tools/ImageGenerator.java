/*
 * Copyright (C) 2007-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
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
package com.bioinceptionlabs.reactionblast.tools;

import java.awt.Image;
import java.awt.image.BufferedImage;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Stub image generator. The full rendering engine has been removed
 * to reduce codebase size. These methods return blank images or no-op.
 * For full rendering, use CDK's DepictionGenerator directly.
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class ImageGenerator {

    private static final ILoggingTool LOGGER = createLoggingTool(ImageGenerator.class);

    /** Return a blank image of the given dimensions. */
    public static Image getBlankImage(int width, int height) {
        return new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    }

    /** No-op: full rendering engine removed. */
    public static void LeftToRightReactionLayoutImageSmall(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public static void LeftToRightReactionCenterImageSmall(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public static void TopToBottomReactionLayoutImageSmall(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public static void LeftToRightReactionLayoutImage(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public static void LeftToRightReactionCenterImage(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public static void TopToBottomReactionLayoutImage(
            IReaction reaction, String name, String dir) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op instance method: full rendering engine removed. */
    public void drawLeftToRightReactionLayout(Object fileOrDir, IReaction reaction, String name) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op instance method: full rendering engine removed. */
    public void drawTopToBottomReactionLayout(Object fileOrDir, IReaction reaction, String name) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public void addImages(org.openscience.cdk.interfaces.IAtomContainer query,
            org.openscience.cdk.interfaces.IAtomContainer target, String label,
            org.openscience.smsd.AtomAtomMapping mapping) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }

    /** No-op: full rendering engine removed. */
    public void createImage(String path, String queryLabel, String targetLabel) {
        LOGGER.debug("Image generation disabled (graphics module removed)");
    }
}
