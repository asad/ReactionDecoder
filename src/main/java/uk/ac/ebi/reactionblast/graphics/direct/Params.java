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

package uk.ac.ebi.reactionblast.graphics.direct;

import java.awt.Color;
import static java.awt.Color.BLUE;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import static uk.ac.ebi.reactionblast.graphics.direct.Params.ArrowType.FORWARD;
import static uk.ac.ebi.reactionblast.graphics.direct.Params.BondStrokeCap.BUTT;
import static uk.ac.ebi.reactionblast.graphics.direct.Params.BondStrokeJoin.MITRE;
import static uk.ac.ebi.reactionblast.graphics.direct.Params.MoleculeAlignMethod.MAX_AXIS;

/**
 *
 * @author asad
 */
public class Params {
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Params.class);

    /**
     *
     */
    public BondStrokeCap bondStrokeCap = BUTT;

    /**
     *
     */
    public BondStrokeJoin bondStrokeJoin = MITRE;

    /**
     *
     */
    public XAlign leftRightAlignment = XAlign.CENTER;

    /**
     *
     */
    public YAlign topBottomAlignment = YAlign.CENTER;

    /**
     *
     */
    public int bondLength = 30;

    /**
     *
     */
    public int borderX = 20;

    /**
     *
     */
    public int borderY = 20;

    /**
     *
     */
    public int plusGap = 20;

    /**
     *
     */
    public int arrowLength = 30;

    /**
     *
     */
    public int arrowGap = 10;

    /**
     *
     */
    public int arrowHeadLength = 10;

    /**
     *
     */
    public boolean drawBounds = false;

    /**
     *
     */
    public boolean drawCarbons = false;

    /**
     *
     */
    public boolean drawExplicitHydrogens = true;

    /**
     *
     */
    public boolean drawImplicitHydrogens = true;

    /**
     *
     */
    public boolean drawTerminalCarbons = true;

    /**
     *
     */
    public int atomSymbolFontSize = 10;

    /**
     *
     */
    public int plusFontSize = 14;

    /**
     *
     */
    public boolean drawMappings = true;

    /**
     *
     */
    public int subgraphBoxXBorder = 1;

    /**
     *
     */
    public int subgraphBoxYBorder = 2;

    /**
     *
     */
    public double doubleBondGap = 2;

    /**
     *
     */
    public int subscriptHeight = 2;

    /**
     *
     */
    public int subscriptTextSize = 9;

    /**
     *
     */
    public boolean drawAromaticCircles = true;

    /**
     *
     */
    public double ringProportion = 0.75;

    /**
     *
     */
    public float bondStrokeWidth = 1.1f;

    /**
     *
     */
    public double offsetBondDistanceProportion = 0.75;

    /**
     *
     */
    public int filledWedgeWidth = 6;

    /**
     *
     */
    public double wiggleLineWidth = 4;

    /**
     *
     */
    public boolean drawAtomID = false;

    /**
     *
     */
    public int atomIDFontSize = 7;

    /**
     *
     */
    public double labelYGap = 10;

    /**
     *
     */
    public int moleculeLabelFontSize = 7;

    /**
     *
     */
    public int leftToRightMoleculeLabelFontSize = 9;

    /**
     *
     */
    public int topToBottomMoleculeLabelFontSize = 8;

    /**
     *
     */
    public boolean drawMoleculeID = true;

    /**
     *
     */
    public boolean drawLonePairs = true;

    /**
     *
     */
    public double electronRadius = 1.0;

    /**
     *
     */
    public double bondMarkLength = 6;

    /**
     *
     */
    public boolean drawSubgraphBoxes = true;

    /**
     *
     */
    public double doubleMarkGap = 1;

    /**
     *
     */
    public int lonePairSeparation = 4;

    /**
     *
     */
    public boolean drawHighlights = true;

    /**
     *
     */
    public double highlightRadius = 8;

    /**
     *
     */
    public Color highlightColor = BLUE;

    /**
     *
     */
    public boolean highlightsAbove = true;

    /**
     *
     */
    public boolean highlightsBelow = false;

    /**
     *
     */
    public float highlightAlpha = 0.15f;

    /**
     *
     */
    public float highlightBondStroke = 4.0f;

    /**
     *
     */
    public boolean drawSubgraphMappingLines = false;

    /**
     *
     */
    public boolean colorSubgraphBoxes = true;

    /**
     *
     */
    public boolean drawReactionID = false;

    /**
     *
     */
    public boolean layoutLeftToRight = true;

    /**
     *
     */
    public boolean highlightSubgraphs = false;

    /**
     *
     */
    public boolean drawBondStereoChanges = true;

    /**
     *
     */
    public double arrowHeadAngle = 45;

    /**
     *
     */
    public double circularHighlightBorder = 5;

    /**
     *
     */
    public boolean useCircularHighlight = false;

    /**
     *
     */
    public double circularHighlightMinRadius = 10;

    /**
     *
     */
    public boolean circularHighlightIsConcentric = true;

    /**
     *
     */
    public boolean circularHighlightTransparentFilled = false;

    /**
     *
     */
    public boolean useAntialias = true;

    /**
     *
     */
    public double tripleBondGap = 2.5;

    /**
     *
     */
    public boolean drawRS = false;

    /**
     *
     */
    public int chiralSymbolFontSize = 9;

    /**
     *
     */
    public float dashedWedgeStroke = 1.0f;

    /**
     *
     */
    public double dashedGapFactor = 0.1;

    /**
     *
     */
    public double dashedWidthFactor = 0.2;

    /**
     *
     */
    public double dashedWedgeWidth = 6;

    /**
     *
     */
    public int arrowHeadIndent = 5;

    /**
     *
     */
    public int arrowBodyWidth = 5;

    /**
     *
     */
    public boolean drawFatArrow = false;

    /**
     *
     */
    public boolean drawArrowFilled = false;

    /**
     *
     */
    public ArrowType arrowType = FORWARD;

    /**
     *
     */
    public boolean alignMolecules = false;

    /**
     *
     */
    public MoleculeAlignMethod moleculeAlignMethod = MAX_AXIS;

    /**
     *
     */
    public boolean circularHighlightShowAtoms = true;

    /**
     *
     */
    public boolean drawBondFormedCleavedMarks = true;

    /**
     *
     */
    public boolean drawBondOrderChangedMarks = true;

    /**
     *
     */
    public boolean drawLabelPanel = false;

    /**
     *
     */
    public String labelPanelFont = "ROMAN";

    /**
     *
     */
    public int labelPanelFontSize = 14;

    /**
     *
     */
    public boolean shouldCrop = true;

    /**
     *
     */
    public double labelPanelHeight = 20;

    /**
     *
     */
    public double labelGap = 10;;

    /**
     *
     */
    public enum BondStrokeCap {

        /**
         *
         */
        BUTT, 

        /**
         *
         */
        ROUND, 

        /**
         *
         */
        SQUARE
    }

    /**
     *
     */
    public enum BondStrokeJoin {

        /**
         *
         */
        BEVEL, 

        /**
         *
         */
        MITRE, 

        /**
         *
         */
        ROUND
    }

    /**
     *
     */
    public enum XAlign {

        /**
         *
         */
        LEFT, 

        /**
         *
         */
        CENTER, 

        /**
         *
         */
        RIGHT
    }

    /**
     *
     */
    public enum YAlign {

        /**
         *
         */
        TOP, 

        /**
         *
         */
        CENTER, 

        /**
         *
         */
        BOTTOM
    }

    /**
     *
     */
    public enum ArrowType {

        /**
         *
         */
        FORWARD, 

        /**
         *
         */
        BACKWARD, 

        /**
         *
         */
        BIDIRECTIONAL
    }

    /**
     *
     */
    public enum MoleculeAlignMethod {

        /**
         *
         */
        MAX_AXIS, 

        /**
         *
         */
        MIN_AREA
    }
}
