/*
 * ChemicalFileIO - consolidated MDL chemical file readers and writers.
 * Merged from: MDLV2000Reader, MDLV2000Writer, MDLRXNV2000Reader, MDLV2000RXNWriter, MDLValence
 */
package com.bioinceptionlabs.reactionblast.tools;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.DefaultChemObjectWriter;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.formats.MDLFormat;
import org.openscience.cdk.io.formats.MDLRXNFormat;
import org.openscience.cdk.io.formats.MDLV2000Format;
import org.openscience.cdk.io.setting.BooleanIOSetting;
import org.openscience.cdk.io.setting.IOSetting;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupBracket;
import org.openscience.cdk.sgroup.SgroupKey;
import org.openscience.cdk.sgroup.SgroupType;
import org.openscience.cdk.stereo.StereoElementFactory;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import static java.text.NumberFormat.getNumberInstance;
import static java.util.Locale.ENGLISH;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.REMARK;
import static org.openscience.cdk.CDKConstants.TITLE;
import static org.openscience.cdk.geometry.GeometryUtil.has2DCoordinates;
import static org.openscience.cdk.io.formats.MDLFormat.getInstance;
import static org.openscience.cdk.isomorphism.matchers.Expr.Type.ALIPHATIC_ORDER;
import static org.openscience.cdk.isomorphism.matchers.Expr.Type.IS_AROMATIC;
import static org.openscience.cdk.isomorphism.matchers.Expr.Type.OR;
import static org.openscience.cdk.isomorphism.matchers.Expr.Type.ORDER;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Container class for MDL chemical file I/O operations.
 * Contains static inner classes for reading and writing MDL V2000 mol and RXN files.
 */
public final class ChemicalFileIO {


    /**
     * Adds implicit hydrogens and specifies valency using the MDL valence model.
     *
     * @author John May
     * @cdk.module io
     * @see
     * <a href="http://nextmovesoftware.com/blog/2013/02/27/explicit-and-implicit-hydrogens-taking-liberties-with-valence/">Explicit
     * and Implicit Hydrogens: taking liberties with valence</a>
     */
    public static final class MDLValence {

        private MDLValence() {
        }

        /**
         * Apply the MDL valence model to the provided atom container.
         *
         * @param container an atom container loaded from an MDL format
         * @return the container (for convenience)
         */
        static IAtomContainer apply(IAtomContainer container) {

            int n = container.getAtomCount();

            int[] valences = new int[n];

            Map<IAtom, Integer> atomToIndex = Maps.newHashMapWithExpectedSize(n);
            for (IAtom atom : container.atoms()) {
                atomToIndex.put(atom, atomToIndex.size());
            }

            // compute the bond order sums
            for (IBond bond : container.bonds()) {
                int u = atomToIndex.get(bond.getAtom(0));
                int v = atomToIndex.get(bond.getAtom(1));

                int bondOrder = bond.getOrder().numeric();

                valences[u] += bondOrder;
                valences[v] += bondOrder;
            }

            for (int i = 0; i < n; i++) {

                IAtom atom = container.getAtom(i);
                Integer charge = atom.getFormalCharge();
                Integer element = atom.getAtomicNumber();

                if (element == null) {
                    continue;
                }

                // unset = 0 in this case
                charge = charge == null ? 0 : charge;

                int explicit = valences[i];

                // if there was a valence read from the mol file use that otherwise
                // use the default value from the valence model to set the correct
                // number of implied hydrogens
                if (atom.getValency() != null) {
                    atom.setImplicitHydrogenCount(atom.getValency() - explicit);
                } else {
                    int implicit = implicitValence(element, charge, valences[i]);
                    atom.setImplicitHydrogenCount(implicit - explicit);
                    atom.setValency(implicit);
                }
            }

            return container;
        }

        /**
         * Given an element (atomic number) its charge and the explicit valence
         * (bond order sum) return the implicit valence for that atom. This valence
         * is from the MDL valence model which was decoded by NextMove Software and
         * licenced as below.
         *
         * <blockquote> $Id: MDLValence.h 2288 2012-11-26 03:39:27Z glandrum $
         *
         * Copyright (C) 2012 NextMove Software
         *
         * @@ All Rights Reserved @@ This file is part of the RDKit. The contents
         * are covered by the terms of the BSD license which is included in the file
         * license.txt, found at the root of the RDKit source tree. </blockquote>
         * @see
         * <a href="http://nextmovesoftware.com/blog/2013/02/27/explicit-and-implicit-hydrogens-taking-liberties-with-valence/">Explicit
         * and Implicit Hydrogens taking liberties with valence</a>
         */
        static int implicitValence(int elem, int q, int val) {
            switch (elem) {
                case 1: // H
                case 3: // Li
                case 11: // Na
                case 19: // K
                case 37: // Rb
                case 55: // Cs
                case 87: // Fr
                    if (q == 0 && val <= 1) {
                        return 1;
                    }
                    break;

                case 4: // Be
                case 12: // Mg
                case 20: // Ca
                case 38: // Sr
                case 56: // Ba
                case 88: // Ra
                    switch (q) {
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 1:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 5: // B
                    switch (q) {
                        case -4:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                        case -3:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case -2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case -1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 2:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 6: // C
                    switch (q) {
                        case -3:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                        case -2:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case -1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 0:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 2:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 3:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 7: // N
                    switch (q) {
                        case -2:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 3:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 4:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 8: // O
                    switch (q) {
                        case -1:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 2:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 3:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 4:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 5:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 9: // F
                    switch (q) {
                        case 0:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 3:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 4:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 5:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 6:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 13: // Al
                    switch (q) {
                        case -4:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -3:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case -1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 2:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 14: // Si
                    switch (q) {
                        case -3:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -2:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 0:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 2:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 3:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 15: // P
                    switch (q) {
                        case -2:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 3:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 4:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 16: // S
                    switch (q) {
                        case -1:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 2:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 3:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 4:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 5:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 17: // Cl
                    switch (q) {
                        case 0:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 3:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 4:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 5:
                            if (val <= 2) {
                                return 2;
                            }
                            break;
                        case 6:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 31: // Ga
                    switch (q) {
                        case -4:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -3:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case -1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 2:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 32: // Ge
                    switch (q) {
                        case -3:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -2:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 0:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 3:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 33: // As
                    switch (q) {
                        case -2:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 1:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 4:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 34: // Se
                    switch (q) {
                        case -1:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 2:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 3:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 5:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 35: // Br
                    switch (q) {
                        case 0:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 3:
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 4:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 6:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 49: // In
                    switch (q) {
                        case -4:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -3:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 2:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 50: // Sn
                case 82: // Pb
                    switch (q) {
                        case -3:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -2:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 3:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 51: // Sb
                case 83: // Bi
                    switch (q) {
                        case -2:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 0:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 4:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 52: // Te
                case 84: // Po
                    switch (q) {
                        case -1:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 0:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 1:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 2:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 3:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 5:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 53: // I
                case 85: // At
                    switch (q) {
                        case 0:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case 1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case 2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case 3:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 4:
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                        case 6:
                            if (val <= 1) {
                                return 1;
                            }
                            break;
                    }
                    break;

                case 81: // Tl
                    switch (q) {
                        case -4:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            if (val <= 7) {
                                return 7;
                            }
                            break;
                        case -3:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            if (val <= 6) {
                                return 6;
                            }
                            break;
                        case -2:
                            if (val <= 3) {
                                return 3;
                            }
                            if (val <= 5) {
                                return 5;
                            }
                            break;
                        case -1:
                            if (val <= 2) {
                                return 2;
                            }
                            if (val <= 4) {
                                return 4;
                            }
                            break;
                        case 0:
                            if (val <= 1) {
                                return 1;
                            }
                            if (val <= 3) {
                                return 3;
                            }
                            break;
                    }
                    break;

            }
            return val;
        }
    }



    /**
     * RDT format added
     *
     * Writes MDL molfiles, which contains a single molecule (see {
     *
     * @cdk.cite DAL92}). For writing a MDL molfile you can this code:  <pre>
     * MDLV2000Writer writer = new MDLV2000Writer(
     *   new FileWriter(new File("output.mol"))
     * );
     * writer.write((IAtomContainer)molecule);
     * writer.close();
     * </pre>
     *
     * <p>
     * The writer has two IO settings: one for writing 2D coordinates, even if 3D
     * coordinates are given for the written data; the second writes aromatic bonds
     * as bond type 4, which is, strictly speaking, a query bond type, but my many
     * tools used to reflect aromaticity. The full IO setting API is explained in
     * CDK News {
     *
     * @cdk.cite WILLIGHAGEN2004}. One programmatic option to set the option for
     * writing 2D coordinates looks like:  <pre>
     * Properties customSettings = new Properties();
     * customSettings.setProperty(
     *  "ForceWriteAs2DCoordinates", "true"
     * );
     * PropertiesListener listener =
     *   new PropertiesListener(customSettings);
     * writer.addChemObjectIOListener(listener);
     * </pre>
     *
     * @cdk.module io
     * @cdk.githash
     * @cdk.iooptions
     * @cdk.keyword file format, MDL molfile
     */
    public static class MDLV2000Writer extends DefaultChemObjectWriter {

        private static final String NEW_LINE = System.lineSeparator();
        public static final String OptForceWriteAs2DCoordinates = "ForceWriteAs2DCoordinates";
        public static final String OptWriteMajorIsotopes = "WriteMajorIsotopes";
        public static final String OptWriteAromaticBondTypes = "WriteAromaticBondTypes";
        public static final String OptWriteQueryFormatValencies = "WriteQueryFormatValencies";
        public static final String OptWriteDefaultProperties = "WriteDefaultProperties";

        private final static ILoggingTool logger = LoggingToolFactory.createLoggingTool(MDLV2000Writer.class);

        // regular expression to capture R groups with attached numbers
        private Pattern NUMERED_R_GROUP = Pattern.compile("R(\\d+)");

        /**
         * Enumeration of all valid radical values.
         */
        public enum SPIN_MULTIPLICITY {

            None(0, 0),
            Monovalent(2, 1),
            DivalentSinglet(1, 2),
            DivalentTriplet(3, 2);

            // the radical SDF value
            private final int value;
            // the corresponding number of single electrons
            private final int singleElectrons;

            private SPIN_MULTIPLICITY(int value, int singleElectrons) {
                this.value = value;
                this.singleElectrons = singleElectrons;
            }

            /**
             * Radical value for the spin multiplicity in the properties block.
             *
             * @return the radical value
             */
            public int getValue() {
                return value;
            }

            /**
             * The number of single electrons that correspond to the spin
             * multiplicity.
             *
             * @return the number of single electrons
             */
            public int getSingleElectrons() {
                return singleElectrons;
            }

            /**
             * Create a SPIN_MULTIPLICITY instance for the specified value.
             *
             * @param value input value (in the property block)
             * @return instance
             * @throws CDKException unknown spin multiplicity value
             */
            public static SPIN_MULTIPLICITY ofValue(int value) throws CDKException {
                switch (value) {
                    case 0:
                        return None;
                    case 1:
                        return DivalentSinglet;
                    case 2:
                        return Monovalent;
                    case 3:
                        return DivalentTriplet;
                    default:
                        throw new CDKException("unknown spin multiplicity: " + value);
                }
            }
        }

        // number of entries on line; value = 1 to 8
        private static final int NN8 = 8;
        // spacing between entries on line
        private static final int WIDTH = 3;

        private BooleanIOSetting forceWriteAs2DCoords;

        private BooleanIOSetting writeMajorIsotopes;

        // The next two options are MDL Query format options, not really
        // belonging to the MDLV2000 format, and will be removed when
        // a MDLV2000QueryWriter is written.

        /*
         * Should aromatic bonds be written as bond type 4? If true, this makes the
         * output a query file.
         */
        private BooleanIOSetting writeAromaticBondTypes;

        /* Should atomic valencies be written in the Query format. */
        @Deprecated
        private BooleanIOSetting writeQueryFormatValencies;

        private BooleanIOSetting writeDefaultProps;

        private BufferedWriter writer;

        /**
         * Constructs a new MDLWriter that can write an {@link IAtomContainer} to
         * the MDL molfile format.
         *
         * @param out The Writer to write to
         */
        public MDLV2000Writer(Writer out) {
            if (out instanceof BufferedWriter) {
                writer = (BufferedWriter) out;
            } else {
                writer = new BufferedWriter(out);
            }
            initIOSettings();
        }

        /**
         * Constructs a new MDLWriter that can write an {@link IAtomContainer} to a
         * given OutputStream.
         *
         * @param output The OutputStream to write to
         */
        public MDLV2000Writer(OutputStream output) {
            this(new OutputStreamWriter(output, StandardCharsets.UTF_8));
        }

        public MDLV2000Writer() {
            this(new StringWriter());
        }

        @Override
        public IResourceFormat getFormat() {
            return MDLFormat.getInstance();
        }

        @Override
        public void setWriter(Writer out) throws CDKException {
            if (out instanceof BufferedWriter) {
                writer = (BufferedWriter) out;
            } else {
                writer = new BufferedWriter(out);
            }
        }

        @Override
        public void setWriter(OutputStream output) throws CDKException {
            setWriter(new OutputStreamWriter(output));
        }

        /**
         * Flushes the output and closes this object.
         */
        @Override
        public void close() throws IOException {
            writer.close();
        }

        @Override
        public boolean accepts(Class<? extends IChemObject> classObject) {
            Class<?>[] interfaces = classObject.getInterfaces();
            for (int i = 0; i < interfaces.length; i++) {
                if (IAtomContainer.class.equals(interfaces[i])) {
                    return true;
                }
                if (IChemFile.class.equals(interfaces[i])) {
                    return true;
                }
                if (IChemModel.class.equals(interfaces[i])) {
                    return true;
                }
            }
            if (IAtomContainer.class.equals(classObject)) {
                return true;
            }
            if (IChemFile.class.equals(classObject)) {
                return true;
            }
            if (IChemModel.class.equals(classObject)) {
                return true;
            }
            Class superClass = classObject.getSuperclass();
            if (superClass != null) {
                return this.accepts(superClass);
            }
            return false;
        }

        /**
         * Writes a {@link IChemObject} to the MDL molfile formated output. It can
         * only output ChemObjects of type {@link IChemFile},
         * {@link IChemObject} and {@link IAtomContainer}.
         *
         * @param object {@link IChemObject} to write
         * @see #accepts(Class)
         */
        @Override
        public void write(IChemObject object) throws CDKException {
            customizeJob();
            try {
                if (object instanceof IChemFile) {
                    writeChemFile((IChemFile) object);
                    return;
                } else if (object instanceof IChemModel) {
                    IChemFile file = object.getBuilder().newInstance(IChemFile.class);
                    IChemSequence sequence = object.getBuilder().newInstance(IChemSequence.class);
                    sequence.addChemModel((IChemModel) object);
                    file.addChemSequence(sequence);
                    writeChemFile((IChemFile) file);
                    return;
                } else if (object instanceof IAtomContainer) {
                    writeMolecule((IAtomContainer) object);
                    return;
                }
            } catch (Exception ex) {
                logger.error(ex.getMessage());
                logger.debug(ex);
                throw new CDKException("Exception while writing MDL file: " + ex.getMessage(), ex);
            }
            throw new CDKException("Only supported is writing of IChemFile, " + "IChemModel, and IAtomContainer objects.");
        }

        private void writeChemFile(IChemFile file) throws Exception {
            IAtomContainer bigPile = file.getBuilder().newInstance(IAtomContainer.class);
            for (IAtomContainer container : ChemFileManipulator.getAllAtomContainers(file)) {
                bigPile.add(container);
                if (container.getTitle() != null) {
                    if (bigPile.getTitle() != null) {
                        bigPile.setTitle(bigPile.getTitle() + "; " + container.getTitle());
                    } else {
                        bigPile.setTitle(container.getTitle());
                    }
                }
                if (container.getProperty(CDKConstants.REMARK) != null) {
                    if (bigPile.getProperty(CDKConstants.REMARK) != null) {
                        bigPile.setProperty(CDKConstants.REMARK, bigPile.getProperty(CDKConstants.REMARK) + "; "
                                + container.getProperty(CDKConstants.REMARK));
                    } else {
                        bigPile.setProperty(CDKConstants.REMARK, container.getProperty(CDKConstants.REMARK));
                    }
                }
            }
            writeMolecule(bigPile);
        }

        /**
         * Writes a Molecule to an OutputStream in MDL sdf format.
         *
         * @param container Molecule that is written to an OutputStream
         */
        public void writeMolecule(IAtomContainer container) throws Exception {

            /*
             Check for 2D co-ordinates for EC-BLAST
             */
            if (!has2DCoordinates(container)) {
                try {
                    /*
                     Clone it else it will loose mol ID
                     */
                    IAtomContainer clone = container.clone();
                    StructureDiagramGenerator sdg = new StructureDiagramGenerator(clone);
                    sdg.generateCoordinates();
                    container = sdg.getMolecule();
                } catch (CDKException e) {
                }
            }

            final int dim = getNumberOfDimensions(container);
            StringBuilder line = new StringBuilder();
            Map<Integer, Integer> rgroups = null;
            Map<Integer, String> aliases = null;

            /*
             * Add molecule ID for EC-BLAST
             */
            if (container != null
                    && container.getTitle() == null
                    && container.getID() != null) {
                container.setProperty(CDKConstants.TITLE, container.getID());
                container.setTitle(container.getID());
            }

            // write header block
            // lines get shortened to 80 chars, that's in the spec
            String title = container.getTitle();
            if (title == null) {
                title = "";
            }
            if (title.length() > 80) {
                title = title.substring(0, 80);
            }
            writer.write(title);
            writer.newLine();

            /*
             * From CTX spec This line has the format:
             * IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR (FORTRAN:
             * A2<--A8--><---A10-->A2I2<--F10.5-><---F12.5--><-I6-> ) User's first
             * and last initials (l), program name (P), date/time (M/D/Y,H:m),
             * dimensional codes (d), scaling factors (S, s), energy (E) if modeling
             * program input, internal registry number (R) if input through MDL
             * form. A blank line can be substituted for line 2.
             */
            //Overwitten for EC-BLAST
            writer.write("  RDT     ");
            writer.write(new SimpleDateFormat("MMddyyHHmm").format(System.currentTimeMillis()));
            if (dim != 0) {
                writer.write(Integer.toString(dim));
                writer.write('D');
            }
            writer.newLine();

            String comment = (String) container.getProperty(CDKConstants.REMARK);
            if (comment == null) {
                comment = "";
            }
            if (comment.length() > 80) {
                comment = comment.substring(0, 80);
            }
            writer.write(comment);
            writer.newLine();

            // index stereo elements for setting atom parity values
            Map<IAtom, ITetrahedralChirality> atomstereo = new HashMap<>();
            Map<IAtom, Integer> atomindex = new HashMap<>();
            for (IStereoElement element : container.stereoElements()) {
                if (element instanceof ITetrahedralChirality) {
                    atomstereo.put(((ITetrahedralChirality) element).getChiralAtom(), (ITetrahedralChirality) element);
                }
            }
            for (IAtom atom : container.atoms()) {
                atomindex.put(atom, atomindex.size());
            }

            // write Counts line
            line.append(formatMDLInt(container.getAtomCount(), 3));
            line.append(formatMDLInt(container.getBondCount(), 3));
            line.append("  0  0");
            // we mark all stereochemistry to absolute for now
            line.append(atomstereo.isEmpty() ? "  0" : "  1");
            line.append("  0  0  0  0  0999 V2000");
            writer.write(line.toString());
            writer.newLine();

            // write Atom block
            for (int f = 0; f < container.getAtomCount(); f++) {
                IAtom atom = container.getAtom(f);
                line.setLength(0);
                switch (dim) {
                    case 0:
                        // if no coordinates available, then output a number
                        // of zeros
                        line.append("    0.0000    0.0000    0.0000 ");
                        break;
                    case 2:
                        if (atom.getPoint2d() != null) {
                            line.append(formatMDLFloat((float) atom.getPoint2d().x));
                            line.append(formatMDLFloat((float) atom.getPoint2d().y));
                            line.append("    0.0000 ");
                        } else {
                            line.append("    0.0000    0.0000    0.0000 ");
                        }
                        break;
                    case 3:
                        if (atom.getPoint3d() != null) {
                            line.append(formatMDLFloat((float) atom.getPoint3d().x));
                            line.append(formatMDLFloat((float) atom.getPoint3d().y));
                            line.append(formatMDLFloat((float) atom.getPoint3d().z)).append(" ");
                        } else {
                            line.append("    0.0000    0.0000    0.0000 ");
                        }
                        break;
                }
                if (container.getAtom(f) instanceof IPseudoAtom) {
                    //according to http://www.google.co.uk/url?sa=t&ct=res&cd=2&url=http%3A%2F%2Fwww.mdl.com%2Fdownloads%2Fpublic%2Fctfile%2Fctfile.pdf&ei=MsJjSMbjAoyq1gbmj7zCDQ&usg=AFQjCNGaJSvH4wYy4FTXIaQ5f7hjoTdBAw&sig2=eSfruNOSsdMFdlrn7nhdAw an R group is written as R#
                    IPseudoAtom pseudoAtom = (IPseudoAtom) container.getAtom(f);
                    String label = pseudoAtom.getLabel();
                    if (label == null) // set to empty string if null
                    {
                        label = "";
                    }

                    // firstly check if it's a numbered R group
                    Matcher matcher = NUMERED_R_GROUP.matcher(label);
                    if (pseudoAtom.getSymbol().equals("R") && !label.isEmpty() && matcher.matches()) {

                        line.append("R# ");
                        if (rgroups == null) {
                            // we use a tree map to ensure the output order is always the same
                            rgroups = new TreeMap<Integer, Integer>();
                        }
                        rgroups.put(f + 1, Integer.parseInt(matcher.group(1)));

                    } // not a numbered R group - note the symbol may still be R
                    else {

                        // note: no distinction made between alias and pseudo atoms - normally
                        //       aliases maintain their original symbol while pseudo atoms are
                        //       written with a 'A' in the atom block
                        // if the label is longer then 3 characters we need
                        // to use an alias.
                        if (label.length() > 3) {

                            if (aliases == null) {
                                aliases = new TreeMap<Integer, String>();
                            }

                            aliases.put(f + 1, label); // atom index to alias

                            line.append(formatMDLString(atom.getSymbol(), 3));

                        } else { // label is short enough to fit in the atom block

                            // make sure it's not empty
                            if (!label.isEmpty()) {
                                line.append(formatMDLString(label, 3));
                            } else {
                                line.append(formatMDLString(atom.getSymbol(), 3));
                            }

                        }
                    }

                } else {
                    line.append(formatMDLString(container.getAtom(f).getSymbol(), 3));
                }

                // atom properties
                int[] atomprops = new int[12];
                atomprops[0] = determineIsotope(atom);
                atomprops[1] = determineCharge(container, atom);
                atomprops[2] = determineStereoParity(container, atomstereo, atomindex, atom);
                atomprops[5] = determineValence(container, atom);
                atomprops[9] = determineAtomMap(atom);
                line.append(formatMDLInt(atomprops[0], 2)); // dd (mass-number)
                line.append(formatMDLInt(atomprops[1], 3)); // ccc (charge)
                int last = atomprops.length - 1;
                if (!writeDefaultProps.isSet()) {
                    while (last >= 0) {
                        if (atomprops[last] != 0) {
                            break;
                        }
                        last--;
                    }
                    // matches BIOVIA syntax
                    if (last >= 2 && last < atomprops.length) {
                        last = 5;
                    }
                }
                for (int i = 2; i <= last; i++) {
                    line.append(formatMDLInt(atomprops[i], 3));
                }
                line.append(NEW_LINE);
                writer.write(line.toString());
            }

            // write Bond block
            for (IBond bond : container.bonds()) {
                line.setLength(0);
                if (bond.getAtomCount() != 2) {
                    logger.warn("Skipping bond with more/less than two atoms: " + bond);
                } else {
                    if (bond.getStereo() == IBond.Stereo.UP_INVERTED || bond.getStereo() == IBond.Stereo.DOWN_INVERTED
                            || bond.getStereo() == IBond.Stereo.UP_OR_DOWN_INVERTED) {
                        // turn around atom coding to correct for inv stereo
                        line.append(formatMDLInt(atomindex.get(bond.getEnd()) + 1, 3));
                        line.append(formatMDLInt(atomindex.get(bond.getBegin()) + 1, 3));
                    } else {
                        line.append(formatMDLInt(atomindex.get(bond.getBegin()) + 1, 3));
                        line.append(formatMDLInt(atomindex.get(bond.getEnd()) + 1, 3));
                    }

                    int bondType = 0;

                    if (bond instanceof QueryBond) {
                        QueryBond qbond = ((QueryBond) bond);
                        Expr e = qbond.getExpression();
                        switch (e.type()) {
                            case ALIPHATIC_ORDER:
                            case ORDER:
                                bondType = e.value();
                                break;
                            case IS_AROMATIC:
                                bondType = 4;
                                break;
                            case SINGLE_OR_DOUBLE:
                                bondType = 5;
                                break;
                            case SINGLE_OR_AROMATIC:
                                bondType = 6;
                                break;
                            case DOUBLE_OR_AROMATIC:
                                bondType = 7;
                                break;
                            case TRUE:
                                bondType = 8;
                                break;
                            case OR:
                                // SINGLE_OR_DOUBLE
                                if (e.equals(new Expr(Expr.Type.ALIPHATIC_ORDER, 1).or(new Expr(Expr.Type.ALIPHATIC_ORDER, 2)))
                                        || e.equals(new Expr(Expr.Type.ALIPHATIC_ORDER, 2).or(new Expr(Expr.Type.ALIPHATIC_ORDER, 1)))) {
                                    bondType = 5;
                                } // SINGLE_OR_AROMATIC
                                else if (e.equals(new Expr(Expr.Type.ALIPHATIC_ORDER, 1).or(new Expr(Expr.Type.IS_AROMATIC)))
                                        || e.equals(new Expr(Expr.Type.IS_AROMATIC).or(new Expr(Expr.Type.ALIPHATIC_ORDER, 1)))) {
                                    bondType = 6;
                                } // DOUBLE_OR_AROMATIC
                                else if (e.equals(new Expr(Expr.Type.ALIPHATIC_ORDER, 2).or(new Expr(Expr.Type.IS_AROMATIC)))
                                        || e.equals(new Expr(Expr.Type.IS_AROMATIC).or(new Expr(Expr.Type.ALIPHATIC_ORDER, 2)))) {
                                    bondType = 6;
                                }
                                break;
                            default:
                                throw new IllegalArgumentException("Unsupported bond type!");
                        }
                    } else {
                        if (bond.getOrder() != null) {
                            switch (bond.getOrder()) {
                                case SINGLE:
                                case DOUBLE:
                                case TRIPLE:
                                    if (writeAromaticBondTypes.isSet() && bond.isAromatic()) {
                                        bondType = 4;
                                    } else {
                                        bondType = bond.getOrder().numeric();
                                    }
                                    break;
                                case UNSET:
                                    if (bond.isAromatic()) {
                                        if (!writeAromaticBondTypes.isSet()) {
                                            throw new CDKException("Bond at idx " + container.indexOf(bond) + " was an unspecific aromatic bond which should only be used for querie in Molfiles. These can be written if desired by enabling the option 'WriteAromaticBondTypes'.");
                                        }
                                        bondType = 4;
                                    }
                                    break;
                            }
                        }
                    }

                    if (bondType == 0) {
                        throw new CDKException("Bond at idx=" + container.indexOf(bond) + " is not supported by Molfile, bond=" + bond.getOrder());
                    }

                    line.append(formatMDLInt(bondType, 3));
                    line.append("  ");
                    switch (bond.getStereo()) {
                        case UP:
                            line.append("1");
                            break;
                        case UP_INVERTED:
                            line.append("1");
                            break;
                        case DOWN:
                            line.append("6");
                            break;
                        case DOWN_INVERTED:
                            line.append("6");
                            break;
                        case UP_OR_DOWN:
                            line.append("4");
                            break;
                        case UP_OR_DOWN_INVERTED:
                            line.append("4");
                            break;
                        case E_OR_Z:
                            line.append("3");
                            break;
                        default:
                            line.append("0");
                    }
                    if (writeDefaultProps.isSet()) {
                        line.append("  0  0  0 ");
                    }
                    line.append(NEW_LINE);
                    writer.write(line.toString());
                }
            }

            // Write Atom Value
            for (int i = 0; i < container.getAtomCount(); i++) {
                IAtom atom = container.getAtom(i);
                if (atom.getProperty(CDKConstants.COMMENT) != null
                        && atom.getProperty(CDKConstants.COMMENT) instanceof String
                        && !((String) atom.getProperty(CDKConstants.COMMENT)).trim().equals("")) {
                    writer.write("V  ");
                    writer.write(formatMDLInt(i + 1, 3));
                    writer.write(" ");
                    writer.write((String) atom.getProperty(CDKConstants.COMMENT));
                    writer.newLine();
                }
            }

            // write formal atomic charges
            for (int i = 0; i < container.getAtomCount(); i++) {
                IAtom atom = container.getAtom(i);
                Integer charge = atom.getFormalCharge();
                if (charge != null && charge != 0) {
                    writer.write("M  CHG  1 ");
                    writer.write(formatMDLInt(i + 1, 3));
                    writer.write(" ");
                    writer.write(formatMDLInt(charge, 3));
                    writer.newLine();
                }
            }

            // write radical information
            if (container.getSingleElectronCount() > 0) {
                Map<Integer, SPIN_MULTIPLICITY> atomIndexSpinMap = new LinkedHashMap<Integer, SPIN_MULTIPLICITY>();
                for (int i = 0; i < container.getAtomCount(); i++) {
                    int eCount = container.getConnectedSingleElectronsCount(container.getAtom(i));
                    switch (eCount) {
                        case 0:
                            continue;
                        case 1:
                            atomIndexSpinMap.put(i, SPIN_MULTIPLICITY.Monovalent);
                            break;
                        case 2:
                            // information loss, divalent but singlet or triplet?
                            atomIndexSpinMap.put(i, SPIN_MULTIPLICITY.DivalentSinglet);
                            break;
                        default:
                            logger.debug("Invalid number of radicals found: " + eCount);
                            break;
                    }
                }
                Iterator<Map.Entry<Integer, SPIN_MULTIPLICITY>> iterator = atomIndexSpinMap.entrySet().iterator();
                for (int i = 0; i < atomIndexSpinMap.size(); i += NN8) {
                    if (atomIndexSpinMap.size() - i <= NN8) {
                        writer.write("M  RAD" + formatMDLInt(atomIndexSpinMap.size() - i, WIDTH));
                        writeRadicalPattern(iterator, 0);
                    } else {
                        writer.write("M  RAD" + formatMDLInt(NN8, WIDTH));
                        writeRadicalPattern(iterator, 0);
                    }
                    writer.newLine();
                }
            }

            // write formal isotope information
            for (int i = 0; i < container.getAtomCount(); i++) {
                IAtom atom = container.getAtom(i);
                if (!(atom instanceof IPseudoAtom)) {
                    Integer atomicMass = atom.getMassNumber();
                    if (!writeMajorIsotopes.isSet()
                            && isMajorIsotope(atom)) {
                        atomicMass = null;
                    }
                    if (atomicMass != null) {
                        writer.write("M  ISO  1 ");
                        writer.write(formatMDLInt(i + 1, 3));
                        writer.write(" ");
                        writer.write(formatMDLInt(atomicMass, 3));
                        writer.newLine();
                    }
                }
            }

            //write RGP line (max occurrence is 16 data points per line)
            if (rgroups != null) {
                StringBuilder rgpLine = new StringBuilder();
                int cnt = 0;

                // the order isn't guarantied but as we index with the atom
                // number this isn't an issue
                for (Map.Entry<Integer, Integer> e : rgroups.entrySet()) {
                    rgpLine.append(formatMDLInt(e.getKey(), 4));
                    rgpLine.append(formatMDLInt(e.getValue(), 4));
                    cnt++;
                    if (cnt == 8) {
                        rgpLine.insert(0, "M  RGP" + formatMDLInt(cnt, 3));
                        writer.write(rgpLine.toString());
                        writer.newLine();
                        rgpLine = new StringBuilder();
                        cnt = 0;
                    }
                }
                if (cnt != 0) {
                    rgpLine.insert(0, "M  RGP" + formatMDLInt(cnt, 3));
                    writer.write(rgpLine.toString());
                    writer.newLine();
                }

            }

            // write atom aliases
            if (aliases != null) {

                for (Map.Entry<Integer, String> e : aliases.entrySet()) {

                    writer.write("A" + formatMDLInt(e.getKey(), 5));
                    writer.newLine();

                    String label = e.getValue();

                    // fixed width file - doubtful someone would have a label > 70 but trim if they do
                    if (label.length() > 70) {
                        label = label.substring(0, 70);
                    }

                    writer.write(label);
                    writer.newLine();

                }
            }

            writeSgroups(container, writer, atomindex);

            // close molecule
            writer.write("M  END");
            writer.newLine();
            writer.flush();
        }

        // 0 = uncharged or value other than these, 1 = +3, 2 = +2, 3 = +1,
        // 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
        private int determineCharge(IAtomContainer mol, IAtom atom) {
            Integer q = atom.getFormalCharge();
            if (q == null) {
                q = 0;
            }
            switch (q) {
                case -3:
                    return 7;
                case -2:
                    return 6;
                case -1:
                    return 5;
                case 0:
                    if (mol.getConnectedSingleElectronsCount(atom) == 1) {
                        return 4;
                    }
                    return 0;
                case +1:
                    return 3;
                case +2:
                    return 2;
                case +3:
                    return 1;
            }
            return 0;
        }

        private int determineIsotope(IAtom atom) {
            Integer mass = atom.getMassNumber();
            IIsotope major = null;
            if (mass == null) {
                return 0;
            }
            try {
                major = Isotopes.getInstance().getMajorIsotope(atom.getSymbol());
            } catch (IOException e) {
                // ignored
            }
            if (!writeMajorIsotopes.isSet()
                    && major != null
                    && mass.equals(major.getMassNumber())) {
                mass = null;
            }
            if (mass != null) {
                mass -= major != null ? major.getMassNumber() : 0;
                return mass >= -3 && mass <= 4 ? mass : 0;
            }
            return 0;
        }

        private int determineAtomMap(IAtom atom) {
            Object amap = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
            if (amap == null) {
                return 0;
            }
            if (amap instanceof Integer) {
                return (Integer) amap;
            } else {
                if (amap instanceof String) {
                    try {
                        return Integer.parseInt((String) amap);
                    } catch (NumberFormatException ex) {
                        //ignored
                    }
                }
                logger.warn("Skipping non-integer atom map: " + amap
                        + " type:" + amap);
                return 0;
            }
        }

        private int determineValence(IAtomContainer container, IAtom atom) {
            int explicitValence = (int) AtomContainerManipulator.getBondOrderSum(container, atom);
            int charge = atom.getFormalCharge() == null ? 0 : atom.getFormalCharge();
            Integer element = atom.getAtomicNumber();
            int valence = 0;

            if (element != null) {
                int implied = MDLValence.implicitValence(element, charge, explicitValence);
                int actual;
                if (atom.getImplicitHydrogenCount() != null) {
                    actual = explicitValence + atom.getImplicitHydrogenCount();
                } else if (atom.getValency() != null) {
                    actual = atom.getValency();
                } else {
                    return 0;
                }
                if (implied != actual) {
                    if (actual == 0) {
                        return 15;
                    } else if (actual > 0 && actual < 15) {
                        return actual;
                    }
                }
            }
            return valence;
        }

        private int determineStereoParity(IAtomContainer container,
                Map<IAtom, ITetrahedralChirality> atomstereo,
                Map<IAtom, Integer> atomindex, IAtom atom) {
            final ITetrahedralChirality tc = atomstereo.get(atom);
            if (tc == null) {
                return 0;
            }
            int parity = tc.getStereo() == ITetrahedralChirality.Stereo.CLOCKWISE ? 1 : 2;
            IAtom focus = tc.getChiralAtom();
            IAtom[] carriers = tc.getLigands();

            int hidx = -1;
            for (int i = 0; i < 4; i++) {
                // hydrogen position
                if (carriers[i].equals(focus) || carriers[i].getAtomicNumber() == 1) {
                    if (hidx >= 0) {
                        parity = 0;
                    }
                    hidx = i;
                }
            }

            if (parity != 0) {
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        int a = atomindex.get(carriers[i]);
                        int b = atomindex.get(carriers[j]);
                        if (i == hidx) {
                            a = container.getAtomCount();
                        }
                        if (j == hidx) {
                            b = container.getAtomCount();
                        }
                        if (a > b) {
                            parity ^= 0x3;
                        }
                    }
                }
            }
            return parity;
        }

        private boolean isMajorIsotope(IAtom atom) {
            if (atom.getMassNumber() == null) {
                return false;
            }
            try {
                IIsotope major = Isotopes.getInstance().getMajorIsotope(atom.getSymbol());
                return major != null && major.getMassNumber().equals(atom.getMassNumber());
            } catch (IOException ex) {
                return false;
            }
        }

        private void writeSgroups(IAtomContainer container, BufferedWriter writer, Map<IAtom, Integer> atomidxs) throws IOException {
            List<Sgroup> sgroups = container.getProperty(CDKConstants.CTAB_SGROUPS);
            if (sgroups == null) {
                return;
            }

            // going to modify
            sgroups = new ArrayList<>(sgroups);

            // remove non-ctab Sgroups
            Iterator<Sgroup> iter = sgroups.iterator();
            while (iter.hasNext()) {
                if (iter.next().getType() == SgroupType.ExtMulticenter) {
                    iter.remove();
                }
            }

            for (List<Sgroup> wrapSgroups : wrap(sgroups, 8)) {
                // Declare the SGroup type
                writer.write("M  STY");
                writer.write(formatMDLInt(wrapSgroups.size(), 3));
                for (Sgroup sgroup : wrapSgroups) {
                    writer.write(' ');
                    writer.write(formatMDLInt(1 + sgroups.indexOf(sgroup), 3));
                    writer.write(' ');
                    writer.write(sgroup.getType().getKey());
                }
                writer.newLine();
            }

            // Sgroup output is non-compact for now - but valid
            for (int id = 1; id <= sgroups.size(); id++) {
                Sgroup sgroup = sgroups.get(id - 1);

                // Sgroup Atom List
                for (List<IAtom> atoms : wrap(sgroup.getAtoms(), 15)) {
                    writer.write("M  SAL ");
                    writer.write(formatMDLInt(id, 3));
                    writer.write(formatMDLInt(atoms.size(), 3));
                    for (IAtom atom : atoms) {
                        writer.write(' ');
                        writer.write(formatMDLInt(1 + atomidxs.get(atom), 3));
                    }
                    writer.newLine();
                }

                // Sgroup Bond List
                for (List<IBond> bonds : wrap(sgroup.getBonds(), 15)) {
                    writer.write("M  SBL ");
                    writer.write(formatMDLInt(id, 3));
                    writer.write(formatMDLInt(bonds.size(), 3));
                    for (IBond bond : bonds) {
                        writer.write(' ');
                        writer.write(formatMDLInt(1 + container.indexOf(bond), 3));
                    }
                    writer.newLine();
                }

                // Sgroup Parent List
                for (List<Sgroup> parents : wrap(sgroup.getParents(), 8)) {
                    writer.write("M  SPL");
                    writer.write(formatMDLInt(parents.size(), 3));
                    for (Sgroup parent : parents) {
                        writer.write(' ');
                        writer.write(formatMDLInt(id, 3));
                        writer.write(' ');
                        writer.write(formatMDLInt(1 + sgroups.indexOf(parent), 3));
                    }
                    writer.newLine();
                }

                Set<SgroupKey> attributeKeys = sgroup.getAttributeKeys();
                // TODO order and aggregate attribute keys
                for (SgroupKey key : attributeKeys) {
                    switch (key) {
                        case CtabSubScript:
                            writer.write("M  SMT ");
                            writer.write(formatMDLInt(id, 3));
                            writer.write(' ');
                            writer.write((String) sgroup.getValue(key));
                            writer.newLine();
                            break;
                        case CtabExpansion:
                            final boolean expanded = sgroup.getValue(key);
                            if (expanded) {
                                writer.write("M  SDS EXP");
                                writer.write(formatMDLInt(1, 3));
                                writer.write(' ');
                                writer.write(formatMDLInt(id, 3));
                                writer.newLine();
                            }
                            break;
                        case CtabBracket:
                            final List<SgroupBracket> brackets = sgroup.getValue(key);
                            for (SgroupBracket bracket : brackets) {
                                writer.write("M  SDI ");
                                writer.write(formatMDLInt(id, 3));
                                writer.write(formatMDLInt(4, 3));
                                writer.write(formatMDLFloat((float) bracket.getFirstPoint().x));
                                writer.write(formatMDLFloat((float) bracket.getFirstPoint().y));
                                writer.write(formatMDLFloat((float) bracket.getSecondPoint().x));
                                writer.write(formatMDLFloat((float) bracket.getSecondPoint().y));
                                writer.newLine();
                            }
                            break;
                        case CtabBracketStyle:
                            writer.write("M  SBT");
                            writer.write(formatMDLInt(1, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt(id, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt((int) sgroup.getValue(key), 3));
                            writer.newLine();
                            break;
                        case CtabConnectivity:
                            writer.write("M  SCN");
                            writer.write(formatMDLInt(1, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt(id, 3));
                            writer.write(' ');
                            writer.write(((String) sgroup.getValue(key)).toUpperCase(Locale.ROOT));
                            writer.newLine();
                            break;
                        case CtabSubType:
                            writer.write("M  SST");
                            writer.write(formatMDLInt(1, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt(id, 3));
                            writer.write(' ');
                            writer.write((String) sgroup.getValue(key));
                            writer.newLine();
                            break;
                        case CtabParentAtomList:
                            Set<IAtom> parentAtomList = sgroup.getValue(key);
                            for (List<IAtom> atoms : wrap(parentAtomList, 15)) {
                                writer.write("M  SPA ");
                                writer.write(formatMDLInt(id, 3));
                                writer.write(formatMDLInt(atoms.size(), 3));
                                for (IAtom atom : atoms) {
                                    writer.write(' ');
                                    writer.write(formatMDLInt(1 + atomidxs.get(atom), 3));
                                }
                                writer.newLine();
                            }
                            break;
                        case CtabComponentNumber:
                            Integer compNumber = sgroup.getValue(key);
                            writer.write("M  SNC");
                            writer.write(formatMDLInt(1, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt(id, 3));
                            writer.write(' ');
                            writer.write(formatMDLInt(compNumber, 3));
                            writer.newLine();
                            break;
                    }
                }

            }
        }

        private <T> List<List<T>> wrap(Collection<T> set, int lim) {
            List<List<T>> wrapped = new ArrayList<>();
            List<T> list = new ArrayList<T>(set);
            if (set.size() <= lim) {
                if (!list.isEmpty()) {
                    wrapped.add(list);
                }
            } else {
                int i = 0;
                for (; (i + lim) < set.size(); i += lim) {
                    wrapped.add(list.subList(i, i + lim));
                }
                wrapped.add(list.subList(i, list.size()));
            }
            return wrapped;
        }

        private int getNumberOfDimensions(IAtomContainer mol) {
            for (IAtom atom : mol.atoms()) {
                if (atom.getPoint3d() != null && !forceWriteAs2DCoords.isSet()) {
                    return 3;
                } else if (atom.getPoint2d() != null) {
                    return 2;
                }
            }
            return 0;
        }

        private void writeRadicalPattern(Iterator<Map.Entry<Integer, SPIN_MULTIPLICITY>> iterator, int i)
                throws IOException {

            Map.Entry<Integer, SPIN_MULTIPLICITY> entry = iterator.next();
            writer.write(" ");
            writer.write(formatMDLInt(entry.getKey() + 1, WIDTH));
            writer.write(" ");
            writer.write(formatMDLInt(entry.getValue().getValue(), WIDTH));

            i = i + 1;
            if (i < NN8 && iterator.hasNext()) {
                writeRadicalPattern(iterator, i);
            }
        }

        /**
         * Formats an integer to fit into the connection table and changes it to a
         * String.
         *
         * @param x The int to be formated
         * @param n Length of the String
         * @return The String to be written into the connectiontable
         */
        protected static String formatMDLInt(int x, int n) {
            char[] buf = new char[n];
            Arrays.fill(buf, ' ');
            String val = Integer.toString(x);
            if (val.length() > n) {
                val = "0";
            }
            int off = n - val.length();
            for (int i = 0; i < val.length(); i++) {
                buf[off + i] = val.charAt(i);
            }
            return new String(buf);
        }

        /**
         * Formats a float to fit into the connectiontable and changes it to a
         * String.
         *
         * @param fl The float to be formated
         * @return The String to be written into the connectiontable
         */
        protected static String formatMDLFloat(float fl) {
            String s = "", fs = "";
            int l;
            NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
            nf.setMinimumIntegerDigits(1);
            nf.setMaximumIntegerDigits(4);
            nf.setMinimumFractionDigits(4);
            nf.setMaximumFractionDigits(4);
            nf.setGroupingUsed(false);
            if (Double.isNaN(fl) || Double.isInfinite(fl)) {
                s = "0.0000";
            } else {
                s = nf.format(fl);
            }
            l = 10 - s.length();
            for (int f = 0; f < l; f++) {
                fs += " ";
            }
            fs += s;
            return fs;
        }

        /**
         * Formats a String to fit into the connectiontable.
         *
         * @param s The String to be formated
         * @param le The length of the String
         * @return The String to be written in the connectiontable
         */
        protected static String formatMDLString(String s, int le) {
            s = s.trim();
            if (s.length() > le) {
                return s.substring(0, le);
            }
            int l;
            l = le - s.length();
            for (int f = 0; f < l; f++) {
                s += " ";
            }
            return s;
        }

        /**
         * Initializes IO settings.<br>
         * Please note with regards to "writeAromaticBondTypes": bond type values 4
         * through 8 are for SSS queries only, so a 'query file' is created if the
         * container has aromatic bonds and this settings is true.
         */
        private void initIOSettings() {
            forceWriteAs2DCoords = addSetting(new BooleanIOSetting(OptForceWriteAs2DCoordinates, IOSetting.Importance.LOW,
                    "Should coordinates always be written as 2D?", "false"));
            writeMajorIsotopes = addSetting(new BooleanIOSetting(OptWriteMajorIsotopes, IOSetting.Importance.LOW,
                    "Write atomic mass of any non-null atomic mass including major isotopes (e.g. [12]C)", "true"));
            writeAromaticBondTypes = addSetting(new BooleanIOSetting(OptWriteAromaticBondTypes, IOSetting.Importance.LOW,
                    "Should aromatic bonds be written as bond type 4?", "false"));
            writeQueryFormatValencies = addSetting(new BooleanIOSetting(OptWriteQueryFormatValencies,
                    IOSetting.Importance.LOW, "Should valencies be written in the MDL Query format? (deprecated)", "false"));
            writeDefaultProps = addSetting(new BooleanIOSetting(OptWriteDefaultProperties,
                    IOSetting.Importance.LOW,
                    "Write trailing zero's on atom/bond property blocks even if they're not used.",
                    "true"));
        }

        /**
         * Convenience method to set the option for writing aromatic bond types.
         *
         * @param val the value.
         */
        public void setWriteAromaticBondTypes(boolean val) {
            try {
                writeAromaticBondTypes.setSetting(Boolean.toString(val));
            } catch (CDKException e) {
                // ignored can't happen since we are statically typed here
            }
        }

        public void customizeJob() {
            getSettings().forEach((setting) -> {
                fireIOSettingQuestion(setting);
            });
        }

    }



    /**
     * Reads content from MDL molfiles and SD files. It can read a {@link
     * IAtomContainer} or {@link IChemModel} from an MDL molfile, and a {@link
     * IChemFile} from a SD file, with a {@link IChemSequence} of {@link
     * IChemModel}'s, where each IChemModel will contain one {@link IAtomContainer}.
     *
     * <p>
     * From the Atom block it reads atomic coordinates, element types and formal
     * charges. From the Bond block it reads the bonds and the orders. Additionally,
     * it reads 'M CHG', 'G ', 'M RAD' and 'M ISO' lines from the property block.
     *
     * <p>
     * If all z coordinates are 0.0, then the xy coordinates are taken as 2D,
     * otherwise the coordinates are read as 3D.
     *
     * <p>
     * The title of the MOL file is read and can be retrieved with:
     * <pre>
     *   molecule.getProperty(CDKConstants.TITLE);
     * </pre>
     *
     * <p>
     * RGroups which are saved in the MDL molfile as R#, are renamed according to
     * their appearance, e.g. the first R# is named R1. With PseudAtom.getLabel()
     * "R1" is returned (instead of R#). This is introduced due to the SAR table
     * generation procedure of Scitegics PipelinePilot.
     *
     * @author steinbeck
     * @author Egon Willighagen
     * @cdk.module io
     * @cdk.githash
     * @cdk.iooptions
     * @cdk.created 2000-10-02
     * @cdk.keyword file format, MDL molfile
     * @cdk.keyword file format, SDF
     * @cdk.bug 1587283
     */
    public static class MDLV2000Reader extends DefaultChemObjectReader {

        private static final String NEW_LINE = System.lineSeparator();
        BufferedReader input = null;
        private static ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(MDLV2000Reader.class);

        private BooleanIOSetting forceReadAs3DCoords;
        private BooleanIOSetting interpretHydrogenIsotopes;
        private BooleanIOSetting addStereoElements;

        // Pattern to remove trailing space (String.trim() will remove leading space, which we don't want)
        private static final Pattern TRAILING_SPACE = Pattern.compile("\\s+$");

        /**
         * Delimits Structure-Data (SD) Files.
         */
        private static final String RECORD_DELIMITER = "$$$$";

        /**
         * Valid pseudo labels.
         */
        private static final Set<String> PSEUDO_LABELS = ImmutableSet.<String>builder().add("*").add("A").add("Q")
                .add("L").add("LP").add("R") // XXX: not in spec
                .add("R#").build();

        public MDLV2000Reader() {
            this(new StringReader(""));
        }

        /**
         * Constructs a new MDLReader that can read Molecule from a given
         * InputStream.
         *
         * @param in The InputStream to read from
         */
        public MDLV2000Reader(InputStream in) {
            this(new InputStreamReader(in));
        }

        public MDLV2000Reader(InputStream in, Mode mode) {
            this(new InputStreamReader(in), mode);
        }

        /**
         * Constructs a new MDLReader that can read Molecule from a given Reader.
         *
         * @param in The Reader to read from
         */
        public MDLV2000Reader(Reader in) {
            this(in, Mode.RELAXED);
        }

        public MDLV2000Reader(Reader in, Mode mode) {
            input = new BufferedReader(in);
            initIOSettings();
            super.mode = mode;
        }

        @Override
        public IResourceFormat getFormat() {
            return MDLV2000Format.getInstance();
        }

        @Override
        public void setReader(Reader input) throws CDKException {
            if (input instanceof BufferedReader) {
                this.input = (BufferedReader) input;
            } else {
                this.input = new BufferedReader(input);
            }
        }

        @Override
        public void setReader(InputStream input) throws CDKException {
            setReader(new InputStreamReader(input));
        }

        @SuppressWarnings("unchecked")
        @Override
        public boolean accepts(Class<? extends IChemObject> classObject) {
            Class<?>[] interfaces = classObject.getInterfaces();
            for (Class<?> anInterface : interfaces) {
                if (IChemFile.class.equals(anInterface)) {
                    return true;
                }
                if (IChemModel.class.equals(anInterface)) {
                    return true;
                }
                if (IAtomContainer.class.equals(anInterface)) {
                    return true;
                }
            }
            if (IAtomContainer.class.equals(classObject)) {
                return true;
            }
            if (IChemFile.class.equals(classObject)) {
                return true;
            }
            if (IChemModel.class.equals(classObject)) {
                return true;
            }
            Class superClass = classObject.getSuperclass();
            return superClass != null && this.accepts(superClass);
        }

        /**
         * Takes an object which subclasses IChemObject, e.g.Molecule, and will read
         * this (from file, database, internet etc). If the specific implementation
         * does not support a specific IChemObject it will throw an Exception.
         *
         * @param <T>
         * @param object The object that subclasses IChemObject
         * @return The IChemObject read
         * @throws CDKException
         */
        @SuppressWarnings("unchecked")
        @Override
        public <T extends IChemObject> T read(T object) throws CDKException {
            if (object instanceof IAtomContainer) {
                return (T) readAtomContainer((IAtomContainer) object);
            } else if (object instanceof IChemFile) {
                return (T) readChemFile((IChemFile) object);
            } else if (object instanceof IChemModel) {
                return (T) readChemModel((IChemModel) object);
            } else {
                throw new CDKException("Only supported are ChemFile and Molecule.");
            }
        }

        private IChemModel readChemModel(IChemModel chemModel) throws CDKException {
            IAtomContainerSet setOfMolecules = chemModel.getMoleculeSet();
            if (setOfMolecules == null) {
                setOfMolecules = chemModel.getBuilder().newInstance(IAtomContainerSet.class);
            }
            IAtomContainer m = readAtomContainer(chemModel.getBuilder().newInstance(IAtomContainer.class));
            if (m != null) {
                setOfMolecules.addAtomContainer(m);
            }
            chemModel.setMoleculeSet(setOfMolecules);
            return chemModel;
        }

        /**
         * Read a ChemFile from a file in MDL SDF format.
         *
         * @return The ChemFile that was read from the MDL file.
         */
        private IChemFile readChemFile(IChemFile chemFile) throws CDKException {

            IChemObjectBuilder builder = chemFile.getBuilder();
            IChemSequence sequence = builder.newInstance(IChemSequence.class);

            try {
                IAtomContainer m;
                while ((m = readAtomContainer(builder.newInstance(IAtomContainer.class))) != null) {
                    sequence.addChemModel(newModel(m));
                }
            } catch (CDKException e) {
                throw e;
            } catch (IllegalArgumentException exception) {
                String error = "Error while parsing SDF";
                LOGGER.error(error);
                LOGGER.debug(exception);
                throw new CDKException(error, exception);
            }
            try {
                input.close();
            } catch (Exception exc) {
                String error = "Error while closing file: " + exc.getMessage();
                LOGGER.error(error);
                throw new CDKException(error, exc);
            }

            chemFile.addChemSequence(sequence);
            return chemFile;
        }

        /**
         * Create a new chem model for a single {@link IAtomContainer}.
         *
         * @param container the container to create the model for
         * @return a new {@link IChemModel}
         */
        private static IChemModel newModel(final IAtomContainer container) {

            if (container == null) {
                throw new NullPointerException("cannot create chem model for a null container");
            }

            final IChemObjectBuilder builder = container.getBuilder();
            final IChemModel model = builder.newInstance(IChemModel.class);
            final IAtomContainerSet containers = builder.newInstance(IAtomContainerSet.class);

            containers.addAtomContainer(container);
            model.setMoleculeSet(containers);

            return model;
        }

        /**
         * Read an IAtomContainer from a file in MDL sd format
         *
         * @return The Molecule that was read from the MDL file.
         */
        private IAtomContainer readAtomContainer(IAtomContainer molecule) throws CDKException {

            IAtomContainer outputContainer = null;
            Map<IAtom, Integer> parities = new HashMap<>();

            int linecount = 0;
            String title = null;
            String program = null;
            String remark = null;
            String line = "";

            try {

                line = input.readLine();
                linecount++;
                if (line == null) {
                    return null;
                }

                if (line.startsWith("$$$$")) {
                    return molecule;
                }
                if (line.trim().length() > 0) {
                    title = line;
                }
                line = input.readLine();
                linecount++;
                program = line;
                line = input.readLine();
                linecount++;
                if (line.length() > 0) {
                    remark = line;
                }

                line = input.readLine();
                linecount++;

                // if the line is empty we hav a problem - either a malformed
                // molecule entry or just extra new lines at the end of the file
                if (line.length() == 0) {
                    handleError("Unexpected empty line", linecount, 0, 0);
                    // read till the next $$$$ or EOF
                    while (true) {
                        line = input.readLine();
                        linecount++;
                        if (line == null) {
                            return null;
                        }
                        if (line.startsWith("$$$$")) {
                            return molecule; // an empty molecule
                        }
                    }
                }

                final CTabVersion version = CTabVersion.ofHeader(line);

                // check the CT block version
                if (version == CTabVersion.V3000) {
                    handleError("This file must be read with the MDLV3000Reader.");
                    // even if relaxed we can't read V3000 using the V2000 parser
                    throw new CDKException("This file must be read with the MDLV3000Reader.");
                } else if (version == CTabVersion.UNSPECIFIED) {
                    handleError("This file must be read with the MDLReader.");
                    // okay to read in relaxed mode
                }

                int nAtoms = readMolfileInt(line, 0);
                int nBonds = readMolfileInt(line, 3);

                final IAtom[] atoms = new IAtom[nAtoms];
                final IBond[] bonds = new IBond[nBonds];

                // used for applying the MDL valence model
                int[] explicitValence = new int[nAtoms];

                boolean hasX = false, hasY = false, hasZ = false;

                for (int i = 0; i < nAtoms; i++) {
                    line = input.readLine();
                    linecount++;

                    final IAtom atom = readAtomFast(line, molecule.getBuilder(), parities, linecount);

                    atoms[i] = atom;

                    Point3d p = atom.getPoint3d();
                    hasX = hasX || p.x != 0d;
                    hasY = hasY || p.y != 0d;
                    hasZ = hasZ || p.z != 0d;
                }

                // convert to 2D, if totalZ == 0
                if (!hasX && !hasY && !hasZ) {
                    if (nAtoms == 1) {
                        atoms[0].setPoint2d(new Point2d(0, 0));
                    } else {
                        for (IAtom atomToUpdate : atoms) {
                            atomToUpdate.setPoint3d(null);
                        }
                    }
                } else if (!hasZ) {
                    //'  CDK     09251712073D'
                    // 0123456789012345678901
                    if (is3Dfile(program)) {
                        hasZ = true;
                    } else if (!forceReadAs3DCoords.isSet()) {
                        for (IAtom atomToUpdate : atoms) {
                            Point3d p3d = atomToUpdate.getPoint3d();
                            if (p3d != null) {
                                atomToUpdate.setPoint2d(new Point2d(p3d.x, p3d.y));
                                atomToUpdate.setPoint3d(null);
                            }
                        }
                    }
                }

                boolean hasQueryBonds = false;
                for (int i = 0; i < nBonds; i++) {
                    line = input.readLine();
                    linecount++;

                    bonds[i] = readBondFast(line, molecule.getBuilder(), atoms, explicitValence, linecount);
                    hasQueryBonds = hasQueryBonds
                            || (bonds[i].getOrder() == IBond.Order.UNSET && !bonds[i].getFlag(CDKConstants.ISAROMATIC));
                }

                if (!hasQueryBonds) {
                    outputContainer = molecule;
                } else {
                    outputContainer = new QueryAtomContainer(molecule.getBuilder());
                }

                if (title != null) {
                    outputContainer.setTitle(title);
                }
                if (remark != null) {
                    outputContainer.setProperty(CDKConstants.REMARK, remark);
                }

                // if the container is empty we can simply set the atoms/bonds
                // otherwise we add them to the end
                if (outputContainer.isEmpty()) {
                    outputContainer.setAtoms(atoms);
                    outputContainer.setBonds(bonds);
                } else {
                    for (IAtom atom : atoms) {
                        outputContainer.addAtom(atom);
                    }
                    for (IBond bond : bonds) {
                        outputContainer.addBond(bond);
                    }
                }

                // create 0D stereochemistry
                if (addStereoElements.isSet()) {
                    Parities:
                    for (Map.Entry<IAtom, Integer> e : parities.entrySet()) {
                        int parity = e.getValue();
                        if (parity != 1 && parity != 2) {
                            continue; // 3=unspec
                        }
                        int idx = 0;
                        IAtom focus = e.getKey();
                        IAtom[] carriers = new IAtom[4];
                        int hidx = -1;
                        for (IAtom nbr : outputContainer.getConnectedAtomsList(focus)) {
                            if (idx == 4) {
                                continue Parities; // too many neighbors
                            }
                            if (nbr.getAtomicNumber() == 1) {
                                if (hidx >= 0) {
                                    continue Parities;
                                }
                                hidx = idx;
                            }
                            carriers[idx++] = nbr;
                        }
                        // to few neighbors, or already have a hydrogen defined
                        if (idx < 3 || idx < 4 && hidx >= 0) {
                            continue;
                        }
                        if (idx == 3) {
                            carriers[idx++] = focus;
                        }

                        if (idx == 4) {
                            Stereo winding = parity == 1 ? Stereo.CLOCKWISE : Stereo.ANTI_CLOCKWISE;
                            // H is always at back, even if explicit! At least this seems to be the case.
                            // we adjust the winding as needed
                            if (hidx == 0 || hidx == 2) {
                                winding = winding.invert();
                            }
                            outputContainer.addStereoElement(new TetrahedralChirality(focus, carriers, winding));
                        }
                    }
                }

                // read PROPERTY block
                readPropertiesFast(input, outputContainer, nAtoms);

                // read potential SD file data between M  END and $$$$
                readNonStructuralData(input, outputContainer);

                // note: apply the valence model last so that all fixes (i.e. hydrogen
                // isotopes) are in place we need to use a offset as this atoms
                // could be added to a molecule which already had atoms present
                int offset = outputContainer.getAtomCount() - nAtoms;
                for (int i = offset; i < outputContainer.getAtomCount(); i++) {
                    int valence = explicitValence[i - offset];
                    if (valence < 0) {
                        hasQueryBonds = true; // also counts aromatic bond as query
                    } else {
                        int unpaired = outputContainer.getConnectedSingleElectronsCount(outputContainer.getAtom(i));
                        applyMDLValenceModel(outputContainer.getAtom(i), valence + unpaired, unpaired);
                    }
                }

                // sanity check that we have a decent molecule, query bonds mean we
                // don't have a hydrogen count for atoms and stereo perception isn't
                // currently possible
                if (!hasQueryBonds && addStereoElements.isSet() && hasX && hasY) {
                    if (hasZ) { // has 3D coordinates
                        outputContainer.setStereoElements(StereoElementFactory.using3DCoordinates(outputContainer)
                                .createAll());
                    } else if (!forceReadAs3DCoords.isSet()) { // has 2D coordinates (set as 2D coordinates)
                        outputContainer.setStereoElements(StereoElementFactory.using2DCoordinates(outputContainer)
                                .createAll());
                    }
                }

            } catch (CDKException exception) {
                String error = "Error while parsing line " + linecount + ": " + line + " -> " + exception.getMessage();
                LOGGER.error(error);
                throw exception;
            } catch (IOException exception) {
                String error = "Error while parsing line " + linecount + ": " + line + " -> " + exception.getMessage();
                LOGGER.error(error);
                handleError("Error while parsing line: " + line, linecount, 0, 0, exception);
            }

            /*
             * Set TITLE as ID for molecules EC-BLAST
             */
            String property = outputContainer == null ? ""
                    : (String) outputContainer.getTitle();
            if (outputContainer != null && property != null) {
                outputContainer.setID(property);
            }
            return outputContainer;
        }

        private boolean is3Dfile(String program) {
            return program.length() >= 22 && program.substring(20, 22).equals("3D");
        }

        /**
         * Applies the MDL valence model to atoms using the explicit valence (bond
         * order sum) and charge to determine the correct number of implicit
         * hydrogens. The model is not applied if the explicit valence is less than
         * 0 - this is the case when a query bond was read for an atom.
         *
         * @param atom the atom to apply the model to
         * @param unpaired unpaired electron count
         * @param explicitValence the explicit valence (bond order sum)
         */
        private void applyMDLValenceModel(IAtom atom, int explicitValence, int unpaired) {

            if (atom.getValency() != null) {
                if (atom.getValency() >= explicitValence) {
                    atom.setImplicitHydrogenCount(atom.getValency() - (explicitValence - unpaired));
                } else {
                    atom.setImplicitHydrogenCount(0);
                }
            } else {
                Integer element = atom.getAtomicNumber();
                if (element == null) {
                    element = 0;
                }

                Integer charge = atom.getFormalCharge();
                if (charge == null) {
                    charge = 0;
                }

                int implicitValence = MDLValence.implicitValence(element, charge, explicitValence);
                if (implicitValence < explicitValence) {
                    atom.setValency(explicitValence);
                    atom.setImplicitHydrogenCount(0);
                } else {
                    atom.setValency(implicitValence);
                    atom.setImplicitHydrogenCount(implicitValence - explicitValence);
                }
            }
        }

        private void fixHydrogenIsotopes(IAtomContainer molecule, IsotopeFactory isotopeFactory) {
            for (IAtom atom : AtomContainerManipulator.getAtomArray(molecule)) {
                if (atom instanceof IPseudoAtom) {
                    IPseudoAtom pseudo = (IPseudoAtom) atom;
                    if ("D".equals(pseudo.getLabel())) {
                        IAtom newAtom = molecule.getBuilder().newInstance(IAtom.class, atom);
                        newAtom.setSymbol("H");
                        newAtom.setAtomicNumber(1);
                        isotopeFactory.configure(newAtom, isotopeFactory.getIsotope("H", 2));
                        AtomContainerManipulator.replaceAtomByAtom(molecule, atom, newAtom);
                    } else if ("T".equals(pseudo.getLabel())) {
                        IAtom newAtom = molecule.getBuilder().newInstance(IAtom.class, atom);
                        newAtom.setSymbol("H");
                        newAtom.setAtomicNumber(1);
                        isotopeFactory.configure(newAtom, isotopeFactory.getIsotope("H", 3));
                        AtomContainerManipulator.replaceAtomByAtom(molecule, atom, newAtom);
                    }
                }
            }
        }

        @Override
        public void close() throws IOException {
            input.close();
        }

        private void initIOSettings() {
            forceReadAs3DCoords = addSetting(new BooleanIOSetting("ForceReadAs3DCoordinates", IOSetting.Importance.LOW,
                    "Should coordinates always be read as 3D?", "false"));
            interpretHydrogenIsotopes = addSetting(new BooleanIOSetting("InterpretHydrogenIsotopes",
                    IOSetting.Importance.LOW, "Should D and T be interpreted as hydrogen isotopes?", "true"));
            addStereoElements = addSetting(new BooleanIOSetting("AddStereoElements", IOSetting.Importance.LOW,
                    "Detect and create IStereoElements for the input.", "true"));
        }

        public void customizeJob() {
            getSettings().forEach((setting) -> {
                fireIOSettingQuestion(setting);
            });
        }

        private String removeNonDigits(String input) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < input.length(); i++) {
                char inputChar = input.charAt(i);
                if (Character.isDigit(inputChar)) {
                    sb.append(inputChar);
                }
            }
            return sb.toString();
        }

        IAtom readAtomFast(String line, IChemObjectBuilder builder, int lineNum) throws CDKException, IOException {
            return readAtomFast(line, builder, Collections.<IAtom, Integer>emptyMap(), lineNum);
        }

        /**
         * Parse an atom line from the atom block using the format: {@code
         * xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee}
         * where: <ul> <li>x: x coordinate</li> <li>y: y coordinate</li> <li>z: z
         * coordinate</li> <li>a: atom symbol</li> <li>d: mass difference</li>
         * <li>c: charge</li> <li>s: stereo parity</li> <li>h: hydrogen count + 1
         * (not read - query)</li> <li>b: stereo care (not read - query)</li> <li>v:
         * valence</li> <li>H: H0 designator (not read - query)</li> <li>r: not
         * used</li> <li>i: not used</li> <li>m: atom reaction mapping</li> <li>n:
         * inversion/retention flag</li> <li>e: exact change flag</li> </ul>
         *
         * The parsing is strict and does not allow extra columns (i.e. NMR shifts)
         * malformed input.
         *
         * @param line input line
         * @param builder chem object builder to create the atom
         * @param parities map of atom parities for creation 0D stereochemistry
         * @param lineNum the line number - for printing error messages
         * @return a new atom instance
         */
        IAtom readAtomFast(String line, IChemObjectBuilder builder, Map<IAtom, Integer> parities, int lineNum) throws CDKException, IOException {

            // The line may be truncated and it's checked in reverse at the specified
            // lengths:
            //          1         2         3         4         5         6
            // 123456789012345678901234567890123456789012345678901234567890123456789
            //                                  | |  |  |  |  |  |  |  |  |  |  |  |
            // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
            String symbol;
            double x, y, z;
            int massDiff = 0, charge = 0, parity = 0, valence = 0, mapping = 0;

            int length = length(line);
            if (length > 69) // excess data we should check all fields
            {
                length = 69;
            }

            // given the length we jump to the position and parse all fields
            // that could be present (note - fall through switch)
            switch (length) {
                case 69: // eee: exact charge flag [reaction, query]
                case 66: // nnn: inversion / retention [reaction]
                case 63: // mmm: atom-atom mapping [reaction]
                    mapping = readMolfileInt(line, 60);
                case 60: // iii: not used
                case 57: // rrr: not used
                case 54: // HHH: H0 designation [redundant]
                case 51: // vvv: valence
                    valence = readMolfileInt(line, 48);
                case 48: // bbb: stereo care [query]
                case 45: // hhh: hydrogen count + 1 [query]
                case 42: // sss: stereo parity
                    parity = toInt(line.charAt(41));
                case 39: // ccc: charge
                    charge = toCharge(line.charAt(38));
                case 36: // dd: mass difference
                    massDiff = sign(line.charAt(34)) * toInt(line.charAt(35));
                case 34: // x y z and aaa: atom coordinates and symbol
                case 33: // symbol is left aligned
                case 32:
                    x = readMDLCoordinate(line, 0);
                    y = readMDLCoordinate(line, 10);
                    z = readMDLCoordinate(line, 20);
                    symbol = line.substring(31, 34).trim().intern();
                    break;
                default:
                    handleError("invalid line length", lineNum, 0, 0);
                    throw new CDKException("invalid line length, " + length + ": " + line);
            }

            IAtom atom = createAtom(symbol, builder, lineNum);

            atom.setPoint3d(new Point3d(x, y, z));
            atom.setFormalCharge(charge);
            atom.setStereoParity(parity);
            if (parity != 0) {
                parities.put(atom, parity);
            }

            // if there was a mass difference, set the mass number
            if (massDiff != 0 && atom.getAtomicNumber() > 0) {
                IIsotope majorIsotope = Isotopes.getInstance().getMajorIsotope(atom.getAtomicNumber());
                if (majorIsotope == null) {
                    atom.setMassNumber(-1); // checked after M ISO is processed
                } else {
                    atom.setMassNumber(majorIsotope.getMassNumber() + massDiff);
                }
            }

            if (valence > 0 && valence < 16) {
                atom.setValency(valence == 15 ? 0 : valence);
            }

            if (mapping != 0) {
                atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, mapping);
            }

            return atom;
        }

        /**
         * Read a bond from a line in the MDL bond block. The bond block is
         * formatted as follows, {@code 111222tttsssxxxrrrccc}, where:
         * <ul>
         * <li>111: first atom number</li>
         * <li>222: second atom number</li>
         * <li>ttt: bond type</li>
         * <li>xxx: bond stereo</li>
         * <li>rrr: bond topology</li>
         * <li>ccc: reaction center</li>
         * </ul>
         *
         * @param line the input line
         * @param builder builder to create objects with
         * @param atoms atoms read from the atom block
         * @param explicitValence array to fill with explicit valence
         * @param lineNum the input line number
         * @return a new bond
         * @throws CDKException thrown if the input was malformed or didn't make
         * sense
         */
        IBond readBondFast(String line, IChemObjectBuilder builder, IAtom[] atoms, int[] explicitValence, int lineNum)
                throws CDKException {

            // The line may be truncated and it's checked in reverse at the specified
            // lengths. Absolutely required is atom indices, bond type and stereo.
            //          1         2
            // 123456789012345678901
            //            |  |  |  |
            // 111222tttsssxxxrrrccc
            int length = length(line);
            if (length > 21) {
                length = 21;
            }

            int u, v, type, stereo = 0;

            switch (length) {
                case 21: // ccc: reaction centre status
                case 18: // rrr: bond topology
                case 15: // xxx: not used
                case 12: // sss: stereo
                    stereo = readUInt(line, 9, 3);
                case 9: // 111222ttt: atoms, type and stereo
                    u = readMolfileInt(line, 0) - 1;
                    v = readMolfileInt(line, 3) - 1;
                    type = readMolfileInt(line, 6);
                    break;
                default:
                    throw new CDKException("invalid line length: " + length + " " + line);
            }

            IBond bond = builder.newBond();
            bond.setAtoms(new IAtom[]{atoms[u], atoms[v]});

            switch (type) {
                case 1: // single
                    bond.setOrder(IBond.Order.SINGLE);
                    bond.setStereo(toStereo(stereo, type));
                    break;
                case 2: // double
                    bond.setOrder(IBond.Order.DOUBLE);
                    bond.setStereo(toStereo(stereo, type));
                    break;
                case 3: // triple
                    bond.setOrder(IBond.Order.TRIPLE);
                    break;
                case 4: // aromatic
                    bond.setOrder(IBond.Order.UNSET);
                    bond.setFlag(CDKConstants.ISAROMATIC, true);
                    bond.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
                    atoms[u].setFlag(CDKConstants.ISAROMATIC, true);
                    atoms[v].setFlag(CDKConstants.ISAROMATIC, true);
                    break;
                case 5: // single or double
                    bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.SINGLE_OR_DOUBLE);
                    break;
                case 6: // single or aromatic
                    bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.SINGLE_OR_AROMATIC);
                    break;
                case 7: // double or aromatic
                    bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.DOUBLE_OR_AROMATIC);
                    break;
                case 8: // any
                    bond = new QueryBond(bond.getBegin(), bond.getEnd(), Expr.Type.TRUE);
                    break;
                default:
                    throw new CDKException("unrecognised bond type: " + type + ", " + line);
            }

            if (type < 4) {
                explicitValence[u] += type;
                explicitValence[v] += type;
            } else {
                explicitValence[u] = explicitValence[v] = Integer.MIN_VALUE;
            }

            return bond;
        }

        /**
         * Reads the property block from the {@code input} setting the values in the
         * container.
         *
         * @param input input resource
         * @param container the structure with atoms / bonds present
         * @param nAtoms the number of atoms in the atoms block
         * @throws IOException low-level IO error
         */
        void readPropertiesFast(final BufferedReader input, final IAtomContainer container, final int nAtoms)
                throws IOException, CDKException {
            String line;

            // first atom index in this Molfile, the container may have
            // already had atoms present before reading the file
            int offset = container.getAtomCount() - nAtoms;

            Map<Integer, Sgroup> sgroups = new LinkedHashMap<>();

            LINES:
            while ((line = input.readLine()) != null) {

                int index, count, lnOffset;
                Sgroup sgroup;
                int length = line.length();
                final PropertyKey key = PropertyKey.of(line);
                switch (key) {

                    // A  aaa
                    // x...
                    //
                    // atom alias is stored as label on a pseudo atom
                    case ATOM_ALIAS:
                        index = readMolfileInt(line, 3) - 1;
                        final String label = input.readLine();
                        if (label == null) {
                            return;
                        }
                        label(container, offset + index, label);
                        break;

                    // V  aaa v...
                    //
                    // an atom value is stored as comment on an atom
                    case ATOM_VALUE:
                        index = readMolfileInt(line, 3) - 1;
                        final String comment = line.substring(7);
                        container.getAtom(offset + index).setProperty(CDKConstants.COMMENT, comment);
                        break;

                    // G  aaappp
                    // x...
                    //
                    // Abbreviation is required for compatibility with previous versions of MDL ISIS/Desktop which
                    // allowed abbreviations with only one attachment. The attachment is denoted by two atom
                    // numbers, aaa and ppp. All of the atoms on the aaa side of the bond formed by aaa-ppp are
                    // abbreviated. The coordinates of the abbreviation are the coordinates of aaa. The text of the
                    // abbreviation is on the following line (x...). In current versions of ISIS, abbreviations can have any
                    // number of attachments and are written out using the Sgroup appendixes. However, any ISIS
                    // abbreviations that do have one attachment are also written out in the old style, again for
                    // compatibility with older ISIS versions, but this behavior might not be supported in future
                    // versions.
                    case GROUP_ABBREVIATION:
                        // not supported, existing parsing doesn't do what is
                        // mentioned in the specification above
                        // final int    from  = readMolfileInt(line, 3) - 1;
                        // final int    to    = readMolfileInt(line, 6) - 1;
                        final String group = input.readLine();
                        if (group == null) {
                            return;
                        }
                        break;

                    // M  CHGnn8 aaa vvv ...
                    //
                    // vvv: -15 to +15. Default of 0 = uncharged atom. When present, this property supersedes
                    //      all charge and radical values in the atom block, forcing a 0 charge on all atoms not
                    //      listed in an M CHG or M RAD line.
                    case M_CHG:
                        count = readUInt(line, 6, 3);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            index = readMolfileInt(line, st) - 1;
                            int charge = readMolfileInt(line, st + 4);
                            container.getAtom(offset + index).setFormalCharge(charge);
                        }
                        break;

                    // M  ISOnn8 aaa vvv ...
                    //
                    // vvv: Absolute mass of the atom isotope as a positive integer. When present, this property
                    //      supersedes all isotope values in the atom block. Default (no entry) means natural
                    //      abundance. The difference between this absolute mass value and the natural
                    //      abundance value specified in the PTABLE.DAT file must be within the range of -18
                    //      to +12.
                    case M_ISO:
                        count = readUInt(line, 6, 3);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            index = readMolfileInt(line, st) - 1;
                            int mass = readMolfileInt(line, st + 4);
                            if (mass < 0) {
                                handleError("Absolute mass number should be >= 0, " + line);
                            } else {
                                container.getAtom(offset + index).setMassNumber(mass);
                            }
                        }
                        break;

                    // M  RADnn8 aaa vvv ...
                    //
                    // vvv: Default of 0 = no radical, 1 = singlet (:), 2 = doublet ( . or ^), 3 = triplet (^^). When
                    //      present, this property supersedes all charge and radical values in the atom block,
                    //      forcing a 0 (zero) charge and radical on all atoms not listed in an M CHG or
                    //      M RAD line.
                    case M_RAD:
                        count = readUInt(line, 6, 3);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            index = readMolfileInt(line, st) - 1;
                            int value = readMolfileInt(line, st + 4);
                            MDLV2000Writer.SPIN_MULTIPLICITY multiplicity = MDLV2000Writer.SPIN_MULTIPLICITY.ofValue(value);

                            for (int e = 0; e < multiplicity.getSingleElectrons(); e++) {
                                container.addSingleElectron(offset + index);
                            }
                        }
                        break;

                    // M  RGPnn8 aaa rrr ...
                    //
                    // rrr: Rgroup number, value from 1 to 32 *, labels position of Rgroup on root.
                    //
                    // see also, RGroupQueryReader
                    case M_RGP:
                        count = readUInt(line, 6, 3);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            index = readMolfileInt(line, st) - 1;
                            int number = readMolfileInt(line, st + 4);
                            label(container, offset + index, "R" + number);
                        }
                        break;

                    // M  ZZC aaa c...
                    // 
                    // c: first character of the label, extends to EOL.
                    //
                    // Proprietary atom labels created by ACD/Labs ChemSketch using the Manual Numbering Tool.
                    // This atom property appears to be undocumented, but experimentation leads to the following
                    // specification (tested with ACD/ChemSketch version 12.00 Build 29305, 25 Nov 2008)
                    //
                    // It's not necessary to label any/all atoms but if a label is present, the following applies:
                    //
                    // The atom label(s) consist of an optional prefix, a required numeric label, and optional suffix.
                    //                         
                    // The numeric label is an integer in the range 0 - 999 inclusive.
                    // 
                    // If present, the prefix and suffix can each contain 1 - 50 characters, from the set of printable 
                    // ASCII characters shown here
                    //                            
                    //    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
                    //                    
                    // In addition, both the prefix and suffix may contain leading and/or trailing and/or embedded 
                    // whitespace, included within the limit of 50 characters. These should be preserved when read.
                    //                    
                    // Long labels in the mol/sdfile are not truncated or wrapped onto multiple lines. As a result, the
                    // line could be 114 characters in length (excluding the newline).
                    //
                    // By stopping and restarting the Manual Numbering Tool, it's possible to create non-sequential
                    // or even duplicate numbers or labels. This is reasonable for the intended purpose of the tool - 
                    // labelling the structure as you wish. If unique labels are required, downstream processing will be
                    // necessary to enforce this.
                    //
                    case M_ZZC:
                        if (mode == Mode.STRICT) {
                            throw new CDKException("Atom property ZZC is illegal in STRICT mode");
                        }
                        index = readMolfileInt(line, 7) - 1;
                        String atomLabel = line.substring(11);  // DO NOT TRIM
                        container.getAtom(offset + index).setProperty(CDKConstants.ACDLABS_LABEL, atomLabel);
                        break;

                    // M STYnn8 sss ttt ...
                    //  sss: Sgroup number
                    //  ttt: Sgroup type: SUP = abbreviation Sgroup (formerly called superatom), MUL = multiple group,
                    //                    SRU = SRU type, MON = monomer, MER = Mer type, COP = copolymer, CRO = crosslink,
                    //                    MOD = modification, GRA = graft, COM = component, MIX = mixture,
                    //                    FOR = formulation, DAT = data Sgroup, ANY = any polymer, GEN = generic.
                    //
                    // Note: For a given Sgroup, an STY line giving its type must appear before any other line that
                    //       supplies information about it. For a data Sgroup, an SDT line must describe the data
                    //       field before the SCD and SED lines that contain the data (see Data Sgroup Data below).
                    //       When a data Sgroup is linked to another Sgroup, the Sgroup must already have been defined.
                    //
                    // Sgroups can be in any order on the Sgroup Type line. Brackets are drawn around Sgroups with the
                    // M SDI lines defining the coordinates.
                    case M_STY:
                        count = readMolfileInt(line, 6);
                        for (int i = 0; i < count; i++) {
                            lnOffset = 10 + (i * 8);
                            index = readMolfileInt(line, lnOffset);

                            if (mode == Mode.STRICT && sgroups.containsKey(index)) {
                                handleError("STY line must appear before any other line that supplies Sgroup information");
                            }

                            sgroup = new Sgroup();
                            sgroups.put(index, sgroup);

                            SgroupType type = SgroupType.parseCtabKey(line.substring(lnOffset + 4, lnOffset + 7));
                            if (type != null) {
                                sgroup.setType(type);
                            }
                        }
                        break;

                    // Sgroup Subtype [Sgroup]
                    // M  SSTnn8 sss ttt ...
                    // ttt: Polymer Sgroup subtypes: ALT = alternating, RAN = random, BLO = block
                    case M_SST:
                        count = readMolfileInt(line, 6);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            sgroup = ensureSgroup(sgroups,
                                    readMolfileInt(line, st));
                            if (mode == Mode.STRICT && sgroup.getType() != SgroupType.CtabCopolymer) {
                                handleError("SST (Sgroup Subtype) specified for a non co-polymer group");
                            }

                            String sst = line.substring(st + 4, st + 7);

                            if (mode == Mode.STRICT && !("ALT".equals(sst) || "RAN".equals(sst) || "BLO".equals(sst))) {
                                handleError("Invalid sgroup subtype: " + sst + " expected (ALT, RAN, or BLO)");
                            }

                            sgroup.putValue(SgroupKey.CtabSubType, sst);
                        }
                        break;

                    // Sgroup Atom List [Sgroup]
                    // M   SAL sssn15 aaa ...
                    // aaa: Atoms in Sgroup sss
                    case M_SAL:
                        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
                        count = readMolfileInt(line, 10);
                        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
                            index = readMolfileInt(line, st) - 1;
                            sgroup.addAtom(container.getAtom(offset + index));
                        }
                        break;

                    // Sgroup Bond List [Sgroup]
                    // M  SBL sssn15 bbb ...
                    // bbb: Bonds in Sgroup sss.
                    // (For data Sgroups, bbb’s are the containment bonds, for all other
                    //  Sgroup types, bbb’s are crossing bonds.)
                    case M_SBL:
                        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
                        count = readMolfileInt(line, 10);
                        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
                            index = readMolfileInt(line, st) - 1;
                            sgroup.addBond(container.getBond(offset + index));
                        }
                        break;

                    // Sgroup Hierarchy Information [Sgroup]
                    // M  SPLnn8 ccc ppp ...
                    //   ccc: Sgroup index of the child Sgroup
                    //   ppp: Sgroup index of the parent Sgroup (ccc and ppp must already be defined via an
                    //        STY line prior to encountering this line)
                    case M_SPL:
                        count = readMolfileInt(line, 6);
                        for (int i = 0, st = 10; i < count && st + 6 <= length; i++, st += 8) {
                            sgroup = ensureSgroup(sgroups, readMolfileInt(line, st));
                            sgroup.addParent(ensureSgroup(sgroups, readMolfileInt(line, st + 4)));
                        }
                        break;

                    // Sgroup Connectivity [Sgroup]
                    // M  SCNnn8 sss ttt ...
                    // ttt: HH = head-to-head, HT = head-to-tail, EU = either unknown.
                    // Left justified.
                    case M_SCN:
                        count = readMolfileInt(line, 6);
                        for (int i = 0, st = 10; i < count && st + 6 <= length; i++, st += 8) {
                            sgroup = ensureSgroup(sgroups,
                                    readMolfileInt(line, st));
                            String con = line.substring(st + 4, Math.min(length, st + 7)).trim();
                            if (mode == Mode.STRICT && !("HH".equals(con) || "HT".equals(con) || "EU".equals(con))) {
                                handleError("Unknown SCN type (expected: HH, HT, or EU) was " + con);
                            }
                            sgroup.putValue(SgroupKey.CtabConnectivity,
                                    con);
                        }
                        break;

                    // Sgroup Display Information
                    // M SDI sssnn4 x1 y1 x2 y2
                    // x1,y1, Coordinates of bracket endpoints
                    // x2,y2:
                    case M_SDI:
                        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
                        count = readMolfileInt(line, 10);
                        assert count == 4; // fixed?
                        sgroup.addBracket(new SgroupBracket(readMDLCoordinate(line, 13),
                                readMDLCoordinate(line, 23),
                                readMDLCoordinate(line, 33),
                                readMDLCoordinate(line, 43)));
                        break;

                    // Sgroup subscript
                    // M SMT sss m...
                    // m...: Text of subscript Sgroup sss.
                    // (For multiple groups, m... is the text representation of the multiple group multiplier.
                    //  For abbreviation Sgroups, m... is the text of the abbreviation Sgroup label.)
                    case M_SMT:
                        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
                        sgroup.putValue(SgroupKey.CtabSubScript,
                                line.substring(11).trim());
                        break;

                    // Sgroup Bracket Style
                    // The format for the Sgroup bracket style is as follows:
                    // M  SBTnn8 sss ttt ...
                    // where:
                    //   sss: Index of Sgroup
                    //   ttt: Bracket display style: 0 = default, 1 = curved (parenthetic) brackets
                    // This appendix supports altering the display style of the Sgroup brackets.
                    case M_SBT:
                        count = readMolfileInt(line, 6);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            sgroup = ensureSgroup(sgroups,
                                    readMolfileInt(line, st));
                            sgroup.putValue(SgroupKey.CtabBracketStyle,
                                    readMolfileInt(line, st + 4));
                        }
                        break;

                    // Sgroup Expansion
                    // M  SDS EXPn15 sss ...
                    // sss: Sgroup index of expanded abbreviation Sgroups
                    case M_SDS:

                        if ("EXP".equals(line.substring(7, 10))) {
                            count = readMolfileInt(line, 10);
                            for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
                                sgroup = ensureSgroup(sgroups, readMolfileInt(line, st));
                                sgroup.putValue(SgroupKey.CtabExpansion, true);
                            }
                        } else if (mode == Mode.STRICT) {
                            handleError("Expected EXP to follow SDS tag");
                        }
                        break;

                    // Multiple Group Parent Atom List [Sgroup]
                    // M SPA sssn15 aaa ...
                    // aaa: Atoms in paradigmatic repeating unit of multiple group sss
                    // Note: To ensure that all current molfile readers consistently
                    //       interpret chemical structures, multiple groups are written
                    //       in their fully expanded state to the molfile. The M SPA atom
                    //       list is a subset of the full atom list that is defined by the
                    //       Sgroup Atom List M SAL entry.
                    case M_SPA:
                        sgroup = ensureSgroup(sgroups, readMolfileInt(line, 7));
                        count = readMolfileInt(line, 10);
                        Set<IAtom> parentAtomList = sgroup.getValue(SgroupKey.CtabParentAtomList);
                        if (parentAtomList == null) {
                            sgroup.putValue(SgroupKey.CtabParentAtomList, parentAtomList = new HashSet<IAtom>());
                        }
                        for (int i = 0, st = 14; i < count && st + 3 <= length; i++, st += 4) {
                            index = readMolfileInt(line, st) - 1;
                            parentAtomList.add(container.getAtom(offset + index));
                        }
                        break;

                    // Sgroup Component Numbers [Sgroup]
                    // M  SNCnn8 sss ooo ...
                    // sss: Index of component Sgroup
                    // ooo: Integer component order (1...256). This limit applies only to MACCS-II
                    case M_SNC:
                        count = readMolfileInt(line, 6);
                        for (int i = 0, st = 10; i < count && st + 7 <= length; i++, st += 8) {
                            sgroup = ensureSgroup(sgroups,
                                    readMolfileInt(line, st));
                            sgroup.putValue(SgroupKey.CtabComponentNumber,
                                    readMolfileInt(line, st + 4));
                        }
                        break;

                    // M  END
                    //
                    // This entry goes at the end of the properties block and is required for molfiles which contain a
                    // version stamp in the counts line.
                    case M_END:
                        break LINES;
                }
            }

            // check of ill specified atomic mass
            for (IAtom atom : container.atoms()) {
                if (atom.getMassNumber() != null && atom.getMassNumber() < 0) {
                    handleError("Unstable use of mass delta on " + atom.getSymbol() + " please use M  ISO");
                    atom.setMassNumber(null);
                }
            }

            if (!sgroups.isEmpty()) {
                // load Sgroups into molecule, first we downcast
                List<Sgroup> sgroupOrgList = new ArrayList<>(sgroups.values());
                List<Sgroup> sgroupCpyList = new ArrayList<>(sgroupOrgList.size());
                for (int i = 0; i < sgroupOrgList.size(); i++) {
                    Sgroup cpy = sgroupOrgList.get(i).downcast();
                    sgroupCpyList.add(cpy);
                }
                // update replaced parents
                for (int i = 0; i < sgroupOrgList.size(); i++) {
                    Sgroup newSgroup = sgroupCpyList.get(i);
                    Set<Sgroup> oldParents = new HashSet<>(newSgroup.getParents());
                    newSgroup.removeParents(oldParents);
                    for (Sgroup parent : oldParents) {
                        newSgroup.addParent(sgroupCpyList.get(sgroupOrgList.indexOf(parent)));
                    }
                }
                container.setProperty(CDKConstants.CTAB_SGROUPS, sgroupCpyList);
            }
        }

        private Sgroup ensureSgroup(Map<Integer, Sgroup> map, int idx) throws CDKException {
            Sgroup sgroup = map.get(idx);
            if (sgroup == null) {
                if (mode == Mode.STRICT) {
                    handleError("Sgroups must first be defined by a STY property");
                }
                map.put(idx, sgroup = new Sgroup());
            }
            return sgroup;
        }

        /**
         * Convert an MDL V2000 stereo value to the CDK {@link IBond.Stereo}. The
         * method should only be invoked for single/double bonds. If strict mode is
         * enabled irrational bond stereo/types cause errors (e.g. up double bond).
         *
         * @param stereo stereo value
         * @param type bond type
         * @return bond stereo
         * @throws CDKException the stereo value was invalid (strict mode).
         */
        private IBond.Stereo toStereo(final int stereo, final int type) throws CDKException {
            switch (stereo) {
                case 0:
                    return type == 2 ? IBond.Stereo.E_Z_BY_COORDINATES : IBond.Stereo.NONE;
                case 1:
                    if (mode == Mode.STRICT && type == 2) {
                        throw new CDKException("stereo flag was 'up' but bond order was 2");
                    }
                    return IBond.Stereo.UP;
                case 3:
                    if (mode == Mode.STRICT && type == 1) {
                        throw new CDKException("stereo flag was 'cis/trans' but bond order was 1");
                    }
                    return IBond.Stereo.E_OR_Z;
                case 4:
                    if (mode == Mode.STRICT && type == 2) {
                        throw new CDKException("stereo flag was 'up/down' but bond order was 2");
                    }
                    return IBond.Stereo.UP_OR_DOWN;
                case 6:
                    if (mode == Mode.STRICT && type == 2) {
                        throw new CDKException("stereo flag was 'down' but bond order was 2");
                    }
                    return IBond.Stereo.DOWN;
            }
            if (mode == Mode.STRICT) {
                throw new CDKException("unknown bond stereo type: " + stereo);
            }
            return IBond.Stereo.NONE;
        }

        /**
         * Determine the length of the line excluding trailing whitespace.
         *
         * @param str a string
         * @return the length when trailing white space is removed
         */
        static int length(final String str) {
            int i = str.length() - 1;
            while (i >= 0 && str.charAt(i) == ' ') {
                i--;
            }
            return i + 1;
        }

        /**
         * Create an atom for the provided symbol. If the atom symbol is a periodic
         * element a new 'Atom' is created otherwise if the symbol is an allowed
         * query atom ('R', 'Q', 'A', '*', 'L', 'LP') a new 'PseudoAtom' is created.
         * If the symbol is invalid an exception is thrown.
         *
         * @param symbol input symbol
         * @param builder chem object builder
         * @return a new atom
         * @throws CDKException the symbol is not allowed
         */
        private IAtom createAtom(String symbol, IChemObjectBuilder builder, int lineNum) throws CDKException {
            final Elements elem = Elements.ofString(symbol);
            if (elem != Elements.Unknown) {
                IAtom atom = builder.newAtom();
                atom.setSymbol(elem.symbol());
                atom.setAtomicNumber(elem.number());
                return atom;
            }
            if (symbol.equals("D") && interpretHydrogenIsotopes.isSet()) {
                if (mode == Mode.STRICT) {
                    throw new CDKException("invalid symbol: " + symbol);
                }
                IAtom atom = builder.newInstance(IAtom.class, "H");
                atom.setMassNumber(2);
                return atom;
            }
            if (symbol.equals("T") && interpretHydrogenIsotopes.isSet()) {
                if (mode == Mode.STRICT) {
                    throw new CDKException("invalid symbol: " + symbol);
                }
                IAtom atom = builder.newInstance(IAtom.class, "H");
                atom.setMassNumber(3);
                return atom;
            }

            if (!isPseudoElement(symbol)) {
                handleError("invalid symbol: " + symbol, lineNum, 31, 34);
                // when strict only accept labels from the specification
                if (mode == Mode.STRICT) {
                    throw new CDKException("invalid symbol: " + symbol);
                }
            }

            // will be renumbered later by RGP if R1, R2 etc. if not renumbered then
            // 'R' is a better label than 'R#' if now RGP is specified
            if (symbol.equals("R#")) {
                symbol = "R";
            }

            IAtom atom = builder.newInstance(IPseudoAtom.class, symbol);
            atom.setSymbol(symbol);
            atom.setAtomicNumber(0); // avoid NPE downstream

            return atom;
        }

        /**
         * Is the atom symbol a non-periodic element (i.e. pseudo). Valid pseudo
         * atoms are 'R#', 'A', 'Q', '*', 'L' and 'LP'. We also accept 'R' but this
         * is not listed in the specification.
         *
         * @param symbol a symbol from the input
         * @return the symbol is a valid pseudo element
         */
        static boolean isPseudoElement(final String symbol) {
            return PSEUDO_LABELS.contains(symbol);
        }

        /**
         * Read a coordinate from an MDL input. The MDL V2000 input coordinate has
         * 10 characters, 4 significant figures and is prefixed with whitespace for
         * padding: 'xxxxx.xxxx'. Knowing the format allows us to use an optimised
         * parser which does not consider exponents etc.
         *
         * @param line input line
         * @param offset first character of the coordinate
         * @return the specified value
         * @throws CDKException the coordinates specification was not valid
         */
        double readMDLCoordinate(final String line, int offset) throws CDKException {
            // to be valid the decimal should be at the fifth index (4 sig fig)
            if (line.charAt(offset + 5) != '.') {
                handleError("Bad coordinate format specified, expected 4 decimal places: " + line.substring(offset));
                int start = offset;
                while (line.charAt(start) == ' ' && start < offset + 9) {
                    start++;
                }

                int dot = -1;
                int end = start;
                for (char c = line.charAt(end); c != ' ' && end < offset + 9; c = line.charAt(end), end++) {
                    if (c == '.') {
                        dot = end;
                    }
                }

                if (start == end) {

                    return 0.0;
                } else if (dot != -1) {

                    int sign = sign(line.charAt(start));
                    if (sign < 0) {
                        start++;
                    }

                    int integral = readUInt(line, start, dot - start - 1);
                    int fraction = readUInt(line, dot, end - dot);

                    return sign * (integral * 10000L + fraction) / 10000d;
                } else {

                    return Double.parseDouble(line.substring(start, end));
                }
            } else {
                int start = offset;
                while (line.charAt(start) == ' ') {
                    start++;
                }
                int sign = sign(line.charAt(start));
                if (sign < 0) {
                    start++;
                }
                int integral = readUInt(line, start, (offset + 5) - start);
                int fraction = readUInt(line, offset + 6, 4);
                return sign * (integral * 10000L + fraction) / 10000d;
            }
        }

        /**
         * Convert the a character (from an MDL V2000 input) to a charge value: 1 =
         * +1, 2 = +2, 3 = +3, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3.
         *
         * @param c a character
         * @return formal charge
         */
        private static int toCharge(final char c) {
            switch (c) {
                case '1':
                    return +3;
                case '2':
                    return +2;
                case '3':
                    return +1;
                case '4':
                    return 0; // doublet radical - superseded by M  RAD
                case '5':
                    return -1;
                case '6':
                    return -2;
                case '7':
                    return -3;
            }
            return 0;
        }

        /**
         * Obtain the sign of the character, -1 if the character is '-', +1
         * otherwise.
         *
         * @param c a character
         * @return the sign
         */
        private static int sign(final char c) {
            return c == '-' ? -1 : +1;
        }

        /**
         * Convert a character (ASCII code points) to an integer. If the character
         * was not a digit (i.e. space) the value defaults to 0.
         *
         * @param c a character
         * @return the numerical value
         */
        private static int toInt(final char c) {
            // Character.getNumericalValue allows all of unicode which we don't want
            // or need it - imagine an MDL file with roman numerals!
            return c >= '0' && c <= '9' ? c - '0' : 0;
        }

        /**
         * Read an unsigned int value from the given index with the expected number
         * of digits.
         *
         * @param line input line
         * @param index start index
         * @param digits number of digits (max)
         * @return an unsigned int
         */
        private static int readUInt(final String line, int index, int digits) {
            int result = 0;
            while (digits-- > 0) {
                result = (result * 10) + toInt(line.charAt(index++));
            }
            return result;
        }

        /**
         * Optimised method for reading a integer from 3 characters in a string at a
         * specified index. MDL V2000 Molfile make heavy use of the 3 character ints
         * in the atom/bond and property blocks. The integer may be signed and
         * pre/post padded with white space.
         *
         * @param line input
         * @param index start index
         * @return the value specified in the string
         */
        private static int readMolfileInt(final String line, final int index) {
            int sign = 1;
            int result = 0;
            char c;
            switch ((c = line.charAt(index))) {
                case ' ':
                    break;
                case '-':
                    sign = -1;
                    break;
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    result = (c - '0');
                    break;
                default:
                    return 0;
            }
            switch ((c = line.charAt(index + 1))) {
                case ' ':
                    if (result > 0) {
                        return sign * result;
                    }
                    break;
                case '-':
                    if (result > 0) {
                        return sign * result;
                    }
                    sign = -1;
                    break;
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    result = (result * 10) + (c - '0');
                    break;
                default:
                    return sign * result;
            }
            switch ((c = line.charAt(index + 2))) {
                case ' ':
                    if (result > 0) {
                        return sign * result;
                    }
                    break;
                case '-':
                    if (result > 0) {
                        return sign * result;
                    }
                    sign = -1;
                    break;
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    result = (result * 10) + (c - '0');
                    break;
                default:
                    return sign * result;
            }
            return sign * result;
        }

        /**
         * Labels the atom at the specified index with the provide label. If the
         * atom was not already a pseudo atom then the original atom is replaced.
         *
         * @param container structure
         * @param index atom index to replace
         * @param label the label for the atom
         * @see IPseudoAtom#setLabel(String)
         */
        static void label(final IAtomContainer container, final int index, final String label) {
            final IAtom atom = container.getAtom(index);
            final IPseudoAtom pseudoAtom = atom instanceof IPseudoAtom ? (IPseudoAtom) atom : container.getBuilder()
                    .newInstance(IPseudoAtom.class);
            if (atom.equals(pseudoAtom)) {
                pseudoAtom.setLabel(label);
            } else {
                pseudoAtom.setSymbol(label);
                pseudoAtom.setAtomicNumber(0);
                pseudoAtom.setPoint2d(atom.getPoint2d());
                pseudoAtom.setPoint3d(atom.getPoint3d());
                pseudoAtom.setMassNumber(atom.getMassNumber());
                pseudoAtom.setFormalCharge(atom.getFormalCharge());
                pseudoAtom.setValency(atom.getValency());
                pseudoAtom.setLabel(label);
                // XXX: would be faster to track all replacements and do it all in one
                AtomContainerManipulator.replaceAtomByAtom(container, atom, pseudoAtom);
            }
        }

        /**
         * Reads an atom from the input allowing for non-standard formatting (i.e
         * truncated lines) and chemical shifts.
         *
         * @param line input line
         * @param builder chem object builder
         * @param linecount the current line count
         * @return an atom to add to a container
         * @throws CDKException a CDK error occurred
         * @throws IOException the isotopes file could not be read
         */
        private IAtom readAtomSlow(String line, IChemObjectBuilder builder, int linecount) throws CDKException, IOException {
            IAtom atom;
            Matcher trailingSpaceMatcher = TRAILING_SPACE.matcher(line);
            if (trailingSpaceMatcher.find()) {
                handleError("Trailing space found", linecount, trailingSpaceMatcher.start(), trailingSpaceMatcher.end());
                line = trailingSpaceMatcher.replaceAll("");
            }
            double x = Double.parseDouble(line.substring(0, 10).trim());
            double y = Double.parseDouble(line.substring(10, 20).trim());
            double z = Double.parseDouble(line.substring(20, 30).trim());

            String element = line.substring(31, Math.min(line.length(), 34)).trim();
            if (line.length() < 34) {
                handleError("Element atom type does not follow V2000 format type should of length three"
                        + " and padded with space if required", linecount, 31, 34);
            }

            LOGGER.debug("Atom type: ", element);
            IsotopeFactory isotopeFactory = Isotopes.getInstance();
            if (isotopeFactory.isElement(element)) {
                atom = isotopeFactory.configure(builder.newInstance(IAtom.class, element));
            } else if ("A".equals(element)) {
                atom = builder.newInstance(IPseudoAtom.class, element);
            } else if ("Q".equals(element)) {
                atom = builder.newInstance(IPseudoAtom.class, element);
            } else if ("*".equals(element)) {
                atom = builder.newInstance(IPseudoAtom.class, element);
            } else if ("LP".equals(element)) {
                atom = builder.newInstance(IPseudoAtom.class, element);
            } else if ("L".equals(element)) {
                atom = builder.newInstance(IPseudoAtom.class, element);
            } else if (element.equals("R") || (element.length() > 0 && element.charAt(0) == 'R')) {
                LOGGER.debug("Atom ", element, " is not an regular element. Creating a PseudoAtom.");
                //check if the element is R
                String[] rGroup = element.split("^R");
                if (rGroup.length > 1) {
                    try {
                        element = "R" + Integer.valueOf(rGroup[(rGroup.length - 1)]);
                        atom = builder.newInstance(IPseudoAtom.class, element);

                    } catch (Exception ex) {
                        // This happens for atoms labeled "R#".
                        // The Rnumber may be set later on, using RGP line
                        atom = builder.newInstance(IPseudoAtom.class, "R");
                    }
                } else {
                    atom = builder.newInstance(IPseudoAtom.class, element);
                }
            } else {
                handleError("Invalid element type. Must be an existing " + "element, or one in: A, Q, L, LP, *.",
                        linecount, 32, 35);
                atom = builder.newInstance(IPseudoAtom.class, element);
                atom.setSymbol(element);
            }

            // store as 3D for now, convert to 2D (if totalZ == 0.0) later
            atom.setPoint3d(new Point3d(x, y, z));

            // parse further fields
            if (line.length() >= 36) {
                String massDiffString = line.substring(34, 36).trim();
                LOGGER.debug("Mass difference: ", massDiffString);
                if (!(atom instanceof IPseudoAtom)) {
                    try {
                        int massDiff = Integer.parseInt(massDiffString);
                        if (massDiff != 0) {
                            IIsotope major = Isotopes.getInstance().getMajorIsotope(element);
                            atom.setMassNumber(major.getMassNumber() + massDiff);
                        }
                    } catch (NumberFormatException | IOException exception) {
                        handleError("Could not parse mass difference field.", linecount, 35, 37, exception);
                    }
                } else {
                    LOGGER.error("Cannot set mass difference for a non-element!");
                }
            } else {
                handleError("Mass difference is missing", linecount, 34, 36);
            }

            // set the stereo partiy
            Integer parity = line.length() > 41 ? Character.digit(line.charAt(41), 10) : 0;
            atom.setStereoParity(parity);

            if (line.length() >= 51) {
                String valenceString = removeNonDigits(line.substring(48, 51));
                LOGGER.debug("Valence: ", valenceString);
                if (!(atom instanceof IPseudoAtom)) {
                    try {
                        int valence = Integer.parseInt(valenceString);
                        if (valence != 0) {
                            //15 is defined as 0 in mol files
                            if (valence == 15) {
                                atom.setValency(0);
                            } else {
                                atom.setValency(valence);
                            }
                        }
                    } catch (Exception exception) {
                        handleError("Could not parse valence information field", linecount, 49, 52, exception);
                    }
                } else {
                    LOGGER.error("Cannot set valence information for a non-element!");
                }
            }

            if (line.length() >= 39) {
                String chargeCodeString = line.substring(36, 39).trim();
                LOGGER.debug("Atom charge code: ", chargeCodeString);
                int chargeCode = Integer.parseInt(chargeCodeString);
                if (chargeCode == 0) {
                    // uncharged species
                } else if (chargeCode == 1) {
                    atom.setFormalCharge(+3);
                } else if (chargeCode == 2) {
                    atom.setFormalCharge(+2);
                } else if (chargeCode == 3) {
                    atom.setFormalCharge(+1);
                } else if (chargeCode == 4) {
                } else if (chargeCode == 5) {
                    atom.setFormalCharge(-1);
                } else if (chargeCode == 6) {
                    atom.setFormalCharge(-2);
                } else if (chargeCode == 7) {
                    atom.setFormalCharge(-3);
                }
            } else {
                handleError("Atom charge is missing", linecount, 36, 39);
            }

            try {
                // read the mmm field as position 61-63
                String reactionAtomIDString = line.substring(60, 63).trim();
                LOGGER.debug("Parsing mapping id: ", reactionAtomIDString);
                try {
                    int reactionAtomID = Integer.parseInt(reactionAtomIDString);
                    if (reactionAtomID != 0) {
                        atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, reactionAtomID);
                    }
                } catch (Exception exception) {
                    LOGGER.error("Mapping number ", reactionAtomIDString, " is not an integer.");
                    LOGGER.debug(exception);
                }
            } catch (Exception exception) {
                // older mol files don't have all these fields...
                LOGGER.warn("A few fields are missing. Older MDL MOL file?");
            }

            //shk3: This reads shifts from after the molecule. I don't think this is an official format, but I saw it frequently 80=>78 for alk
            if (line.length() >= 78) {
                double shift = Double.parseDouble(line.substring(69, 80).trim());
                atom.setProperty("first shift", shift);
            }
            if (line.length() >= 87) {
                double shift = Double.parseDouble(line.substring(79, 87).trim());
                atom.setProperty("second shift", shift);
            }

            return atom;
        }

        /**
         * Read a bond line from an MDL V2000 molfile bond block (slow). The
         * explicit valence is also modified.
         *
         * @param line the input from the bond block
         * @param builder chem object builder
         * @param atoms array of atoms
         * @param explicitValence stores the explicit valence of each atom (bond
         * order sum)
         * @param linecount the current line count
         * @return a new bond
         * @throws CDKException the bond line could not be parsed
         */
        private IBond readBondSlow(String line, IChemObjectBuilder builder, IAtom[] atoms, int[] explicitValence,
                int linecount) throws CDKException {
            int atom1 = Integer.parseInt(line.substring(0, 3).trim());
            int atom2 = Integer.parseInt(line.substring(3, 6).trim());
            int order = Integer.parseInt(line.substring(6, 9).trim());
            IBond.Stereo stereo = null;
            if (line.length() >= 12) {
                int mdlStereo = line.length() > 12 ? Integer.parseInt(line.substring(9, 12).trim()) : Integer.parseInt(line
                        .substring(9).trim());
                if (mdlStereo == 1) {
                    // MDL up bond
                    stereo = IBond.Stereo.UP;
                } else if (mdlStereo == 6) {
                    // MDL down bond
                    stereo = IBond.Stereo.DOWN;
                } else if (mdlStereo == 0) {
                    if (order == 2) {
                        // double bond stereo defined by coordinates
                        stereo = IBond.Stereo.E_Z_BY_COORDINATES;
                    } else {
                        // bond has no stereochemistry
                        stereo = IBond.Stereo.NONE;
                    }
                } else if (mdlStereo == 3 && order == 2) {
                    // unknown E/Z stereochemistry
                    stereo = IBond.Stereo.E_OR_Z;
                } else if (mdlStereo == 4) {
                    //MDL bond undefined
                    stereo = IBond.Stereo.UP_OR_DOWN;
                }
            } else {
                handleError("Missing expected stereo field at line: ", linecount, 10, 12);
            }
            if (LOGGER.isDebugEnabled()) {
                LOGGER.debug("Bond: " + atom1 + " - " + atom2 + "; order " + order);
            }
            // interpret CTfile's special bond orders
            IAtom a1 = atoms[atom1 - 1];
            IAtom a2 = atoms[atom2 - 1];
            IBond newBond;
            if (order >= 1 && order <= 3) {
                IBond.Order cdkOrder = IBond.Order.SINGLE;
                if (order == 2) {
                    cdkOrder = IBond.Order.DOUBLE;
                }
                if (order == 3) {
                    cdkOrder = IBond.Order.TRIPLE;
                }
                if (stereo != null) {
                    newBond = builder.newInstance(IBond.class, a1, a2, cdkOrder, stereo);
                } else {
                    newBond = builder.newInstance(IBond.class, a1, a2, cdkOrder);
                }
                explicitValence[atom1 - 1] += cdkOrder.numeric();
                explicitValence[atom2 - 1] += cdkOrder.numeric();
            } else if (order == 4) {
                // aromatic bond
                if (stereo != null) {
                    newBond = builder.newInstance(IBond.class, a1, a2, IBond.Order.UNSET, stereo);
                } else {
                    newBond = builder.newInstance(IBond.class, a1, a2, IBond.Order.UNSET);
                }
                // mark both atoms and the bond as aromatic and raise the SINGLE_OR_DOUBLE-flag
                newBond.setFlag(CDKConstants.SINGLE_OR_DOUBLE, true);
                newBond.setFlag(CDKConstants.ISAROMATIC, true);
                newBond.setIsAromatic(true);
                a1.setFlag(CDKConstants.ISAROMATIC, true);
                a2.setFlag(CDKConstants.ISAROMATIC, true);
                explicitValence[atom1 - 1] = explicitValence[atom2 - 1] = Integer.MIN_VALUE;
            } else {
                newBond = new QueryBond(builder);
                IAtom[] bondAtoms = {a1, a2};
                newBond.setAtoms(bondAtoms);
                switch (order) {
                    case 5:
                        ((QueryBond) newBond).getExpression()
                                .setPrimitive(Expr.Type.SINGLE_OR_DOUBLE);
                        break;
                    case 6:
                        ((QueryBond) newBond).getExpression()
                                .setPrimitive(Expr.Type.SINGLE_OR_AROMATIC);
                        break;
                    case 7:
                        ((QueryBond) newBond).getExpression()
                                .setPrimitive(Expr.Type.DOUBLE_OR_AROMATIC);
                        break;
                    case 8:
                        ((QueryBond) newBond).getExpression()
                                .setPrimitive(Expr.Type.TRUE);
                        break;
                }
                newBond.setStereo(stereo);
                explicitValence[atom1 - 1] = explicitValence[atom2 - 1] = Integer.MIN_VALUE;
            }
            return newBond;
        }

        /**
         * Read the properties from the V2000 block (slow).
         *
         * @param input input source
         * @param container the container with the atoms / bonds loaded
         * @param nAtoms the number of atoms in the atom block
         * @param linecount the line count
         * @throws IOException internal low-level error
         * @throws CDKException the properties block could not be parsed
         */
        private void readPropertiesSlow(BufferedReader input, IAtomContainer container, int nAtoms, int linecount)
                throws IOException, CDKException {
            LOGGER.info("Reading property block");
            String line;
            while (true) {
                line = input.readLine();
                linecount++;
                if (line == null) {
                    handleError("The expected property block is missing!", linecount, 0, 0);
                }
                if (line.startsWith("M  END")) {
                    break;
                }

                boolean lineRead = false;
                if (line.startsWith("M  CHG")) {
                    // FIXME: if this is encountered for the first time, all
                    // atom charges should be set to zero first!
                    int infoCount = Integer.parseInt(line.substring(6, 9).trim());
                    StringTokenizer st = new StringTokenizer(line.substring(9));
                    for (int i = 1; i <= infoCount; i++) {
                        String token = st.nextToken();
                        int atomNumber = Integer.parseInt(token.trim());
                        token = st.nextToken();
                        int charge = Integer.parseInt(token.trim());
                        container.getAtom(atomNumber - 1).setFormalCharge(charge);
                    }
                } else if (line.matches("A\\s{1,4}\\d+")) {
                    // Reads the pseudo atom property from the mol file

                    // The atom number of the to replaced atom
                    int aliasAtomNumber = Integer.parseInt(line.replaceFirst("A\\s{1,4}", ""));
                    String alias = input.readLine();
                    linecount++;
                    IAtom aliasAtom = container.getAtom(aliasAtomNumber - 1);

                    // skip if already a pseudoatom
                    if (aliasAtom instanceof IPseudoAtom) {
                        ((IPseudoAtom) aliasAtom).setLabel(alias);
                        continue;
                    }

                    IAtom newPseudoAtom = container.getBuilder().newInstance(IPseudoAtom.class, alias);
                    if (aliasAtom.getPoint2d() != null) {
                        newPseudoAtom.setPoint2d(aliasAtom.getPoint2d());
                    }
                    if (aliasAtom.getPoint3d() != null) {
                        newPseudoAtom.setPoint3d(aliasAtom.getPoint3d());
                    }
                    AtomContainerManipulator.replaceAtomByAtom(container, aliasAtom, newPseudoAtom);
                } else if (line.startsWith("M  ISO")) {
                    try {
                        String countString = line.substring(6, 10).trim();
                        int infoCount = Integer.parseInt(countString);
                        StringTokenizer st = new StringTokenizer(line.substring(10));
                        for (int i = 1; i <= infoCount; i++) {
                            int atomNumber = Integer.parseInt(st.nextToken().trim());
                            int absMass = Integer.parseInt(st.nextToken().trim());
                            if (absMass != 0) {
                                IAtom isotope = container.getAtom(atomNumber - 1);
                                isotope.setMassNumber(absMass);
                            }
                        }
                    } catch (NumberFormatException exception) {
                        String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": "
                                + line + " in property block.";
                        LOGGER.error(error);
                        handleError("NumberFormatException in isotope information.", linecount, 7, 11, exception);
                    }
                } else if (line.startsWith("M  RAD")) {
                    try {
                        String countString = line.substring(6, 9).trim();
                        int infoCount = Integer.parseInt(countString);
                        StringTokenizer st = new StringTokenizer(line.substring(9));
                        for (int i = 1; i <= infoCount; i++) {
                            int atomNumber = Integer.parseInt(st.nextToken().trim());
                            int rad = Integer.parseInt(st.nextToken().trim());
                            MDLV2000Writer.SPIN_MULTIPLICITY spin = MDLV2000Writer.SPIN_MULTIPLICITY.None;
                            if (rad > 0) {
                                IAtom radical = container.getAtom(atomNumber - 1);
                                spin = MDLV2000Writer.SPIN_MULTIPLICITY.ofValue(rad);
                                for (int j = 0; j < spin.getSingleElectrons(); j++) {
                                    container.addSingleElectron(container.getBuilder().newInstance(ISingleElectron.class,
                                            radical));
                                }
                            }
                        }
                    } catch (NumberFormatException exception) {
                        String error = "Error (" + exception.getMessage() + ") while parsing line " + linecount + ": "
                                + line + " in property block.";
                        LOGGER.error(error);
                        handleError("NumberFormatException in radical information", linecount, 7, 10, exception);
                    }
                } else if (line.startsWith("G  ")) {
                    try {
                        String atomNumberString = line.substring(3, 6).trim();
                        int atomNumber = Integer.parseInt(atomNumberString);
                        //String whatIsThisString = line.substring(6,9).trim();

                        String atomName = input.readLine();

                        // convert Atom into a PseudoAtom
                        IAtom prevAtom = container.getAtom(atomNumber - 1);
                        IPseudoAtom pseudoAtom = container.getBuilder().newInstance(IPseudoAtom.class, atomName);
                        if (prevAtom.getPoint2d() != null) {
                            pseudoAtom.setPoint2d(prevAtom.getPoint2d());
                        }
                        if (prevAtom.getPoint3d() != null) {
                            pseudoAtom.setPoint3d(prevAtom.getPoint3d());
                        }
                        AtomContainerManipulator.replaceAtomByAtom(container, prevAtom, pseudoAtom);
                    } catch (NumberFormatException exception) {
                        String error = "Error (" + exception.toString() + ") while parsing line " + linecount + ": " + line
                                + " in property block.";
                        LOGGER.error(error);
                        handleError("NumberFormatException in group information", linecount, 4, 7, exception);
                    }
                } else if (line.startsWith("M  RGP")) {
                    StringTokenizer st = new StringTokenizer(line);
                    //Ignore first 3 tokens (overhead).
                    st.nextToken();
                    st.nextToken();
                    st.nextToken();
                    //Process the R group numbers as defined in RGP line.
                    while (st.hasMoreTokens()) {
                        Integer position = Integer.valueOf(st.nextToken());
                        int rNumber = Integer.valueOf(st.nextToken());
                        // the container may have already had atoms before the new atoms were read
                        int index = container.getAtomCount() - nAtoms + position - 1;
                        IPseudoAtom pseudoAtom = (IPseudoAtom) container.getAtom(index);
                        if (pseudoAtom != null) {
                            pseudoAtom.setLabel("R" + rNumber);
                        }
                    }
                }
                if (line.startsWith("V  ")) {
                    Integer atomNumber = Integer.valueOf(line.substring(3, 6).trim());
                    IAtom atomWithComment = container.getAtom(atomNumber - 1);
                    atomWithComment.setProperty(CDKConstants.COMMENT, line.substring(7));
                }

                if (!lineRead) {
                    LOGGER.warn("Skipping line in property block: ", line);
                }
            }
        }

        /**
         * Read non-structural data from input and store as properties the provided
         * 'container'. Non-structural data appears in a structure data file (SDF)
         * after an Molfile and before the record deliminator ('$$$$'). The data
         * consists of one or more Data Header and Data blocks, an example is seen
         * below.
         *
         * <pre>{@code
         * > 29 <DENSITY>
         * 0.9132 - 20.0
         *
         * > 29 <BOILING.POINT>
         * 63.0 (737 MM)
         * 79.0 (42 MM)
         *
         * > 29 <ALTERNATE.NAMES>
         * SYLVAN
         *
         * > 29 <DATE>
         * 09-23-1980
         *
         * > 29 <CRC.NUMBER>
         * F-0213
         *
         * }</pre>
         *
         *
         * @param input input source
         * @param container the container
         * @throws IOException an error occur whilst reading the input
         */
        static void readNonStructuralData(final BufferedReader input, final IAtomContainer container) throws IOException {

            String line, header = null;
            boolean wrap = false;

            final StringBuilder data = new StringBuilder(80);

            while (!endOfRecord(line = input.readLine())) {

                final String newHeader = dataHeader(line);

                if (newHeader != null) {

                    if (header != null) {
                        container.setProperty(header, data.toString());
                    }

                    header = newHeader;
                    wrap = false;
                    data.setLength(0);

                } else {

                    if (data.length() > 0 || !line.equals(" ")) {
                        line = line.trim();
                    }

                    if (line.isEmpty()) {
                        continue;
                    }

                    if (!wrap && data.length() > 0) {
                        data.append(NEW_LINE);
                    }
                    data.append(line);

                    wrap = line.length() == 80;
                }
            }

            if (header != null) {
                container.setProperty(header, data.toString());
            }
        }

        /**
         * Obtain the field name from a potential SD data header. If the header does
         * not contain a field name, then null is returned. The method does not
         * currently return field numbers (e.g. DT&lt;n&gt;).
         *
         * @param line an input line
         * @return the field name
         */
        static String dataHeader(final String line) {
            if (line.length() > 2 && line.charAt(0) != '>' && line.charAt(1) != ' ') {
                return null;
            }
            int i = line.indexOf('<', 2);
            if (i < 0) {
                return null;
            }
            int j = line.indexOf('>', i);
            if (j < 0) {
                return null;
            }
            return line.substring(i + 1, j);
        }

        /**
         * Is the line the end of a record. A line is the end of a record if it is
         * 'null' or is the SDF deliminator, '$$$$'.
         *
         * @param line a line from the input
         * @return the line indicates the end of a record was reached
         */
        private static boolean endOfRecord(final String line) {
            return line == null || line.equals(RECORD_DELIMITER);
        }

        /**
         * Enumeration of property keys that can be specified in the V2000 property
         * block.
         */
        enum PropertyKey {

            /**
             * Atom Alias.
             */
            ATOM_ALIAS,
            /**
             * Atom Value.
             */
            ATOM_VALUE,
            /**
             * Group Abbreviation.
             */
            GROUP_ABBREVIATION,
            /**
             * Skip lines.
             */
            SKIP,
            /**
             * Charge [Generic].
             */
            M_CHG,
            /**
             * Radical [Generic].
             */
            M_RAD,
            /**
             * Isotope [Generic].
             */
            M_ISO,
            /**
             * Ring Bond Count [Query].
             */
            M_RBC,
            /**
             * Substitution Count [Query].
             */
            M_SUB,
            /**
             * Unsaturated Atom [Query].
             */
            M_UNS,
            /**
             * Link Atom [Query].
             */
            M_LIN,
            /**
             * Atom List [Query].
             */
            M_ALS,
            /**
             * Attachment Point [Rgroup].
             */
            M_APO,
            /**
             * Atom Attachment Order [Rgroup].
             */
            M_AAL,
            /**
             * Rgroup Label Location [Rgroup].
             */
            M_RGP,
            /**
             * Rgroup Logic, Unsatisfied Sites, Range of Occurrence [Rgroup].
             */
            M_LOG,
            /**
             * Sgroup Type [Sgroup].
             */
            M_STY,
            /**
             * Sgroup Subtype [Sgroup].
             */
            M_SST,
            /**
             * Sgroup Labels [Sgroup].
             */
            M_SLB,
            /**
             * Sgroup Connectivity [Sgroup].
             */
            M_SCN,
            /**
             * Sgroup Expansion [Sgroup].
             */
            M_SDS,
            /**
             * Sgroup Atom List [Sgroup].
             */
            M_SAL,
            /**
             * Sgroup Bond List [Sgroup].
             */
            M_SBL,
            /**
             * Multiple Group Parent Atom List [Sgroup].
             */
            M_SPA,
            /**
             * Sgroup Subscript [Sgroup].
             */
            M_SMT,
            /**
             * Sgroup Correspondence [Sgroup].
             */
            M_CRS,
            /**
             * Sgroup Display Information [Sgroup].
             */
            M_SDI,
            /**
             * Superatom Bond and Vector Information [Sgroup].
             */
            M_SBV,
            /**
             * Data Sgroup Field Description [Sgroup].
             */
            M_SDT,
            /**
             * Data Sgroup Display Information [Sgroup].
             */
            M_SDD,
            /**
             * Data Sgroup Data.
             */
            M_SCD,
            /**
             * Data Sgroup Data.
             */
            M_SED,
            /**
             * Sgroup Hierarchy Information.
             */
            M_SPL,
            /**
             * Sgroup Component Numbers.
             */
            M_SNC,
            /**
             * Sgroup Bracket Style.
             */
            M_SBT,
            /**
             * 3D Feature Properties.
             */
            M_$3D,
            /**
             * ACDLabs Atom Label
             */
            M_ZZC,
            /**
             * End of Block.
             */
            M_END,
            /**
             * Non-property header.
             */
            UNKNOWN;

            /**
             * Index of 'M XXX' properties for quick lookup.
             */
            private static final Map<String, PropertyKey> mSuffix = new HashMap<String, PropertyKey>(60);

            static {
                for (PropertyKey p : values()) {
                    if (p.name().charAt(0) == 'M') {
                        mSuffix.put(p.name().substring(2, 5), p);
                    }
                }
            }

            /**
             * Determine the property key of the provided line.
             *
             * @param line an property line
             * @return the key (defaults to {@link #UNKNOWN})
             */
            static PropertyKey of(final String line) {
                if (line.length() < 5) {
                    return UNKNOWN;
                }
                switch (line.charAt(0)) {
                    case 'A':
                        if (line.charAt(1) == ' ' && line.charAt(2) == ' ') {
                            return ATOM_ALIAS;
                        }
                        return UNKNOWN;
                    case 'G':
                        if (line.charAt(1) == ' ' && line.charAt(2) == ' ') {
                            return GROUP_ABBREVIATION;
                        }
                        return UNKNOWN;
                    case 'S':
                        if (line.charAt(1) == ' ' && line.charAt(2) == ' ') {
                            return SKIP;
                        }
                        return UNKNOWN;
                    case 'V':
                        if (line.charAt(1) == ' ' && line.charAt(2) == ' ') {
                            return ATOM_VALUE;
                        }
                        return UNKNOWN;
                    case 'M':
                        if (line.charAt(1) != ' ' || line.charAt(2) != ' ') {
                            return UNKNOWN;
                        }
                        PropertyKey property = mSuffix.get(line.substring(3, 6));
                        if (property != null) {
                            return property;
                        }
                        return UNKNOWN;
                }
                return UNKNOWN;
            }

        }

        /**
         * Defines the version of the CTab.
         */
        enum CTabVersion {
            V2000, V3000, UNSPECIFIED;

            /**
             * Given a CTab header, what version was specified. The version is
             * identifier in the by the presence of 'V[2|3]000'. If not version tag
             * is present the version is unspecified.
             *
             * <pre>  5  5  0  0  0  0            999 V2000</prev>
             * <pre>  0  0  0  0  0  0            999 V3000</prev>
             *
             * @param header input line (non-null)
             * @return the CTab version
             */
            static CTabVersion ofHeader(String header) {
                if (header.length() < 39) {
                    return UNSPECIFIED;
                }
                char c = header.charAt(34);
                if (c != 'v' && c != 'V') {
                    return UNSPECIFIED;
                }
                if (header.charAt(35) == '2') // could check for '000'
                {
                    return V2000;
                }
                if (header.charAt(35) == '3') // could check for '000'
                {
                    return V3000;
                }
                return UNSPECIFIED;
            }
        }

    }



    /**
     * Reads a molecule from an MDL RXN file {
     *
     * @cdk.cite DAL92}. This MDL RXN reader uses the MDLV2000 reader to read each
     * mol file
     * @cdk.module io
     * @cdk.githash
     * @cdk.iooptions
     *
     * @author Egon Willighagen
     * @author Thomas Kuhn
     * @cdk.created 2003-07-24
     *
     * @cdk.keyword file format, MDL RXN
     * @cdk.bug 1849923
     */
    public static class MDLRXNV2000Reader extends DefaultChemObjectReader {

        BufferedReader input = null;
        private static ILoggingTool logger = LoggingToolFactory.createLoggingTool(MDLRXNV2000Reader.class);

        /**
         * Constructs a new MDLReader that can read Molecule from a given Reader.
         *
         * @param in The Reader to read from
         */
        public MDLRXNV2000Reader(Reader in) {
            this(in, IChemObjectReader.Mode.RELAXED);
        }

        public MDLRXNV2000Reader(Reader in, IChemObjectReader.Mode mode) {
            if (in instanceof BufferedReader) {
                input = (BufferedReader) in;
            } else {
                input = new BufferedReader(in);
            }
            super.mode = mode;
        }

        public MDLRXNV2000Reader(InputStream input) {
            this(input, IChemObjectReader.Mode.RELAXED);
        }

        public MDLRXNV2000Reader(InputStream input, IChemObjectReader.Mode mode) {
            this(new InputStreamReader(input), mode);
        }

        public MDLRXNV2000Reader() {
            this(new StringReader(""));
        }

        @Override
        public IResourceFormat getFormat() {
            return MDLRXNFormat.getInstance();
        }

        @Override
        public void setReader(Reader input) throws CDKException {
            if (input instanceof BufferedReader) {
                this.input = (BufferedReader) input;
            } else {
                this.input = new BufferedReader(input);
            }
        }

        @Override
        public void setReader(InputStream input) throws CDKException {
            setReader(new InputStreamReader(input));
        }

        @Override
        public boolean accepts(Class<? extends IChemObject> classObject) {
            if (IChemFile.class.equals(classObject)) {
                return true;
            }
            if (IChemModel.class.equals(classObject)) {
                return true;
            }
            if (IReaction.class.equals(classObject)) {
                return true;
            }
            Class<?>[] interfaces = classObject.getInterfaces();
            for (Class<?> intf : interfaces) {
                if (IChemModel.class.equals(intf)) {
                    return true;
                }
                if (IChemFile.class.equals(intf)) {
                    return true;
                }
                if (IReaction.class.equals(intf)) {
                    return true;
                }
            }
            Class superClass = classObject.getSuperclass();
            if (superClass != null) {
                return this.accepts(superClass);
            }
            return false;
        }

        /**
         * Takes an object which subclasses IChemObject, e.g.Molecule, and will read
         * this (from file, database, internet etc).If the specific implementation
         * does not support a specific IChemObject it will throw an Exception.
         *
         * @param <T>
         * @param object The object that subclasses IChemObject
         * @return The IChemObject read
         * @exception CDKException
         */
        @Override
        public <T extends IChemObject> T read(T object) throws CDKException {
            if (object instanceof IReaction) {
                return (T) readReaction(object.getBuilder());
            } else if (object instanceof IReactionSet) {
                IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
                reactionSet.addReaction(readReaction(object.getBuilder()));
                return (T) reactionSet;
            } else if (object instanceof IChemModel) {
                IChemModel model = object.getBuilder().newInstance(IChemModel.class);
                IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
                reactionSet.addReaction(readReaction(object.getBuilder()));
                model.setReactionSet(reactionSet);
                return (T) model;
            } else if (object instanceof IChemFile) {
                IChemFile chemFile = object.getBuilder().newInstance(IChemFile.class);
                IChemSequence sequence = object.getBuilder().newInstance(IChemSequence.class);
                sequence.addChemModel((IChemModel) read(object.getBuilder().newInstance(IChemModel.class)));
                chemFile.addChemSequence(sequence);
                return (T) chemFile;
            } else {
                throw new CDKException("Only supported are Reaction and ChemModel, and not " + object.getClass().getName()
                        + ".");
            }
        }

        public boolean accepts(IChemObject object) {
            if (object instanceof IReaction) {
                return true;
            } else if (object instanceof IChemModel) {
                return true;
            } else if (object instanceof IChemFile) {
                return true;
            } else if (object instanceof IReactionSet) {
                return true;
            }
            return false;
        }

        /**
         * Read a Reaction from a file in MDL RXN format
         *
         * @return The Reaction that was read from the MDL file.
         */
        private IReaction readReaction(IChemObjectBuilder builder) throws CDKException {
            IReaction reaction = builder.newInstance(IReaction.class);
            try {
                input.readLine(); // first line should be $RXN
                input.readLine(); // second line
                input.readLine(); // third line
                input.readLine(); // fourth line
            } catch (IOException exception) {
                logger.debug(exception);
                throw new CDKException("Error while reading header of RXN file", exception);
            }

            int numReactans = 0;
            int numProducts = 0;
            int agentCount = 0;
            try {
                String countsLine = input.readLine();
                /*
                 * this line contains the number of reactants and products
                 */
                StringTokenizer tokenizer = new StringTokenizer(countsLine);
                numReactans = Integer.valueOf(tokenizer.nextToken());
                logger.info("Expecting " + numReactans + " reactants in file");
                numProducts = Integer.valueOf(tokenizer.nextToken());
                if (tokenizer.hasMoreTokens()) {
                    agentCount = Integer.valueOf(tokenizer.nextToken());
                    // ChemAxon extension, technically BIOVIA now support this but
                    // not documented yet
                    if (mode == IChemObjectReader.Mode.STRICT && agentCount > 0) {
                        throw new CDKException("RXN files uses agent count extension");
                    }
                }
                logger.info("Expecting " + numProducts + " products in file");
            } catch (IOException | NumberFormatException exception) {
                logger.debug(exception);
                throw new CDKException("Error while counts line of RXN file", exception);
            }

            // now read the molecules
            try {
                String line = input.readLine();
                if (line == null || !line.startsWith("$MOL")) {
                    throw new CDKException("Expected $MOL to start, was" + line);
                }

                List<IAtomContainer> components = new ArrayList<>();

                StringBuilder sb = new StringBuilder();
                while ((line = input.readLine()) != null) {
                    if (line.startsWith("$MOL")) {
                        processMol(builder.newAtomContainer(), components, sb);
                        sb.setLength(0);
                    } else {
                        sb.append(line).append('\n');
                    }
                }

                // last record
                if (sb.length() > 0) {
                    processMol(builder.newAtomContainer(), components, sb);
                }

                for (IAtomContainer component : components.subList(0, numReactans)) {
                    reaction.addReactant(component);
                }
                for (IAtomContainer component : components.subList(numReactans,
                        numReactans + numProducts)) {
                    reaction.addProduct(component);
                }
                for (IAtomContainer component : components.subList(numReactans + numProducts,
                        components.size())) {
                    reaction.addAgent(component);
                }

            } catch (CDKException exception) {
                // rethrow exception from MDLReader
                throw exception;
            } catch (IOException | IllegalArgumentException exception) {
                logger.debug(exception);
                throw new CDKException("Error while reading reactant", exception);
            }

            // now try to map things, if wanted
            logger.info("Reading atom-atom mapping from file");
            // distribute all atoms over two AtomContainer's
            IAtomContainer reactingSide = builder.newInstance(IAtomContainer.class);
            Iterator<IAtomContainer> molecules = reaction.getReactants().atomContainers().iterator();
            while (molecules.hasNext()) {
                reactingSide.add(molecules.next());
            }
            IAtomContainer producedSide = builder.newInstance(IAtomContainer.class);
            molecules = reaction.getProducts().atomContainers().iterator();
            while (molecules.hasNext()) {
                producedSide.add(molecules.next());
            }

            // map the atoms
            int mappingCount = 0;
            //        IAtom[] reactantAtoms = reactingSide.getAtoms();
            //        IAtom[] producedAtoms = producedSide.getAtoms();
            for (int i = 0; i < reactingSide.getAtomCount(); i++) {
                for (int j = 0; j < producedSide.getAtomCount(); j++) {
                    IAtom eductAtom = reactingSide.getAtom(i);
                    IAtom productAtom = producedSide.getAtom(j);
                    if (eductAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null
                            && eductAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING).equals(
                                    productAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING))) {
                        reaction.addMapping(builder.newInstance(IMapping.class, eductAtom, productAtom));
                        mappingCount++;
                        break;
                    }
                }
            }
            logger.info("Mapped atom pairs: " + mappingCount);

            return reaction;
        }

        private void processMol(IAtomContainer mol, List<IAtomContainer> components, StringBuilder sb) throws CDKException, IOException {
            try (MDLV2000Reader reader = new MDLV2000Reader(new StringReader(sb.toString()), super.mode)) {
                components.add(reader.read(mol));
            }
        }

        @Override
        public void close() throws IOException {
            input.close();
        }
    }



    /**
     * Writes a reaction to a MDL rxn or SDF file. Attention: Stoichiometric
     * coefficients have to be natural numbers.
     *
     * <pre>
     * MDLRXNWriter writer = new MDLRXNWriter(new FileWriter(new File("output.mol")));
     * writer.write((AtomContainer)molecule);
     * writer.close();
     * </pre>
     *
     * See {
     *
     * @cdk.cite DAL92}.
     *
     * @cdk.module io
     *
     *
     * @cdk.keyword file format, MDL RXN file
     */
    public static class MDLV2000RXNWriter extends DefaultChemObjectWriter {

        private static ILoggingTool LOGGER = createLoggingTool(MDLV2000RXNWriter.class);
        private BufferedWriter writer;
        private int reactionNumber;

        /**
         *
         */
        public Map rdFields = null;

        /**
         * Constructs a new MDLWriter that can write an array of Molecules to a
         * Writer.
         *
         * @param out The Writer to write to
         */
        public MDLV2000RXNWriter(Writer out) {
            try {
                if (out instanceof BufferedWriter) {
                    writer = (BufferedWriter) out;
                } else {
                    writer = new BufferedWriter(out);
                }
            } catch (Exception ex) {
                LOGGER.error(ex);
            }
            this.reactionNumber = 1;
        }

        /**
         * Constructs a new MDLWriter that can write an array of Molecules to a
         * given OutputStream.
         *
         * @param output The OutputStream to write to
         */
        public MDLV2000RXNWriter(OutputStream output) {
            this(new OutputStreamWriter(output));
        }

        /**
         *
         */
        public MDLV2000RXNWriter() {
            this(new StringWriter());
        }

        /**
         *
         * @return
         */
        @Override
        public IResourceFormat getFormat() {
            return getInstance();
        }

        /**
         *
         * @param out
         * @throws CDKException
         */
        @Override
        public void setWriter(Writer out) throws CDKException {
            if (out instanceof BufferedWriter) {
                writer = (BufferedWriter) out;
            } else {
                writer = new BufferedWriter(out);
            }
        }

        /**
         *
         * @param output
         * @throws CDKException
         */
        @Override
        public void setWriter(OutputStream output) throws CDKException {
            setWriter(new OutputStreamWriter(output));
        }

        /**
         * Here you can set a map which will be used to build rd fields in the file.
         * The entries will be translated to rd fields like this:<br>
         * &gt; &lt;key&gt;<br>
         * &gt; value<br>
         * empty line<br>
         *
         * @param map The map to be used, map of String-String pairs
         */
        public void setRdFields(Map map) {
            rdFields = map;
        }

        /**
         * Flushes the output and closes this object.
         *
         * @throws java.io.IOException
         */
        @Override
        public void close() throws IOException {
            writer.close();
        }

        /**
         *
         * @param classObject
         * @return
         */
        @Override
        public boolean accepts(Class classObject) {
            Class[] interfaces = classObject.getInterfaces();
            for (Class intf : interfaces) {
                if (IReaction.class.equals(intf)) {
                    return true;
                }
                if (IReactionSet.class.equals(intf)) {
                    return true;
                }
            }
            Class superClass = classObject.getSuperclass();
            if (superClass != null) {
                return this.accepts(superClass);
            }
            return false;
        }

        /**
         * Writes a IChemObject to the MDL RXN file formated output. It can only
         * output ChemObjects of type Reaction
         *
         * @param object class must be of type AtomContainer or MoleculeSet.
         * @throws org.openscience.cdk.exception.CDKException
         *
         * @see org.openscience.cdk.ChemFile
         */
        @Override
        public void write(IChemObject object) throws CDKException {
            if (object instanceof IReactionSet) {
                writeReactionSet((IReactionSet) object);
            } else if (object instanceof IReaction) {
                writeReaction((IReaction) object);
            } else {
                throw new CDKException("Only supported is writing ReactionSet, Reaction objects.");
            }
        }

        /**
         * Writes an array of Reaction to an OutputStream in MDL rdf format.
         *
         * @param som Array of Reactions that is written to an OutputStream
         */
        private void writeReactionSet(IReactionSet reactions) throws CDKException {

            for (Iterator<IReaction> it = reactions.reactions().iterator(); it.hasNext();) {
                writeReaction(it.next());
            }
        }

        /**
         * Writes a Reaction to an OutputStream in MDL sdf format.
         *
         * @param reaction A Reaction that is written to an OutputStream
         */
        private void writeReaction(IReaction reaction) throws CDKException {

            /*Fixed correct reactant product count*/
            int reactantCount = 0;
            for (IAtomContainer e : reaction.getReactants().atomContainers()) {
                reactantCount += reaction.getReactantCoefficient(e).intValue();
            }
            int productCount = 0;
            for (IAtomContainer p : reaction.getProducts().atomContainers()) {
                productCount += reaction.getProductCoefficient(p).intValue();
            }

            if (reactantCount <= 0 || productCount <= 0) {
                throw new CDKException("Either no reactants or no products present.");
            }

            try {
                // taking care of the $$$$ signs:
                // we do not write such a sign at the end of the first reaction, thus we have to write on BEFORE the second reaction
                if (reactionNumber == 2) {
                    writer.write("$$$$");
                    writer.newLine();
                }
                writer.write("$RXN");
                writer.newLine();

                // reaction name
                String line = (String) reaction.getProperty(TITLE);
                if (line == null) {
                    String rid = reaction.getID() == null ? "" : reaction.getID();
                    line = "  " + "EC-BLAST" + "     " + rid;
                }
                if (line.length() > 80) {
                    line = line.substring(0, 80);
                }
                writer.newLine();
                writer.write(line);
                // user/program/date&time/reaction registry no. line
                writer.newLine();
                // comment line
                line = (String) reaction.getProperty(REMARK);
                if (line == null) {
                    line = "";
                }
                if (line.length() > 80) {
                    line = line.substring(0, 80);
                }
                writer.write(line);
                writer.newLine();

                line = "";
                line += formatMDLInt(reactantCount, 3);
                line += formatMDLInt(productCount, 3);
                writer.write(line);
                writer.newLine();

                int i = 0;
                for (IMapping mapping : reaction.mappings()) {
                    Iterator<IChemObject> it = mapping.relatedChemObjects().iterator();
                    /*
                     Do not overwrite the existing labels
                     */

                    if (it.next().getProperty(ATOM_ATOM_MAPPING) == null) {
                        it.next().setProperty(ATOM_ATOM_MAPPING, i + 1);
                        it.next().setProperty(ATOM_ATOM_MAPPING, i + 1);
                        i++;
                    }
                }
                writeMoleculeSet(reaction.getReactants());
                writeMoleculeSet(reaction.getProducts());

                //write sdfields, if any
                if (rdFields != null) {
                    Set set = rdFields.keySet();
                    Iterator iterator = set.iterator();
                    while (iterator.hasNext()) {
                        Object element = iterator.next();
                        writer.write("> <" + element + ">");
                        writer.newLine();
                        writer.write(rdFields.get(element).toString());
                        writer.newLine();
                        writer.newLine();
                    }
                }
                // taking care of the $$$$ signs:
                // we write such a sign at the end of all except the first molecule
                if (reactionNumber != 1) {
                    writer.write("$$$$");
                    writer.newLine();
                }
                reactionNumber++;

            } catch (IOException ex) {
                LOGGER.error(ex.getMessage());
                LOGGER.debug(ex);
                throw new CDKException("Exception while writing MDL file: " + ex.getMessage(), ex);
            }
        }

        /**
         * Writes a MoleculeSet to an OutputStream for the reaction.
         *
         * @param som The MoleculeSet that is written to an OutputStream
         */
        private void writeMoleculeSet(IAtomContainerSet som) throws IOException, CDKException {

            for (int i = 0; i < som.getAtomContainerCount(); i++) {
                IAtomContainer mol = som.getAtomContainer(i);
                for (int j = 0; j < som.getMultiplier(i); j++) {
                    StringWriter sw = new StringWriter();
                    writer.write("$MOL");
                    writer.newLine();
                    MDLV2000Writer mdlwriter = null;
                    try {
                        mdlwriter = new MDLV2000Writer(sw);
                        // GMT : added to force writing of aromatic bond types
                        // Asad: Swith off aromatic bond writing
                        // mdlwriter.getIOSettings()[1].setSetting("true");
                    } catch (Exception ex) {
                        LOGGER.error(ex.getMessage());
                        LOGGER.debug(ex);
                        throw new CDKException("Exception while creating MDLWriter: " + ex.getMessage(), ex);
                    }
                    mdlwriter.write(mol);
                    writer.write(sw.toString());
                }
            }
        }

        /**
         * Formats an int to fit into the connection table and changes it to a
         * String.
         *
         * @param i The int to be formated
         * @param l Length of the String
         * @return The String to be written into the connection table
         */
        private String formatMDLInt(int i, int l) {
            String s = "", fs = "";
            NumberFormat nf = getNumberInstance(ENGLISH);
            nf.setParseIntegerOnly(true);
            nf.setMinimumIntegerDigits(1);
            nf.setMaximumIntegerDigits(l);
            nf.setGroupingUsed(false);
            s = nf.format(i);
            l -= s.length();
            for (int f = 0; f < l; f++) {
                fs += " ";
            }
            fs += s;
            return fs;
        }
    }


}
