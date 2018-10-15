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
import static java.awt.Color.RGBtoHSB;
import static java.awt.Color.getHSBColor;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple utility class to generate a 'ramp' of colors between two values.
 *
 * @author maclean
 *
 */
public class ColorRamp {

    /**
     * Get N colors as a list.
     *
     * @param number the number of colors to generate
     * @return a list of colors.
     */
    public static List<Color> getColors(int number) {
        List<Color> colors = new ArrayList<>();
        for (int i = 0; i < number; i++) {
            colors.add(colorRamp(i, 0, number));
        }
        return colors;
    }

    /**
     * Get a color for a value 'v' between vmin and vmax.
     *
     * @param v the point on the ramp to make a color for
     * @param vmin the minimum value in the range
     * @param vmax the maximum value in the range
     * @return the color for v
     */
    public static Color colorRamp(int v, int vmin, int vmax) {
        double r = 1.0;
        double g = 1.0;
        double b = 1.0;
        if (v < vmin) {
            v = vmin;
        }
        if (v > vmax) {
            v = vmax;
        }
        int dv = vmax - vmin;

        try {
            if (v < (vmin + 0.25 * dv)) {
                r = 0.0;
                g = 4.0 * (v - vmin) / dv;
            } else if (v < (vmin + 0.5 * dv)) {
                r = 0.0;
                b = 1.0 + 4.0 * (vmin + 0.25 * dv - v) / dv;
            } else if (v < (vmin + 0.75 * dv)) {
                r = 4.0 * (v - vmin - 0.5 * dv) / dv;
                b = 0.0;
            } else {
                g = 1.0 + 4.0 * (vmin + 0.75 * dv - v) / dv;
                b = 0.0;
            }
            float[] hsb = RGBtoHSB(
                    (int) (r * 255), (int) (g * 255), (int) (b * 255), null);
            return getHSBColor(hsb[0], hsb[1], hsb[2]);
        } catch (ArithmeticException zde) {
            float[] hsb = RGBtoHSB(0, 0, 0, null);
            return getHSBColor(hsb[0], hsb[1], hsb[2]);
        }
    }

    private ColorRamp() {
    }

}
