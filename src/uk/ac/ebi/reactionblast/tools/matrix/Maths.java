package uk.ac.ebi.reactionblast.tools.matrix;

import java.util.logging.Logger;

public class Maths {

    private static final Logger LOG = Logger.getLogger(Maths.class.getName());

    /**
     * sqrt(a^2 + b^2) without under/overflow. * * @param a
     *
     *
     * @param b
     * @return
     * @param a
     */
    public static Double hypot(Double a, Double b) {
        Double r;
        if (Math.abs(a) > Math.abs(b)) {
            r = b / a;
            r = Math.abs(a) * Math.sqrt(1 + r * r);
        } else if (b != 0) {
            r = a / b;
            r = Math.abs(b) * Math.sqrt(1 + r * r);
        } else {
            r = 0.0;
        }
        return r;
    }
}
