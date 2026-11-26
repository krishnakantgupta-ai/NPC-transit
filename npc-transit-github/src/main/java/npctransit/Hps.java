package npctransit;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

/**
 * Residue parameters for the HPS (hydrophobicity scale) coarse-grained model, Kapcha-Rossky-based lambda values
 * as used by Dignon et al. 2018 (PLoS Comput Biol 14:e1005941). VERIFY these against the original paper before
 * publication; the charge model is formal (Lys/Arg +1, Asp/Glu -1, His 0).
 */
public class Hps {
    static final String AA = "ARNDCQEGHILKMFPSTWYV";
    static final double[] LAMBDA = {0.730,0.000,0.432,0.378,0.595,0.514,0.459,0.649,0.514,0.973,0.973,0.514,0.838,1.000,1.000,0.595,0.676,0.946,0.865,0.892};
    static final double[] CHARGE = {0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0};

    public static double lambda(char c) { int i = AA.indexOf(c); return i < 0 ? 0.5 : LAMBDA[i]; }
    public static double charge(char c) { int i = AA.indexOf(c); return i < 0 ? 0 : CHARGE[i]; }
    public static double meanLambda(CharSequence s) { double t = 0; for (int i = 0; i < s.length(); i++) t += lambda(s.charAt(i)); return t / s.length(); }

    /** Reads FASTA; the FIRST residue of each sequence is the grafted (anchored) end. */
    public static List<String> readFasta(Path p) throws IOException {
        List<String> out = new ArrayList<>(); StringBuilder cur = null;
        for (String l : Files.readAllLines(p)) {
            l = l.trim();
            if (l.startsWith(">")) { if (cur != null && cur.length() > 0) out.add(cur.toString()); cur = new StringBuilder(); }
            else if (!l.isEmpty() && cur != null) cur.append(l.toUpperCase().replaceAll("[^A-Z]", ""));
        }
        if (cur != null && cur.length() > 0) out.add(cur.toString());
        if (out.isEmpty()) throw new IOException("No sequences in " + p);
        return out;
    }
}
