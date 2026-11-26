package npctransit;

import java.util.*;
import java.util.regex.Pattern;

/** Sequence/structure-derived descriptors computed from CA atoms. */
public class Features {
    /** Features used by the passive-permeability model (order matters). */
    public static final String[] NAMES = {"massKDa", "rgA", "hydrophobicFraction", "absNetChargePerResidue"};
    static final Set<String> HYDROPHOBIC = Set.of("ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO");
    static final Pattern NLS = Pattern.compile("K[KR].[KR]|[KR]{4,}");

    public double massKDa, rg, stokesApprox, hydFrac, netCharge, absChargePerRes;
    public int residues, chains;
    public boolean nlsLike;

    public double get(String n) {
        switch (n) {
            case "massKDa": return massKDa; case "rgA": return rg; case "stokesA": return stokesApprox;
            case "hydrophobicFraction": return hydFrac; case "absNetChargePerResidue": return absChargePerRes; case "residues": return residues;
            default: throw new IllegalArgumentException("Unknown feature: " + n + " (massKDa, rgA, stokesA, hydrophobicFraction, absNetChargePerResidue, residues)");
        }
    }
    public double[] vector() { return new double[]{massKDa, rg, hydFrac, absChargePerRes}; }

    public static Features compute(PdbParser.Structure s) {
        Features f = new Features();
        int n = s.ca.size(), pos = 0, neg = 0, hyd = 0;
        double mass = 18.02;
        double[] c = new double[3];
        for (double[] p : s.ca) for (int i = 0; i < 3; i++) c[i] += p[i] / n;
        double rg2 = 0;
        for (double[] p : s.ca) for (int i = 0; i < 3; i++) rg2 += Math.pow(p[i] - c[i], 2) / n;
        for (String r : s.resNames) {
            mass += PdbParser.MASS.get(r);
            if (HYDROPHOBIC.contains(r)) hyd++;
            if (r.equals("LYS") || r.equals("ARG")) pos++;
            if (r.equals("ASP") || r.equals("GLU")) neg++;
        }
        f.residues = n; f.chains = s.chains;
        f.massKDa = mass / 1000.0;
        f.rg = Math.sqrt(rg2);
        f.stokesApprox = f.rg / 0.775;   // Rh ~ Rg/0.775 for compact globular proteins (approximation)
        f.hydFrac = (double) hyd / n;
        f.netCharge = pos - neg;         // formal, pH-independent count (His neutral)
        f.absChargePerRes = Math.abs(f.netCharge) / n;
        f.nlsLike = NLS.matcher(s.sequence).find();
        return f;
    }
}
