package npctransit;

import java.util.*;

/**
 * Heuristic FG-motif probe docking. A phenylalanine-ring-sized probe (sphere) is scanned over a 1 A grid
 * around the protein. Each grid point is scored by hydrophobic contacts (side-chain C / S) minus a penalty
 * for polar contacts (N / O), with a steric-clash filter. This approximates the Phe-pocket binding mode seen
 * in FG-binding proteins (NTF2-like, importin-beta) but is NOT a validated docking engine: scores are
 * relative, unitless, and must be benchmarked on known FG-binding vs. non-binding proteins.
 */
public class FgProbeDock {
    public static class Site {
        public double[] pos; public double score; public int hyd, polar; public List<String> residues = new ArrayList<>();
    }

    static boolean hydrophobic(String el, String name) {
        return (el.equals("C") && !name.equals("C") && !name.equals("CA")) || el.equals("S");
    }

    public static List<Site> run(PdbParser.Structure s, int topN) {
        int n = s.xyz.size();
        double cell = 8.0;
        double[] lo = {1e9, 1e9, 1e9}, hi = {-1e9, -1e9, -1e9};
        for (double[] p : s.xyz) for (int i = 0; i < 3; i++) { lo[i] = Math.min(lo[i], p[i]); hi[i] = Math.max(hi[i], p[i]); }
        Map<Long, List<Integer>> grid = new HashMap<>();
        for (int a = 0; a < n; a++) grid.computeIfAbsent(key(s.xyz.get(a), lo, cell), k -> new ArrayList<>()).add(a);
        boolean[] hyd = new boolean[n], pol = new boolean[n];
        for (int a = 0; a < n; a++) { hyd[a] = hydrophobic(s.elem.get(a), s.atomName.get(a)); pol[a] = s.elem.get(a).equals("N") || s.elem.get(a).equals("O"); }

        List<Site> cand = new ArrayList<>();
        for (double x = lo[0] - 6; x <= hi[0] + 6; x += 1.0)
        for (double y = lo[1] - 6; y <= hi[1] + 6; y += 1.0)
        for (double z = lo[2] - 6; z <= hi[2] + 6; z += 1.0) {
            double[] q = {x, y, z};
            int ci = (int) Math.floor((x - lo[0]) / cell), cj = (int) Math.floor((y - lo[1]) / cell), ck = (int) Math.floor((z - lo[2]) / cell);
            double minD = 1e9; int h = 0, pl = 0, bur = 0;
            for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++) for (int dk = -1; dk <= 1; dk++) {
                List<Integer> L = grid.get(pack(ci + di, cj + dj, ck + dk)); if (L == null) continue;
                for (int a : L) {
                    double[] p = s.xyz.get(a);
                    double d = Math.sqrt((p[0]-x)*(p[0]-x) + (p[1]-y)*(p[1]-y) + (p[2]-z)*(p[2]-z));
                    minD = Math.min(minD, d);
                    if (d < 8.0) bur++;
                    if (hyd[a] && d >= 3.4 && d <= 5.0) h++;
                    if (pol[a] && d >= 3.2 && d <= 4.0) pl++;
                }
            }
            if (minD < 3.4 || bur < 35 || h < 4) continue;   // clash / too exposed / too few hydrophobic contacts
            Site st = new Site(); st.pos = q; st.hyd = h; st.polar = pl; st.score = -h + 0.5 * pl;
            cand.add(st);
        }
        cand.sort(Comparator.comparingDouble(t -> t.score));
        List<Site> out = new ArrayList<>();
        for (Site c : cand) {
            boolean ok = true;
            for (Site o : out) if (dist(o.pos, c.pos) < 5.0) { ok = false; break; }
            if (!ok) continue;
            Set<String> res = new TreeSet<>();
            for (int a = 0; a < n; a++) if (dist(s.xyz.get(a), c.pos) <= 5.0) res.add(s.resLabel.get(a));
            c.residues.addAll(res); out.add(c);
            if (out.size() >= topN) break;
        }
        return out;
    }
    static double dist(double[] a, double[] b) { return Math.sqrt(Math.pow(a[0]-b[0],2)+Math.pow(a[1]-b[1],2)+Math.pow(a[2]-b[2],2)); }
    static long key(double[] p, double[] lo, double cell) {
        return pack((int) Math.floor((p[0]-lo[0])/cell), (int) Math.floor((p[1]-lo[1])/cell), (int) Math.floor((p[2]-lo[2])/cell));
    }
    static long pack(int i, int j, int k) { return ((long) (i + 512) << 40) | ((long) (j + 512) << 20) | (k + 512); }
}
