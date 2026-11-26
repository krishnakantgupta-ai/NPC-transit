package npctransit;

import java.util.*;
import java.util.concurrent.*;

/**
 * Coarse-grained FG-Nup brush in a cylindrical pore (1 bead/residue, reduced units: kT = 1, sigma = 0.6 nm).
 * Two modes:
 *  - generic: FXFG motif chains, F-F attraction epsFG (qualitative);
 *  - HPS: user-supplied FG-Nup sequences, Ashbaugh-Hatch hydropathy interactions + Debye-Hueckel electrostatics.
 * A rigid spherical cargo is held on the axis at successive z (thermodynamic integration of the mean axial force).
 * Independent replicas (different seeds) run in parallel; errors are across replicas.
 */
public class FgBrushSim {
    public static final double EPS0 = 0.336, RC = 4.0, WCA = Math.pow(2, 1.0 / 6);

    public static class Params {
        public int chains = 16, len = 20, equilFirst = 15000, equil = 4000, prod = 15000, windows = 15;
        public double anchorZ = 3.0;
        public double pipeR = 8.0, zMin = -18, zMax = 18, epsFG = 1.5, dt = 0.002;
        public double epsScale = 1.0, xi = 1.33, lB = 1.17;     // HPS: eps multiplier, Debye length & Bjerrum length in sigma
        public List<String> seqs = null;                          // non-null => HPS mode
        public double cargoEps = 1.5, cargoLambda = 0.6, cargoCharge = 0;
        public long seed = 1;
        public Params copy() {
            Params q = new Params();
            q.chains = chains; q.len = len; q.equilFirst = equilFirst; q.equil = equil; q.prod = prod; q.windows = windows;
            q.anchorZ = anchorZ; q.pipeR = pipeR; q.zMin = zMin; q.zMax = zMax; q.epsFG = epsFG; q.dt = dt; q.epsScale = epsScale; q.xi = xi; q.lB = lB;
            q.seqs = seqs; q.cargoEps = cargoEps; q.cargoLambda = cargoLambda; q.cargoCharge = cargoCharge; q.seed = seed; return q;
        }
    }
    public static class Result {
        public boolean fits = true; public double[] z, force, pmf; public double barrier, barrierSd, closure;
    }
    public static class Agg {
        public boolean fits = true; public double[] z, meanPmf; public double[] barriers, closures;
        public double mean, sd, sem, ci95; public int badClosure;
    }

    final Params P; final Random rnd; final boolean hps; final int N;
    final double[][] x, f; final boolean[] isH; final double[] lam, q; final int[] chainOf, idxIn;
    final double bond; int[][] nbr; int[] nbrCount; double fCargoZ; int[] head = new int[0], nxt;

    FgBrushSim(Params p) {
        P = p; rnd = new Random(p.seed); hps = p.seqs != null; bond = hps ? 0.633 : 1.0;
        List<String> chainSeq = new ArrayList<>();
        for (int c = 0; c < p.chains; c++) chainSeq.add(hps ? p.seqs.get(c % p.seqs.size()) : null);
        int n = 0; for (String s : chainSeq) n += hps ? s.length() : p.len; N = n;
        x = new double[N][3]; f = new double[N][3]; isH = new boolean[N]; lam = new double[N]; q = new double[N];
        chainOf = new int[N]; idxIn = new int[N]; int b = 0;
        for (int c = 0; c < p.chains; c++) {
            int L = hps ? chainSeq.get(c).length() : p.len;
            double th = 2 * Math.PI * c / p.chains + 0.1 * rnd.nextGaussian(), z0 = -p.anchorZ + 2 * p.anchorZ * rnd.nextDouble();
            double[] cur = {(p.pipeR - 0.5) * Math.cos(th), (p.pipeR - 0.5) * Math.sin(th), z0};
            for (int i = 0; i < L; i++, b++) {
                chainOf[b] = c; idxIn[b] = i;
                if (hps) { lam[b] = Hps.lambda(chainSeq.get(c).charAt(i)); q[b] = Hps.charge(chainSeq.get(c).charAt(i)); }
                else isH[b] = i >= 1 && ((i - 1) % 8 == 0 || (i - 1) % 8 == 2);
                x[b] = cur.clone();
                double rho = Math.hypot(cur[0], cur[1]) + 1e-9;
                double[] d = {-0.6 * cur[0] / rho + 0.8 * rnd.nextGaussian(), -0.6 * cur[1] / rho + 0.8 * rnd.nextGaussian(), 0.8 * rnd.nextGaussian()};
                double nd = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                for (int k = 0; k < 3; k++) cur[k] += bond * d[k] / nd;
                double rr = Math.hypot(cur[0], cur[1]); if (rr > p.pipeR - 0.7) { cur[0] *= (p.pipeR - 0.7) / rr; cur[1] *= (p.pipeR - 0.7) / rr; }
            }
        }
        nbr = new int[N][128]; nbrCount = new int[N];
    }

    /** Cell-list neighbor build, O(N). */
    void buildNeighbors() {
        double cut = (hps ? RC : 2.5) + 0.8, c2 = cut * cut;
        double[] lo = {1e18, 1e18, 1e18}, hi = {-1e18, -1e18, -1e18};
        for (double[] v : x) for (int k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], v[k]); hi[k] = Math.max(hi[k], v[k]); }
        int[] n3 = new int[3]; for (int k = 0; k < 3; k++) n3[k] = Math.max(1, Math.min(200, (int) ((hi[k] - lo[k]) / cut) + 1));
        int nc = n3[0] * n3[1] * n3[2];
        if (head.length < nc) head = new int[nc];
        Arrays.fill(head, 0, nc, -1); if (nxt == null) nxt = new int[N];
        int[] cx = new int[N], cy = new int[N], cz = new int[N];
        for (int i = 0; i < N; i++) {
            cx[i] = Math.min(n3[0] - 1, (int) ((x[i][0] - lo[0]) / cut)); cy[i] = Math.min(n3[1] - 1, (int) ((x[i][1] - lo[1]) / cut)); cz[i] = Math.min(n3[2] - 1, (int) ((x[i][2] - lo[2]) / cut));
            int c = (cx[i] * n3[1] + cy[i]) * n3[2] + cz[i]; nxt[i] = head[c]; head[c] = i;
        }
        for (int i = 0; i < N; i++) {
            int cnt = 0;
            for (int a = Math.max(0, cx[i] - 1); a <= Math.min(n3[0] - 1, cx[i] + 1); a++)
            for (int b = Math.max(0, cy[i] - 1); b <= Math.min(n3[1] - 1, cy[i] + 1); b++)
            for (int c = Math.max(0, cz[i] - 1); c <= Math.min(n3[2] - 1, cz[i] + 1); c++)
                for (int j = head[(a * n3[1] + b) * n3[2] + c]; j != -1; j = nxt[j]) {
                    if (j <= i || (chainOf[i] == chainOf[j] && j - i == 1)) continue;
                    double dx = x[i][0]-x[j][0], dy = x[i][1]-x[j][1], dz = x[i][2]-x[j][2];
                    if (dx*dx + dy*dy + dz*dz < c2) { if (cnt == nbr[i].length) nbr[i] = Arrays.copyOf(nbr[i], cnt * 2); nbr[i][cnt++] = j; }
                }
            nbrCount[i] = cnt;
        }
    }

    static double ljForce(double r, double eps) { double ir = 1.0 / r, s6 = Math.pow(ir, 6); return 24 * eps * (2 * s6 * s6 - s6) * ir; }

    void forces(double cargoZ, double a) {
        for (double[] v : f) { v[0] = v[1] = v[2] = 0; }
        fCargoZ = 0;
        double eps = P.epsScale * EPS0;
        for (int i = 0; i < N; i++) {
            if (idxIn[i] > 0) {
                int j = i - 1; double[] d = {x[i][0]-x[j][0], x[i][1]-x[j][1], x[i][2]-x[j][2]};
                double r = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) + 1e-12, fm = -100 * (r - bond);
                for (int k = 0; k < 3; k++) { f[i][k] += fm * d[k] / r; f[j][k] -= fm * d[k] / r; }
            } else {
                double th = Math.atan2(x[i][1], x[i][0]);
                f[i][0] -= 100 * (x[i][0] - (P.pipeR - 0.5) * Math.cos(th)); f[i][1] -= 100 * (x[i][1] - (P.pipeR - 0.5) * Math.sin(th));
            }
            for (int m = 0; m < nbrCount[i]; m++) {
                int j = nbr[i][m];
                double[] d = {x[i][0]-x[j][0], x[i][1]-x[j][1], x[i][2]-x[j][2]};
                double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2], rr = Math.sqrt(r2), fm;
                if (hps) {
                    if (rr >= RC) continue;
                    double r = Math.max(rr, 0.4), flj = ljForce(r, eps);
                    fm = r < WCA ? flj : 0.5 * (lam[i] + lam[j]) * flj;
                    if (q[i] != 0 && q[j] != 0) { double E = P.lB * q[i] * q[j] * Math.exp(-r / P.xi) / r; fm += E * (1 / r + 1 / P.xi); }
                } else {
                    boolean hh = isH[i] && isH[j]; double rc = hh ? 2.5 : WCA;
                    if (r2 >= rc * rc) continue;
                    fm = ljForce(Math.max(rr, 0.4), hh ? P.epsFG : 1.0);
                }
                for (int k = 0; k < 3; k++) { double g = fm * d[k] / (rr + 1e-12); f[i][k] += g; f[j][k] -= g; }
            }
            double rho = Math.hypot(x[i][0], x[i][1]);
            if (rho > P.pipeR - 0.5) { double fm = -100 * (rho - (P.pipeR - 0.5)); f[i][0] += fm * x[i][0] / rho; f[i][1] += fm * x[i][1] / rho; }
            double zb = P.zMax + 6;
            if (x[i][2] > zb) f[i][2] -= 50 * (x[i][2] - zb); else if (x[i][2] < -zb) f[i][2] -= 50 * (x[i][2] + zb);
            double[] d = {x[i][0], x[i][1], x[i][2] - cargoZ};
            double r = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) + 1e-12, de = r - a, fm = 0;
            if (hps) {
                if (de < RC) {
                    double dd = Math.max(de, 0.4), flj = ljForce(dd, eps);
                    fm = dd < WCA ? flj : 0.5 * (P.cargoLambda + lam[i]) * flj;
                    if (P.cargoCharge != 0 && q[i] != 0) { double E = P.lB * P.cargoCharge * q[i] * Math.exp(-de / P.xi) / (1 + a / P.xi) / r; fm += E * (1 / r + 1 / P.xi); }
                }
            } else {
                boolean att = isH[i] && P.cargoEps > 0; double rc = att ? 2.5 : WCA;
                if (de < rc) fm = ljForce(Math.max(de, 0.4), att ? P.cargoEps : 1.0);
            }
            if (fm != 0) { for (int k = 0; k < 3; k++) f[i][k] += fm * d[k] / r; fCargoZ -= fm * d[2] / r; }
        }
    }

    void step(double cargoZ, double a) {
        forces(cargoZ, a);
        double s = Math.sqrt(2 * P.dt);
        for (int i = 0; i < N; i++) for (int k = 0; k < 3; k++) {
            double dx = f[i][k] * P.dt; if (dx > 0.05) dx = 0.05; else if (dx < -0.05) dx = -0.05;
            x[i][k] += dx + s * rnd.nextGaussian();
        }
    }

    public static Result run(double a, Params p) {
        Result r = new Result(); int W = p.windows;
        r.z = new double[W]; r.force = new double[W]; r.pmf = new double[W];
        if (a >= p.pipeR - 1.0) { r.fits = false; return r; }
        FgBrushSim sim = new FgBrushSim(p);
        int B = 4; double[][] fb = new double[B][W];
        for (int w = 0; w < W; w++) {
            double z = p.zMin + (p.zMax - p.zMin) * w / (W - 1); r.z[w] = z;
            int eq = (w == 0) ? p.equilFirst : p.equil;
            for (int t = 0; t < eq; t++) { if (t % 10 == 0) sim.buildNeighbors(); sim.step(z, a); }
            double[] acc = new double[B]; int per = Math.max(1, p.prod / B);
            for (int t = 0; t < p.prod; t++) {
                if (t % 10 == 0) sim.buildNeighbors();
                sim.step(z, a); acc[Math.min(B - 1, t / per)] += sim.fCargoZ;
            }
            for (int b = 0; b < B; b++) { fb[b][w] = acc[b] / per; r.force[w] += fb[b][w] / B; }
        }
        double dz = (p.zMax - p.zMin) / (W - 1);
        for (int w = 1; w < W; w++) r.pmf[w] = r.pmf[w-1] - 0.5 * (r.force[w] + r.force[w-1]) * dz;
        double mx = 0; for (double v : r.pmf) mx = Math.max(mx, v);
        r.barrier = mx; r.closure = Math.abs(r.pmf[W - 1]) / Math.max(mx, 1.0);
        return r;
    }

    static final double[] T975 = {12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31, 2.26};

    /** Independent replicas in parallel; error bars are across replicas (95% CI via Student t). */
    public static Agg runReplicas(double a, Params p, int replicas) {
        Agg g = new Agg(); g.barriers = new double[replicas]; g.closures = new double[replicas];
        if (a >= p.pipeR - 1.0) { g.fits = false; g.mean = 1e6; return g; }
        ExecutorService ex = Executors.newFixedThreadPool(Math.max(1, Math.min(replicas, Runtime.getRuntime().availableProcessors())));
        List<Future<Result>> fs = new ArrayList<>();
        for (int r = 0; r < replicas; r++) { Params q = p.copy(); q.seed = p.seed + 7919L * r; fs.add(ex.submit(() -> run(a, q))); }
        try {
            for (int r = 0; r < replicas; r++) {
                Result res = fs.get(r).get(); g.barriers[r] = res.barrier; g.closures[r] = res.closure;
                if (res.closure > 0.25) g.badClosure++;
                if (g.meanPmf == null) { g.meanPmf = new double[res.pmf.length]; g.z = res.z; }
                for (int w = 0; w < res.pmf.length; w++) g.meanPmf[w] += res.pmf[w] / replicas;
            }
        } catch (Exception e) { throw new RuntimeException(e); } finally { ex.shutdown(); }
        g.mean = Arrays.stream(g.barriers).average().orElse(0);
        if (replicas > 1) {
            double v = 0; for (double b : g.barriers) v += (b - g.mean) * (b - g.mean) / (replicas - 1);
            g.sd = Math.sqrt(v); g.sem = g.sd / Math.sqrt(replicas);
            g.ci95 = (replicas - 2 < T975.length ? T975[replicas - 2] : 2.0) * g.sem;
        }
        return g;
    }
}
