package npctransit;

import java.nio.file.*;
import java.util.*;

public class Main {
    static final String LABEL = "[MODEL-PREDICTED passive-permeability probability; not experimental evidence]";

    public static void main(String[] a) throws Exception {
        if (a.length == 0) { usage(); return; }
        switch (a[0]) {
            case "features": features(Path.of(a[1])); break;
            case "train": train(Path.of(a[1]), Path.of(a[2]), a.length > 3 ? Double.parseDouble(a[3]) : 0.1); break;
            case "predict": predict(Path.of(a[1]), Path.of(a[2])); break;
            case "dock": dock(Path.of(a[1])); break;
            case "simulate": simulate(Path.of(a[1]), flags(a, 2)); break;
            case "calibrate": calibrate(Path.of(a[1]), flags(a, 2)); break;
            case "regress": regress(Path.of(a[1]), flags(a, 2)); break;
            case "predictq": predictq(Path.of(a[1]), Path.of(a[2])); break;
            case "barrier": barrierJob(Path.of(a[1]), flags(a, 2)); break;
            default: usage();
        }
    }
    static void usage() {
        System.out.println("NPC-Transit 0.4\n"
            + "Analysis commands (validated on a small dataset; see README):\n"
            + "  features <file.pdb>                          structure-derived descriptors\n"
            + "  regress <dataset.csv> [--features stokesA] [--lambda 0.1] [--out model.properties]\n"
            + "                                               ridge regression of ln(q), leave-one-out validation\n"
            + "  predictq <model.properties> <file.pdb>       predicted ln(q) with error band and domain warning\n"
            + "Classification variant: train <dataset.csv> <model.properties> [lambda]; predict <model.properties> <file.pdb>\n"
            + "Experimental modules (NOT validated; development only):\n"
            + "  dock <file.pdb>\n"
            + "  simulate <file.pdb> [flags]   calibrate <dataset.csv> [flags]   barrier <file.pdb> --id X --eps E --cache c.csv\n"
            + "  flags: --fasta f.fasta --replicas N --quality tiny|smoke|quick|standard|production --eps X --pipe R --chains N\n"
            + "         --seglen L --anchorz Z --equilscale K --seed S --calib file --cache c.csv --out report.txt");
    }
    static Features feats(Path pdb) throws Exception { return Features.compute(PdbParser.parse(pdb)); }

    static void features(Path pdb) throws Exception {
        Features f = feats(pdb);
        System.out.printf("residues=%d chains=%d mass=%.1f kDa Rg=%.1f A Stokes~%.1f A hydrophobic=%.3f netCharge=%.0f NLS-like=%b%n",
                f.residues, f.chains, f.massKDa, f.rg, f.stokesApprox, f.hydFrac, f.netCharge, f.nlsLike);
    }

    /** Dataset CSV columns: id,pdb_file,label,source_doi,assay,notes  (label 1 = passes, 0 = excluded; lines starting # ignored). */
    static void train(Path csv, Path out, double lambda) throws Exception {
        List<double[]> X = new ArrayList<>(); List<Integer> Y = new ArrayList<>(); List<String> ids = new ArrayList<>();
        Path base = csv.toAbsolutePath().getParent();
        for (String line : Files.readAllLines(csv)) {
            if (line.isBlank() || line.startsWith("#") || line.startsWith("id,")) continue;
            String[] t = line.split(",", -1);
            if (t[2].isBlank()) { System.err.println("Skipping unlabeled entry: " + t[0]); continue; }
            Path pdb = base.resolve(t[1]);
            if (!Files.exists(pdb)) { System.err.println("Missing PDB file, skipping: " + t[1]); continue; }
            X.add(feats(pdb).vector()); Y.add(Integer.parseInt(t[2].trim())); ids.add(t[0]);
        }
        int n = X.size(), pos = (int) Y.stream().filter(v -> v == 1).count();
        System.out.printf("Usable entries: %d (passes=%d, excluded=%d)%n", n, pos, n - pos);
        if (pos < 5 || n - pos < 5) throw new IllegalStateException("Need >=5 examples per class; refusing to fit.");
        double[][] Xa = X.toArray(new double[0][]); int[] y = Y.stream().mapToInt(Integer::intValue).toArray();
        int k = Math.min(5, Math.min(pos, n - pos));
        double[] cv = Evaluation.crossValidate(Xa, y, k, 20, lambda, 42L);
        double[] ci = Evaluation.bootstrapAucCI(cv, y, 2000, 7L);
        System.out.printf("%d-fold CV x20: AUC=%.3f (95%% bootstrap CI %.3f-%.3f), Brier=%.3f%n", k, Evaluation.auc(cv, y), ci[0], ci[1], Evaluation.brier(cv, y));
        Model m = Model.fit(Xa, y, lambda); m.save(out);
        for (int j = 0; j < Features.NAMES.length; j++) System.out.printf("  coef[%s] = %.3f (standardised)%n", Features.NAMES[j], m.coef[j]);
        System.out.println("Model saved: " + out);
        try (var w = Files.newBufferedWriter(Path.of(out + ".cv_predictions.csv"))) {
            w.write("id,label,cv_probability\n");
            for (int i = 0; i < n; i++) w.write(ids.get(i) + "," + y[i] + "," + String.format("%.4f", cv[i]) + "\n");
        }
    }

    static final double MAX_KDA = 500, BAND = 0.10;

    static void predict(Path modelFile, Path pdb) throws Exception {
        Model m = Model.load(modelFile); Features f = feats(pdb);
        features(pdb);
        if (f.massKDa > MAX_KDA) {
            System.out.printf("DECISION: NOT ASSESSED (%.0f kDa exceeds the %.0f kDa scope of this tool; large assemblies/capsids are out of scope)%n", f.massKDa, MAX_KDA);
            return;
        }
        double p = m.probability(f.vector());
        List<String> ood = m.outOfDomain(f.vector());
        String decision = !ood.isEmpty() ? "NOT ASSESSED (outside training domain: " + ood + ")"
                : p >= 0.5 + BAND ? "YES (likely crosses by passive diffusion)"
                : p <= 0.5 - BAND ? "NO (likely excluded)" : "INDETERMINATE (probability within +/-" + (int) (BAND * 100) + "% of 50%)";
        System.out.printf("P(passes by passive diffusion) = %.1f %%   (trained on n=%d)%n", 100 * p, m.nTrain);
        System.out.println("DECISION: " + decision);
        if (f.nlsLike) System.out.println("NOTE: NLS-like motif present. Facilitated import needs importin binding; not modelled here.");
        if (f.chains > 1) System.out.println("NOTE: multi-chain file; confirm the biological assembly is what you intend.");
        System.out.println(LABEL);
    }

    static void dock(Path pdb) throws Exception {
        PdbParser.Structure s = PdbParser.parse(pdb);
        List<FgProbeDock.Site> sites = FgProbeDock.run(s, 5);
        System.out.println("FG-probe docking (heuristic; relative scores, lower = more FG-pocket-like)");
        if (sites.isEmpty()) System.out.println("No hydrophobic pocket candidates found.");
        int k = 1;
        for (FgProbeDock.Site st : sites)
            System.out.printf("Site %d: score=%.1f hydContacts=%d polarContacts=%d at (%.1f, %.1f, %.1f)  residues: %s%n",
                    k++, st.score, st.hyd, st.polar, st.pos[0], st.pos[1], st.pos[2], String.join(" ", st.residues));
        System.out.println("[heuristic probe docking; unvalidated, not experimental evidence]");
    }

    // ---------------------------------------------------------------- simulation / calibration
    static Map<String, String> flags(String[] a, int from) {
        Map<String, String> m = new HashMap<>();
        for (int i = from; i + 1 < a.length; i += 2) m.put(a[i].replaceFirst("^--", ""), a[i + 1]);
        return m;
    }
    static FgBrushSim.Params params(Map<String, String> o) throws Exception {
        FgBrushSim.Params p = new FgBrushSim.Params();
        String q = o.getOrDefault("quality", "quick");
        switch (q) {
            case "tiny": p.windows = 5; p.equilFirst = 1000; p.equil = 300; p.prod = 600; break;
            case "smoke": p.windows = 7; p.equilFirst = 3000; p.equil = 1000; p.prod = 2000; break;
            case "quick": break;
            case "standard": p.windows = 25; p.equilFirst = 40000; p.equil = 20000; p.prod = 200000; break;
            case "production": p.windows = 37; p.equilFirst = 200000; p.equil = 50000; p.prod = 1000000; break;
            default: throw new IllegalArgumentException("quality must be tiny|smoke|quick|standard|production");
        }
        if (o.containsKey("fasta")) { p.seqs = Hps.readFasta(Path.of(o.get("fasta"))); p.chains = Integer.parseInt(o.getOrDefault("chains", "16")); }
        else if (o.containsKey("chains")) p.chains = Integer.parseInt(o.get("chains"));
        p.pipeR = Double.parseDouble(o.getOrDefault("pipe", "8.0"));
        p.epsFG = Double.parseDouble(o.getOrDefault("eps", "1.5"));
        p.epsScale = Double.parseDouble(o.getOrDefault("eps", "1.0"));
        p.seed = Long.parseLong(o.getOrDefault("seed", "1"));
        if (o.containsKey("seglen") && p.seqs != null) {   // keep only the first L residues (anchored end) of each FG sequence
            int L = Integer.parseInt(o.get("seglen")); List<String> t = new ArrayList<>();
            for (String sq : p.seqs) t.add(sq.length() > L ? sq.substring(0, L) : sq); p.seqs = t;
        }
        if (o.containsKey("anchorz")) p.anchorZ = Double.parseDouble(o.get("anchorz"));
        if (o.containsKey("equilscale")) { double k = Double.parseDouble(o.get("equilscale")); p.equilFirst = (int) (p.equilFirst * k); p.equil = (int) (p.equil * k); }
        if (o.containsKey("calib") && !o.containsKey("eps")) {
            Properties c = new Properties(); try (var in = Files.newInputStream(Path.of(o.get("calib")))) { c.load(in); }
            if (c.containsKey("epsScale")) { p.epsScale = Double.parseDouble(c.getProperty("epsScale")); p.epsFG = 1.5 * p.epsScale; }
        }
        return p;
    }
    /** Map a protein structure onto cargo parameters (sigma = 6 A). */
    static double cargo(Features f, PdbParser.Structure s, FgBrushSim.Params p) {
        p.cargoEps = p.epsFG * Math.max(0, Math.min(1.0, f.hydFrac / 0.40));
        p.cargoLambda = Hps.meanLambda(s.sequence);
        p.cargoCharge = f.netCharge;
        return f.stokesApprox / 6.0;
    }

    static void simulate(Path pdb, Map<String, String> o) throws Exception {
        PdbParser.Structure s = PdbParser.parse(pdb); Features f = Features.compute(s);
        List<String> out = new ArrayList<>();
        if (f.massKDa > MAX_KDA) { System.out.println("NOT SIMULATED: exceeds " + MAX_KDA + " kDa scope."); return; }
        FgBrushSim.Params p = params(o); int reps = Integer.parseInt(o.getOrDefault("replicas", "4"));
        double a = cargo(f, s, p);
        out.add("NPC-Transit 0.3 simulate | " + new Date());
        out.add(String.format("mode=%s quality=%s replicas=%d seed=%d windows=%d prod=%d chains=%d pipeR=%.1f sigma epsScale=%.2f epsFG=%.2f",
                p.seqs == null ? "generic-FXFG" : "HPS(FASTA:" + o.get("fasta") + ")", o.getOrDefault("quality", "quick"), reps, p.seed, p.windows, p.prod, p.chains, p.pipeR, p.epsScale, p.epsFG));
        out.add(String.format("cargo: R=%.2f sigma (%.1f A) lambda=%.3f charge=%.0f eps_c=%.2f", a, f.stokesApprox, p.cargoLambda, p.cargoCharge, p.cargoEps));
        FgBrushSim.Agg g = FgBrushSim.runReplicas(a, p, reps);
        if (!g.fits) out.add("Cargo does not fit the modelled pore (steric exclusion).");
        else {
            out.add("z(sigma)  mean PMF(kT)");
            for (int i = 0; i < g.z.length; i++) out.add(String.format("%7.1f  %9.3f", g.z[i], g.meanPmf[i]));
            out.add(String.format("Barrier = %.2f kT, replica SD %.2f, 95%% CI +/- %.2f (n=%d replicas)", g.mean, g.sd, g.ci95, reps));
            if (reps < 3) out.add("WARNING: <3 replicas; no meaningful uncertainty. Use >=5 for publication.");
            if (g.badClosure > 0) out.add(String.format("WARNING: PMF did not close (return to ~0) in %d/%d replicas. Use a higher quality preset / longer runs; do not interpret.", g.badClosure, reps));
            if (o.containsKey("calib")) {
                Properties c = new Properties(); try (var in = Files.newInputStream(Path.of(o.get("calib")))) { c.load(in); }
                if (c.containsKey("slope")) {
                    double pr = Double.parseDouble(c.getProperty("intercept")) + Double.parseDouble(c.getProperty("slope")) * g.mean, rm = Double.parseDouble(c.getProperty("loo_rmse"));
                    out.add(String.format("Calibrated estimate: ln(q) = %.2f (q = %.3g relative to the reference probe), +/- %.2f ln-units (nested LOO RMSE); epsScale = %s", pr, Math.exp(pr), rm, c.getProperty("epsScale")));
                } else {
                double T = Double.parseDouble(c.getProperty("threshold_kT"));
                out.add(String.format("Calibrated decision (threshold %.2f kT, epsScale from calibration file): %s", T,
                        g.mean + g.ci95 < T ? "PASSES (barrier below threshold, CI included)" : g.mean - g.ci95 > T ? "EXCLUDED" : "INDETERMINATE (CI overlaps threshold)"));
                }
            } else out.add("No --calib file: barrier is NOT interpretable as pass/exclude.");
        }
        out.add("[coarse-grained reduced-unit model; interpret only after calibration + external validation on experimental reference cargo]");
        for (String l : out) System.out.println(l);
        if (o.containsKey("out")) Files.write(Path.of(o.get("out")), out);
    }

    static final Object CACHE_LOCK = new Object();
    /** Barrier for (protein, epsScale), read from / appended to a CSV cache (acts as checkpointing). */
    static double[] cachedBarrier(String id, Path pdb, double eps, FgBrushSim.Params base, int reps, Map<String, String> o) throws Exception {
        String tag = String.format("%s;r%d;c%d;p%.1f;L%s;Z%s;E%s;%s;s%d", o.getOrDefault("quality", "quick"), reps, base.chains, base.pipeR, o.getOrDefault("seglen", "full"), o.getOrDefault("anchorz", "3"), o.getOrDefault("equilscale", "1"), o.containsKey("fasta") ? Path.of(o.get("fasta")).getFileName() : "generic", base.seed);
        Path cache = o.containsKey("cache") ? Path.of(o.get("cache")) : null;
        String prefix = id + "," + eps + "," + tag + ",";
        if (cache != null && Files.exists(cache))
            for (String l : Files.readAllLines(cache)) if (l.startsWith(prefix)) { String[] t = l.substring(prefix.length()).split(","); return new double[]{Double.parseDouble(t[0]), Double.parseDouble(t[1])}; }
        PdbParser.Structure s = PdbParser.parse(pdb); Features f = Features.compute(s);
        FgBrushSim.Params p = base.copy(); p.epsScale = eps; p.epsFG = 1.5 * eps;
        double a = cargo(f, s, p); FgBrushSim.Agg g = FgBrushSim.runReplicas(a, p, reps);
        if (cache != null) synchronized (CACHE_LOCK) {
            Files.writeString(cache, prefix + g.mean + "," + g.sd + "," + g.badClosure + "\n", java.nio.file.StandardOpenOption.CREATE, java.nio.file.StandardOpenOption.APPEND);
        }
        return new double[]{g.mean, g.sd};
    }
    static void barrierJob(Path pdb, Map<String, String> o) throws Exception {
        double[] r = cachedBarrier(o.getOrDefault("id", pdb.getFileName().toString()), pdb, Double.parseDouble(o.getOrDefault("eps", "1.0")), params(o), Integer.parseInt(o.getOrDefault("replicas", "5")), o);
        System.out.printf("barrier=%.3f kT sd=%.3f%n", r[0], r[1]);
    }

    static double[] youden(double[] bar, int[] y) {   // returns {threshold, J}; pass if barrier <= T
        double[] s = bar.clone(); Arrays.sort(s); double bestT = s[0] - 1, bestJ = -2;
        List<Double> cands = new ArrayList<>(); cands.add(s[0] - 1);
        for (int i = 0; i + 1 < s.length; i++) cands.add(0.5 * (s[i] + s[i + 1])); cands.add(s[s.length - 1] + 1);
        int P = 0, Nn = 0; for (int v : y) if (v == 1) P++; else Nn++;
        for (double T : cands) {
            int tp = 0, tn = 0;
            for (int i = 0; i < bar.length; i++) { boolean pass = bar[i] <= T; if (pass && y[i] == 1) tp++; if (!pass && y[i] == 0) tn++; }
            double J = (double) tp / P + (double) tn / Nn - 1; if (J > bestJ) { bestJ = J; bestT = T; }
        }
        return new double[]{bestT, bestJ};
    }

    /** Fits epsScale (grid) and a barrier threshold on experimental labels; nested leave-one-out for honest accuracy. */
    static void calibrate(Path csv, Map<String, String> o) throws Exception {
        if (isRegression(csv)) { calibrateRegression(csv, o); return; }
        FgBrushSim.Params base = params(o); int reps = Integer.parseInt(o.getOrDefault("replicas", "4"));
        double[] grid = Arrays.stream(o.getOrDefault("grid", "0.5,0.75,1.0,1.5,2.0").split(",")).mapToDouble(Double::parseDouble).toArray();
        List<String> ids = new ArrayList<>(); List<Integer> Y = new ArrayList<>(); List<double[]> row = new ArrayList<>();
        Path baseDir = csv.toAbsolutePath().getParent();
        for (String line : Files.readAllLines(csv)) {
            if (line.isBlank() || line.startsWith("#") || line.startsWith("id,")) continue;
            String[] t = line.split(",", -1);
            if (t[2].isBlank() || !Files.exists(baseDir.resolve(t[1]))) { System.err.println("Skipping " + t[0]); continue; }
            if (Features.compute(PdbParser.parse(baseDir.resolve(t[1]))).massKDa > MAX_KDA) { System.err.println("Skipping (scope) " + t[0]); continue; }
            double[] b = new double[grid.length];
            for (int g = 0; g < grid.length; g++) b[g] = cachedBarrier(t[0], baseDir.resolve(t[1]), grid[g], base, reps, o)[0];
            ids.add(t[0]); Y.add(Integer.parseInt(t[2].trim())); row.add(b);
            System.out.printf("%s: barriers over epsScale grid = %s%n", t[0], Arrays.toString(b));
        }
        int n = ids.size(), pos = (int) Y.stream().filter(v -> v == 1).count();
        if (pos < 5 || n - pos < 5) throw new IllegalStateException("Need >=5 examples per class.");
        int[] y = Y.stream().mapToInt(Integer::intValue).toArray();
        double[] auc = new double[grid.length];
        for (int g = 0; g < grid.length; g++) { double[] sc = new double[n]; for (int i = 0; i < n; i++) sc[i] = -row.get(i)[g]; auc[g] = Evaluation.auc(sc, y); }
        int best = 0; for (int g = 1; g < grid.length; g++) if (auc[g] > auc[best]) best = g;
        double[] bb = new double[n]; for (int i = 0; i < n; i++) bb[i] = row.get(i)[best];
        double[] yt = youden(bb, y);
        int tp = 0, tn = 0;   // nested LOO: choose grid value + threshold without the held-out protein
        for (int h = 0; h < n; h++) {
            int bg = 0; double ba = -1;
            for (int g = 0; g < grid.length; g++) {
                List<Double> sc = new ArrayList<>(); List<Integer> yy = new ArrayList<>();
                for (int i = 0; i < n; i++) if (i != h) { sc.add(-row.get(i)[g]); yy.add(y[i]); }
                double v = Evaluation.auc(sc.stream().mapToDouble(Double::doubleValue).toArray(), yy.stream().mapToInt(Integer::intValue).toArray());
                if (v > ba) { ba = v; bg = g; }
            }
            List<Double> bs = new ArrayList<>(); List<Integer> yy = new ArrayList<>();
            for (int i = 0; i < n; i++) if (i != h) { bs.add(row.get(i)[bg]); yy.add(y[i]); }
            double T = youden(bs.stream().mapToDouble(Double::doubleValue).toArray(), yy.stream().mapToInt(Integer::intValue).toArray())[0];
            boolean pass = row.get(h)[bg] <= T; if (pass && y[h] == 1) tp++; if (!pass && y[h] == 0) tn++;
        }
        List<String> out = new ArrayList<>();
        out.add("NPC-Transit 0.3 calibrate | " + new Date() + String.format(" | n=%d (pass=%d, excluded=%d) replicas=%d quality=%s", n, pos, n - pos, reps, o.getOrDefault("quality", "quick")));
        for (int g = 0; g < grid.length; g++) out.add(String.format("epsScale=%.2f  AUC(-barrier)=%.3f", grid[g], auc[g]));
        out.add(String.format("Selected epsScale=%.2f, threshold=%.2f kT (in-sample Youden J=%.2f; optimistic)", grid[best], yt[0], yt[1]));
        out.add(String.format("Nested leave-one-out: sensitivity=%.2f (%d/%d), specificity=%.2f (%d/%d), accuracy=%.2f  <- report THIS, not in-sample", (double) tp / pos, tp, pos, (double) tn / (n - pos), tn, n - pos, (double) (tp + tn) / n));
        out.add("Also validate on an external set (different lab/study). Small n => wide uncertainty.");
        for (String l : out) System.out.println(l);
        Properties pr = new Properties();
        pr.setProperty("epsScale", String.valueOf(grid[best])); pr.setProperty("threshold_kT", String.valueOf(yt[0]));
        pr.setProperty("n", String.valueOf(n)); pr.setProperty("quality", o.getOrDefault("quality", "quick"));
        String outp = o.getOrDefault("out", "calibration.properties");
        try (var w = Files.newOutputStream(Path.of(outp))) { pr.store(w, "NPC-Transit calibration"); }
        Files.write(Path.of(outp + ".report.txt"), out);
    }

    // ---------------------------------------------------------------- regression mode
    static class RegData { List<String> ids = new ArrayList<>(); List<Path> pdb = new ArrayList<>(); List<Double> y = new ArrayList<>(), rexp = new ArrayList<>(); }
    static int idx(String[] h, String n) { for (int i = 0; i < h.length; i++) if (h[i].trim().equals(n)) return i; return -1; }
    static boolean isRegression(Path csv) throws Exception {
        for (String l : Files.readAllLines(csv)) { if (l.isBlank() || l.startsWith("#")) continue; return idx(l.split(",", -1), "y") >= 0; }
        return false;
    }
    /** Regression dataset columns: id,pdb_file,y,... (y = ln of the measured flux ratio q); optional stokes_nm. */
    static RegData loadReg(Path csv) throws Exception {
        RegData d = new RegData(); Path base = csv.toAbsolutePath().getParent(); String[] hdr = null;
        for (String line : Files.readAllLines(csv)) {
            if (line.isBlank() || line.startsWith("#")) continue;
            String[] t = line.split(",", -1);
            if (hdr == null) { hdr = t; continue; }
            int iy = idx(hdr, "y"), ir = idx(hdr, "stokes_nm");
            if (t[iy].isBlank()) continue;
            Path p = base.resolve(t[1]);
            if (!Files.exists(p)) { System.err.println("Missing PDB file, skipping: " + t[0] + " (" + t[1] + ")"); continue; }
            d.ids.add(t[0]); d.pdb.add(p); d.y.add(Double.parseDouble(t[iy].trim()));
            d.rexp.add(ir >= 0 && !t[ir].isBlank() ? Double.parseDouble(t[ir].trim()) : Double.NaN);
        }
        return d;
    }
    static void metrics(String name, double[] p, double[] y) {
        double[] cr = Regression.bootstrapCI(p, y, 2000, 11L, Regression::rmse), cs = Regression.bootstrapCI(p, y, 2000, 12L, Regression::spearman);
        System.out.printf("%-26s RMSE=%.3f (95%% CI %.3f-%.3f) MAE=%.3f R2=%.3f Pearson=%.3f Spearman=%.3f (95%% CI %.3f-%.3f)%n", name,
                Regression.rmse(p, y), cr[0], cr[1], Regression.mae(p, y), Regression.r2(p, y), Regression.pearson(p, y), Regression.spearman(p, y), cs[0], cs[1]);
    }
    static void regress(Path csv, Map<String, String> o) throws Exception {
        RegData d = loadReg(csv); int n = d.ids.size();
        if (n < 6) throw new IllegalStateException("Need >= 6 proteins with structures (have " + n + ").");
        String fs = o.getOrDefault("features", "stokesA"); String[] feats = fs.isEmpty() ? new String[0] : fs.split(",");
        double lambda = Double.parseDouble(o.getOrDefault("lambda", "0.1"));
        double[][] X = new double[n][feats.length]; double[] y = new double[n], sa = new double[n], re = new double[n]; boolean haveR = true;
        System.out.println("id, y=ln(q), " + String.join(",", feats) + ", Stokes~(A)");
        for (int i = 0; i < n; i++) {
            Features f = feats(d.pdb.get(i)); y[i] = d.y.get(i); sa[i] = f.stokesApprox; re[i] = d.rexp.get(i) * 10.0; if (Double.isNaN(re[i])) haveR = false;
            StringBuilder sb = new StringBuilder(); for (int j = 0; j < feats.length; j++) { X[i][j] = f.get(feats[j]); sb.append(String.format(" %.3f", X[i][j])); }
            System.out.printf("%s, %.3f,%s, %.1f%n", d.ids.get(i), y[i], sb, f.stokesApprox);
        }
        if (haveR) System.out.printf("Structure-derived Stokes radius vs measured: Pearson=%.3f, mean abs diff=%.1f A%n", Regression.pearson(sa, re), Regression.mae(sa, re));
        double[] p = Regression.loo(X, y, lambda, feats), p0 = Regression.loo(new double[n][0], y, lambda, new String[0]);
        System.out.printf("n=%d proteins, features=[%s], lambda=%.2f, leave-one-out:%n", n, fs, lambda);
        metrics("model", p, y); metrics("baseline (mean of others)", p0, y);
        Regression.Fit fit = Regression.fit(X, y, lambda, feats);
        String outp = o.getOrDefault("out", "regression.properties"); Properties pr = new Properties();
        pr.setProperty("features", fs); pr.setProperty("intercept", String.valueOf(fit.intercept)); pr.setProperty("n", String.valueOf(n));
        pr.setProperty("mean", Model.join(fit.mean)); pr.setProperty("sd", Model.join(fit.sd)); pr.setProperty("coef", Model.join(fit.coef));
        pr.setProperty("lo", Model.join(fit.lo)); pr.setProperty("hi", Model.join(fit.hi)); pr.setProperty("loo_rmse", String.valueOf(Regression.rmse(p, y)));
        try (var w = Files.newOutputStream(Path.of(outp))) { pr.store(w, "NPC-Transit regression model"); }
        try (var w = Files.newBufferedWriter(Path.of(outp + ".loo_predictions.csv"))) {
            w.write("id,y_true,y_loo_pred\n"); for (int i = 0; i < n; i++) w.write(d.ids.get(i) + "," + y[i] + "," + String.format("%.4f", p[i]) + "\n");
        }
        System.out.println("Model saved: " + outp + "  (pre-specify the feature set before looking at results; do not tune it on LOO)");
        System.out.println("[model-predicted ln(q) relative to the assay's reference probe; not experimental evidence]");
    }
    static void predictq(Path modelFile, Path pdb) throws Exception {
        Properties pr = new Properties(); try (var in = Files.newInputStream(modelFile)) { pr.load(in); }
        String fs = pr.getProperty("features"); String[] feats = fs.isEmpty() ? new String[0] : fs.split(",");
        Regression.Fit f = new Regression.Fit(); f.feats = feats; f.intercept = Double.parseDouble(pr.getProperty("intercept"));
        f.mean = feats.length == 0 ? new double[0] : Model.split(pr.getProperty("mean")); f.sd = feats.length == 0 ? new double[0] : Model.split(pr.getProperty("sd"));
        f.coef = feats.length == 0 ? new double[0] : Model.split(pr.getProperty("coef")); f.lo = feats.length == 0 ? new double[0] : Model.split(pr.getProperty("lo")); f.hi = feats.length == 0 ? new double[0] : Model.split(pr.getProperty("hi"));
        Features ft = feats(pdb); features(pdb);
        if (ft.massKDa > MAX_KDA) { System.out.println("NOT ASSESSED: exceeds the " + MAX_KDA + " kDa scope."); return; }
        double[] x = new double[feats.length]; for (int j = 0; j < feats.length; j++) x[j] = ft.get(feats[j]);
        double pred = f.predict(x), rm = Double.parseDouble(pr.getProperty("loo_rmse"));
        System.out.printf("Predicted ln(q) = %.2f  (q = %.3g relative to the reference probe; x%.2f to x%.2f at +/-1 LOO RMSE)%n", pred, Math.exp(pred), Math.exp(pred - rm), Math.exp(pred + rm));
        List<String> ood = f.outOfDomain(x); if (!ood.isEmpty()) System.out.println("WARNING outside training domain for: " + ood + " -> unreliable");
        if (ft.chains > 1) System.out.println("NOTE: multi-chain file; confirm the biological assembly.");
        System.out.println("[model-predicted; not experimental evidence]");
    }
    static void calibrateRegression(Path csv, Map<String, String> o) throws Exception {
        FgBrushSim.Params base = params(o); int reps = Integer.parseInt(o.getOrDefault("replicas", "4"));
        double[] grid = {0.5, 0.75, 1.0, 1.5, 2.0}; RegData d = loadReg(csv);
        List<String> ids = new ArrayList<>(); List<Double> Y = new ArrayList<>(); List<double[]> row = new ArrayList<>();
        for (int i = 0; i < d.ids.size(); i++) {
            if (feats(d.pdb.get(i)).massKDa > MAX_KDA) continue;
            double[] b = new double[grid.length];
            for (int g = 0; g < grid.length; g++) b[g] = cachedBarrier(d.ids.get(i), d.pdb.get(i), grid[g], base, reps, o)[0];
            ids.add(d.ids.get(i)); Y.add(d.y.get(i)); row.add(b);
        }
        int n = ids.size(); if (n < 6) throw new IllegalStateException("Need >= 6 proteins.");
        double[] y = Y.stream().mapToDouble(Double::doubleValue).toArray(); double[] rho = new double[grid.length]; int best = 0;
        for (int g = 0; g < grid.length; g++) { double[] bb = new double[n]; for (int i = 0; i < n; i++) bb[i] = row.get(i)[g]; rho[g] = Regression.spearman(bb, y); if (rho[g] < rho[best]) best = g; }
        double[] bb = new double[n]; for (int i = 0; i < n; i++) bb[i] = row.get(i)[best]; double[] ab = Regression.linfit(bb, y);
        double[] p = new double[n];
        for (int h = 0; h < n; h++) {
            int bg = 0; double br = 9;
            for (int g = 0; g < grid.length; g++) { List<Double> xs = new ArrayList<>(), ys = new ArrayList<>(); for (int i = 0; i < n; i++) if (i != h) { xs.add(row.get(i)[g]); ys.add(y[i]); }
                double r = Regression.spearman(xs.stream().mapToDouble(Double::doubleValue).toArray(), ys.stream().mapToDouble(Double::doubleValue).toArray()); if (r < br) { br = r; bg = g; } }
            List<Double> xs = new ArrayList<>(), ys = new ArrayList<>(); for (int i = 0; i < n; i++) if (i != h) { xs.add(row.get(i)[bg]); ys.add(y[i]); }
            double[] l = Regression.linfit(xs.stream().mapToDouble(Double::doubleValue).toArray(), ys.stream().mapToDouble(Double::doubleValue).toArray()); p[h] = l[0] + l[1] * row.get(h)[bg];
        }
        List<String> out = new ArrayList<>();
        out.add("NPC-Transit 0.4 calibrate (regression) | " + new Date() + String.format(" | n=%d replicas=%d quality=%s", n, reps, o.getOrDefault("quality", "quick")));
        for (int g = 0; g < grid.length; g++) out.add(String.format("epsScale=%.2f  Spearman(barrier, ln q)=%.3f  (expect negative)", grid[g], rho[g]));
        out.add(String.format("Selected epsScale=%.2f; ln q = %.3f + %.4f * barrier (in-sample, optimistic)", grid[best], ab[0], ab[1]));
        double rm = Regression.rmse(p, y);
        out.add(String.format("Nested leave-one-out: RMSE=%.3f ln-units, MAE=%.3f, R2=%.3f, Pearson=%.3f, Spearman=%.3f  <- report THIS", rm, Regression.mae(p, y), Regression.r2(p, y), Regression.pearson(p, y), Regression.spearman(p, y)));
        out.add("Validate on an external set (e.g. Timney et al. yeast in vivo data) with the frozen calibration.");
        for (String l : out) System.out.println(l);
        Properties pr = new Properties(); pr.setProperty("epsScale", String.valueOf(grid[best])); pr.setProperty("intercept", String.valueOf(ab[0])); pr.setProperty("slope", String.valueOf(ab[1]));
        pr.setProperty("loo_rmse", String.valueOf(rm)); pr.setProperty("n", String.valueOf(n)); pr.setProperty("quality", o.getOrDefault("quality", "quick"));
        String outp = o.getOrDefault("out", "calibration.properties");
        try (var w = Files.newOutputStream(Path.of(outp))) { pr.store(w, "NPC-Transit regression calibration"); }
        Files.write(Path.of(outp + ".report.txt"), out);
        try (var w = Files.newBufferedWriter(Path.of(outp + ".loo_predictions.csv"))) { w.write("id,y_true,y_loo_pred\n"); for (int i = 0; i < n; i++) w.write(ids.get(i) + "," + y[i] + "," + String.format("%.4f", p[i]) + "\n"); }
    }
}
