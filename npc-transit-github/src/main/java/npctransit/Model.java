package npctransit;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/** L2-regularised logistic regression on standardised features, with applicability-domain bounds. */
public class Model {
    public double[] mean, sd, coef, lo, hi;
    public double intercept;
    public int nTrain;

    public double probability(double[] x) {
        double z = intercept;
        for (int i = 0; i < x.length; i++) z += coef[i] * (x[i] - mean[i]) / sd[i];
        return 1.0 / (1.0 + Math.exp(-z));
    }

    /** Names of features lying outside the range seen in training. */
    public List<String> outOfDomain(double[] x) {
        List<String> out = new ArrayList<>();
        for (int i = 0; i < x.length; i++) if (x[i] < lo[i] || x[i] > hi[i]) out.add(Features.NAMES[i]);
        return out;
    }

    public static Model fit(double[][] X, int[] y, double lambda) {
        int n = X.length, d = X[0].length;
        Model m = new Model();
        m.nTrain = n; m.mean = new double[d]; m.sd = new double[d]; m.lo = new double[d]; m.hi = new double[d];
        Arrays.fill(m.lo, Double.MAX_VALUE); Arrays.fill(m.hi, -Double.MAX_VALUE);
        for (double[] r : X) for (int j = 0; j < d; j++) {
            m.mean[j] += r[j] / n; m.lo[j] = Math.min(m.lo[j], r[j]); m.hi[j] = Math.max(m.hi[j], r[j]);
        }
        for (double[] r : X) for (int j = 0; j < d; j++) m.sd[j] += Math.pow(r[j] - m.mean[j], 2) / n;
        for (int j = 0; j < d; j++) { m.sd[j] = Math.sqrt(m.sd[j]); if (m.sd[j] < 1e-9) m.sd[j] = 1; }
        double[][] Z = new double[n][d];
        for (int i = 0; i < n; i++) for (int j = 0; j < d; j++) Z[i][j] = (X[i][j] - m.mean[j]) / m.sd[j];
        double[] w = new double[d]; double b = 0, lr = 0.1;
        for (int it = 0; it < 20000; it++) {
            double[] gw = new double[d]; double gb = 0;
            for (int i = 0; i < n; i++) {
                double z = b; for (int j = 0; j < d; j++) z += w[j] * Z[i][j];
                double e = 1.0 / (1.0 + Math.exp(-z)) - y[i];
                gb += e / n; for (int j = 0; j < d; j++) gw[j] += e * Z[i][j] / n;
            }
            for (int j = 0; j < d; j++) w[j] -= lr * (gw[j] + lambda * w[j]);
            b -= lr * gb;
        }
        m.coef = w; m.intercept = b;
        return m;
    }

    public void save(Path p) throws IOException {
        Properties pr = new Properties();
        pr.setProperty("features", String.join(",", Features.NAMES));
        pr.setProperty("intercept", String.valueOf(intercept));
        pr.setProperty("nTrain", String.valueOf(nTrain));
        pr.setProperty("mean", join(mean)); pr.setProperty("sd", join(sd)); pr.setProperty("coef", join(coef));
        pr.setProperty("lo", join(lo)); pr.setProperty("hi", join(hi));
        try (OutputStream o = Files.newOutputStream(p)) { pr.store(o, "NPC-Transit trained model"); }
    }

    public static Model load(Path p) throws IOException {
        Properties pr = new Properties();
        try (InputStream i = Files.newInputStream(p)) { pr.load(i); }
        if (!pr.getProperty("features").equals(String.join(",", Features.NAMES)))
            throw new IOException("Model feature set does not match this software version");
        Model m = new Model();
        m.intercept = Double.parseDouble(pr.getProperty("intercept"));
        m.nTrain = Integer.parseInt(pr.getProperty("nTrain"));
        m.mean = split(pr.getProperty("mean")); m.sd = split(pr.getProperty("sd")); m.coef = split(pr.getProperty("coef"));
        m.lo = split(pr.getProperty("lo")); m.hi = split(pr.getProperty("hi"));
        return m;
    }
    public static String join(double[] a) { StringBuilder s = new StringBuilder(); for (double v : a) s.append(s.length() > 0 ? "," : "").append(v); return s.toString(); }
    public static double[] split(String s) { String[] t = s.split(","); double[] a = new double[t.length]; for (int i = 0; i < t.length; i++) a[i] = Double.parseDouble(t[i]); return a; }
}
