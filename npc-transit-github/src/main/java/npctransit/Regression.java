package npctransit;

import java.util.*;
import java.util.function.ToDoubleBiFunction;

/** Ridge regression on standardised features, leave-one-out validation and summary statistics. */
public class Regression {
    public static double[] solve(double[][] A, double[] b) {
        int n = b.length; double[][] M = new double[n][n + 1];
        for (int i = 0; i < n; i++) { System.arraycopy(A[i], 0, M[i], 0, n); M[i][n] = b[i]; }
        for (int c = 0; c < n; c++) {
            int p = c; for (int r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
            double[] t = M[c]; M[c] = M[p]; M[p] = t; double d = M[c][c]; if (Math.abs(d) < 1e-12) d = 1e-12;
            for (int r = c + 1; r < n; r++) { double f = M[r][c] / d; for (int k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) { double s = M[i][n]; for (int k = i + 1; k < n; k++) s -= M[i][k] * x[k]; double d = M[i][i]; if (Math.abs(d) < 1e-12) d = 1e-12; x[i] = s / d; }
        return x;
    }

    public static class Fit {
        public double[] mean, sd, coef, lo, hi; public double intercept, lambda; public int n; public String[] feats;
        public double predict(double[] x) { double z = intercept; for (int j = 0; j < coef.length; j++) z += coef[j] * (x[j] - mean[j]) / sd[j]; return z; }
        public List<String> outOfDomain(double[] x) { List<String> o = new ArrayList<>(); for (int j = 0; j < x.length; j++) if (x[j] < lo[j] || x[j] > hi[j]) o.add(feats[j]); return o; }
    }

    public static Fit fit(double[][] X, double[] y, double lambda, String[] feats) {
        int n = X.length, d = feats.length; Fit f = new Fit(); f.n = n; f.lambda = lambda; f.feats = feats;
        f.mean = new double[d]; f.sd = new double[d]; f.lo = new double[d]; f.hi = new double[d];
        Arrays.fill(f.lo, Double.MAX_VALUE); Arrays.fill(f.hi, -Double.MAX_VALUE);
        for (double[] r : X) for (int j = 0; j < d; j++) { f.mean[j] += r[j] / n; f.lo[j] = Math.min(f.lo[j], r[j]); f.hi[j] = Math.max(f.hi[j], r[j]); }
        for (double[] r : X) for (int j = 0; j < d; j++) f.sd[j] += Math.pow(r[j] - f.mean[j], 2) / n;
        for (int j = 0; j < d; j++) { f.sd[j] = Math.sqrt(f.sd[j]); if (f.sd[j] < 1e-9) f.sd[j] = 1; }
        double ybar = Arrays.stream(y).average().orElse(0);
        f.intercept = ybar; f.coef = new double[d];
        if (d > 0) {
            double[][] A = new double[d][d]; double[] b = new double[d];
            for (int i = 0; i < n; i++) for (int j = 0; j < d; j++) {
                double zj = (X[i][j] - f.mean[j]) / f.sd[j]; b[j] += zj * (y[i] - ybar) / n;
                for (int k = 0; k < d; k++) A[j][k] += zj * (X[i][k] - f.mean[k]) / f.sd[k] / n;
            }
            for (int j = 0; j < d; j++) A[j][j] += lambda;
            f.coef = solve(A, b);
        }
        return f;
    }

    public static double[] loo(double[][] X, double[] y, double lambda, String[] feats) {
        int n = X.length; double[] p = new double[n];
        for (int h = 0; h < n; h++) {
            double[][] Xt = new double[n - 1][]; double[] yt = new double[n - 1]; int c = 0;
            for (int i = 0; i < n; i++) if (i != h) { Xt[c] = X[i]; yt[c++] = y[i]; }
            p[h] = fit(Xt, yt, lambda, feats).predict(X[h]);
        }
        return p;
    }

    public static double rmse(double[] p, double[] y) { double s = 0; for (int i = 0; i < p.length; i++) s += Math.pow(p[i] - y[i], 2); return Math.sqrt(s / p.length); }
    public static double mae(double[] p, double[] y) { double s = 0; for (int i = 0; i < p.length; i++) s += Math.abs(p[i] - y[i]); return s / p.length; }
    public static double r2(double[] p, double[] y) {
        double m = Arrays.stream(y).average().orElse(0), a = 0, b = 0;
        for (int i = 0; i < p.length; i++) { a += Math.pow(y[i] - p[i], 2); b += Math.pow(y[i] - m, 2); }
        return 1 - a / b;
    }
    public static double pearson(double[] a, double[] b) {
        int n = a.length; double ma = Arrays.stream(a).average().orElse(0), mb = Arrays.stream(b).average().orElse(0), sab = 0, saa = 0, sbb = 0;
        for (int i = 0; i < n; i++) { sab += (a[i] - ma) * (b[i] - mb); saa += Math.pow(a[i] - ma, 2); sbb += Math.pow(b[i] - mb, 2); }
        return sab / Math.sqrt(saa * sbb + 1e-300);
    }
    static double[] ranks(double[] a) {
        int n = a.length; Integer[] idx = new Integer[n]; for (int i = 0; i < n; i++) idx[i] = i;
        Arrays.sort(idx, Comparator.comparingDouble(i -> a[i])); double[] r = new double[n];
        for (int i = 0; i < n; ) { int j = i; while (j + 1 < n && a[idx[j + 1]] == a[idx[i]]) j++; for (int k = i; k <= j; k++) r[idx[k]] = (i + j) / 2.0 + 1; i = j + 1; }
        return r;
    }
    public static double spearman(double[] a, double[] b) { return pearson(ranks(a), ranks(b)); }
    public static double[] linfit(double[] x, double[] y) {
        double mx = Arrays.stream(x).average().orElse(0), my = Arrays.stream(y).average().orElse(0), sxy = 0, sxx = 0;
        for (int i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += Math.pow(x[i] - mx, 2); }
        double b = sxy / (sxx + 1e-300); return new double[]{my - b * mx, b};
    }
    public static double[] bootstrapCI(double[] p, double[] y, int B, long seed, ToDoubleBiFunction<double[], double[]> m) {
        Random rnd = new Random(seed); int n = p.length; List<Double> v = new ArrayList<>();
        for (int b = 0; b < B; b++) {
            double[] pp = new double[n], yy = new double[n]; for (int i = 0; i < n; i++) { int j = rnd.nextInt(n); pp[i] = p[j]; yy[i] = y[j]; }
            double x = m.applyAsDouble(pp, yy); if (!Double.isNaN(x)) v.add(x);
        }
        Collections.sort(v); if (v.isEmpty()) return new double[]{Double.NaN, Double.NaN};
        return new double[]{v.get((int) (0.025 * v.size())), v.get(Math.min(v.size() - 1, (int) (0.975 * v.size())))};
    }
}
