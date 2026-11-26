package npctransit;

import java.util.*;

/** Stratified repeated k-fold cross-validation, AUC, Brier score, bootstrap CI for AUC. */
public class Evaluation {
    public static double auc(double[] p, int[] y) {
        double wins = 0; long pairs = 0;
        for (int i = 0; i < p.length; i++) if (y[i] == 1) for (int j = 0; j < p.length; j++) if (y[j] == 0) {
            pairs++; wins += p[i] > p[j] ? 1 : (p[i] == p[j] ? 0.5 : 0);
        }
        return pairs == 0 ? Double.NaN : wins / pairs;
    }
    public static double brier(double[] p, int[] y) {
        double s = 0; for (int i = 0; i < p.length; i++) s += Math.pow(p[i] - y[i], 2); return s / p.length;
    }
    /** Out-of-fold predictions, averaged over repeats. */
    public static double[] crossValidate(double[][] X, int[] y, int k, int repeats, double lambda, long seed) {
        int n = X.length; double[] pred = new double[n];
        Random rnd = new Random(seed);
        for (int r = 0; r < repeats; r++) {
            int[] fold = new int[n];
            for (int cls = 0; cls <= 1; cls++) {
                List<Integer> idx = new ArrayList<>();
                for (int i = 0; i < n; i++) if (y[i] == cls) idx.add(i);
                Collections.shuffle(idx, rnd);
                for (int t = 0; t < idx.size(); t++) fold[idx.get(t)] = t % k;
            }
            for (int f = 0; f < k; f++) {
                List<Integer> tr = new ArrayList<>(), te = new ArrayList<>();
                for (int i = 0; i < n; i++) (fold[i] == f ? te : tr).add(i);
                double[][] Xt = new double[tr.size()][]; int[] yt = new int[tr.size()];
                for (int i = 0; i < tr.size(); i++) { Xt[i] = X[tr.get(i)]; yt[i] = y[tr.get(i)]; }
                Model m = Model.fit(Xt, yt, lambda);
                for (int i : te) pred[i] += m.probability(X[i]) / repeats;
            }
        }
        return pred;
    }
    public static double[] bootstrapAucCI(double[] p, int[] y, int B, long seed) {
        Random rnd = new Random(seed); int n = p.length; List<Double> a = new ArrayList<>();
        for (int b = 0; b < B; b++) {
            double[] pp = new double[n]; int[] yy = new int[n];
            for (int i = 0; i < n; i++) { int j = rnd.nextInt(n); pp[i] = p[j]; yy[i] = y[j]; }
            double v = auc(pp, yy); if (!Double.isNaN(v)) a.add(v);
        }
        Collections.sort(a);
        if (a.isEmpty()) return new double[]{Double.NaN, Double.NaN};
        return new double[]{a.get((int) (0.025 * a.size())), a.get(Math.min(a.size() - 1, (int) (0.975 * a.size())))};
    }
}
