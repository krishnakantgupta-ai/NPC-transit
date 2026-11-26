package npctransit;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

/** Minimal PDB parser: first model, standard residues, altLoc A/blank only, heavy atoms + CA trace. */
public class PdbParser {
    public static final Map<String, Double> MASS = new HashMap<>();
    public static final Map<String, Character> ONE = new HashMap<>();
    static {
        String[] r = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
        char[] o = "ARNDCQEGHILKMFPSTWYV".toCharArray();
        double[] m = {71.08,156.19,114.10,115.09,103.14,128.13,129.12,57.05,137.14,113.16,113.16,128.17,131.19,147.18,97.12,87.08,101.10,186.21,163.18,99.13};
        for (int i = 0; i < r.length; i++) { ONE.put(r[i], o[i]); MASS.put(r[i], m[i]); }
    }

    public static class Structure {
        public final List<double[]> ca = new ArrayList<>();
        public final List<String> resNames = new ArrayList<>();
        public final StringBuilder sequence = new StringBuilder();
        public int chains;
        // heavy-atom data (for FG-probe docking)
        public final List<double[]> xyz = new ArrayList<>();
        public final List<String> elem = new ArrayList<>(), atomName = new ArrayList<>(), resLabel = new ArrayList<>();
    }

    public static Structure parse(Path file) throws IOException {
        Structure s = new Structure();
        Set<String> seen = new HashSet<>(), chainSet = new HashSet<>();
        for (String line : Files.readAllLines(file)) {
            if (line.startsWith("ENDMDL")) break;
            if (!line.startsWith("ATOM") || line.length() < 54) continue;
            char alt = line.charAt(16);
            if (alt != ' ' && alt != 'A') continue;
            String res = line.substring(17, 20).trim();
            if (!MASS.containsKey(res)) continue;
            String name = line.substring(12, 16).trim();
            String el = line.length() >= 78 ? line.substring(76, 78).trim() : "";
            if (el.isEmpty()) el = name.replaceAll("[^A-Za-z]", "").substring(0, 1);
            el = el.toUpperCase();
            if (el.equals("H") || el.equals("D")) continue;
            double[] p = {Double.parseDouble(line.substring(30, 38).trim()),
                          Double.parseDouble(line.substring(38, 46).trim()),
                          Double.parseDouble(line.substring(46, 54).trim())};
            s.xyz.add(p); s.elem.add(el); s.atomName.add(name);
            s.resLabel.add(res + line.charAt(21) + line.substring(22, 26).trim());
            if (!name.equals("CA")) continue;
            String key = line.charAt(21) + line.substring(22, 27);
            if (!seen.add(key)) continue;
            chainSet.add(String.valueOf(line.charAt(21)));
            s.ca.add(p); s.resNames.add(res); s.sequence.append(ONE.get(res));
        }
        s.chains = chainSet.size();
        if (s.ca.isEmpty()) throw new IOException("No standard residues parsed from " + file);
        return s;
    }
}
