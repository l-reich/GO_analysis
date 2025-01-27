import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
//import java.math.BigDecimal;
//import java.math.MathContext;

import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
//import org.apache.commons.math3.distribution.HypergeometricDistribution;
//import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

//root is 0008150

public class DAGHandler {
    private OboParser oboParser;
    private MappingParser mappingParser;
    private String root;
    private String o;
    private int minsize;
    private int maxsize;
    private boolean overlap;
    private String overlapOut;
    private HashMap<String, GOEntry> dag;
    private HashMap<String, String> mapping;
    private HashMap<String, Gene> geneMap;
    private HashSet<String> enrichGenes;
    private HashSet<String> mappingGenes;
    private HashSet<String> relevantGenes;
    private HashSet<IntPair> overlapPairs = new HashSet<>(); // Store all pairs
    private HashSet<IdPair> idPairs = new HashSet<>();
    HashSet<IdPair> relPairs = new HashSet<>();
    int significantGenes;
    GOEntry rootEntry;
    HashMap<String, List<String>> geneToGoMap = new HashMap<>();

    public DAGHandler(String obo, String root, String mapping, boolean mapGo, String enrich, String o, int minsize, int maxsize, boolean overlap, String overlapOut) {
        this.oboParser = new OboParser(obo, root);
        this.root = root;
        this.o = o;
        this.minsize = minsize;
        this.maxsize = maxsize;
        this.overlap = overlap;
        this.overlapOut = overlapOut;
        this.dag = oboParser.getDag();
        rootEntry = oboParser.getRootEntry();
        this.geneMap = new HashMap<>();
        this.enrichGenes = new HashSet<>();
        this.mappingGenes = new HashSet<>();
        if (!overlap) {
            parseEnrichFile(enrich);
        }
        if (mapGo) {
            if (overlap) {
                overlapParseGoMappingFile(mapping);
            } else {
                parseGoMappingFile(mapping);
            }
            parseGoMappingFile(mapping);
        } else {
            if (overlap) {
                overlapParseEnsMappingFile(mapping);
            } else {
                parseEnsMappingFile(mapping);
            }
        }
        propagateGenes();

        if (overlap) {
            for (GOEntry goEntry : dag.values()) {
                for (String gene : goEntry.getMappedGenes()) {
                    for (String goID : geneToGoMap.get(gene)) {
                        if (dag.containsKey(goID)) {
                            if (!goID.equals(goEntry.getId())) {
                                relPairs.add(new IdPair(goEntry.getId(), goID));
                            }
                        }
                    }
                }
            }
        }

        filterGenesInGO();
        relevantGenes = new HashSet<>(enrichGenes);
        relevantGenes.retainAll(mappingGenes);
        calculateNoverlap();

//        int validPairsCount = 0;
//
//        for (IdPair pair : relPairs) {
//
//            // Retrieve GOEntry objects for the terms
//            GOEntry entry1 = dag.get(pair.getFirst());
//            GOEntry entry2 = dag.get(pair.getSecond());
//
//            // Check if both entries exist and satisfy minsize and maxsize constraints
//            if (entry1 != null && entry2 != null &&
//                    entry1.getMappedGenes().size() >= minsize && entry1.getMappedGenes().size() <= maxsize &&
//                    entry2.getMappedGenes().size() >= minsize && entry2.getMappedGenes().size() <= maxsize) {
//                validPairsCount++;
//            } else {
//                // Remove the invalid pair
//                // idPairs.remove(pair);
//            }
//        }
//
//// Print the size of the valid pairs
//        System.out.println("j");
//        System.out.println(validPairsCount);

        //calculate the number of significant genes in relevant genes
        for (String geneId : relevantGenes) {
            Gene gene = geneMap.get(geneId);
            if (gene != null && gene.isSignif()) {
                significantGenes++;
            }
        }

//        if (!overlap) {
            if (false){
            calculateHgPvals();
            //newCalculateHgPvals();
            applyBenjaminiHochbergCorrection(minsize, maxsize);
            calculateFejPvals(minsize, maxsize);
            applyFejBenjaminiHochbergCorrection(minsize, maxsize);
            calculateKsStats(minsize, maxsize);
            applyKsBenjaminiHochbergCorrection(minsize, maxsize);


            // Collect all true GO entries in the DAG
            List<GOEntry> trueEntries = new ArrayList<>();
            for (GOEntry entry : dag.values()) {
                if (entry.isTrue()) {
                    trueEntries.add(entry);
                }
            }

            for (GOEntry goEntry : dag.values()) {
                if (goEntry.getMappedGenes().size() >= minsize && goEntry.getMappedGenes().size() <= maxsize) {
                    String shortestPath = newCalculateShortestPathToTrue(goEntry, trueEntries);
                    goEntry.setShortestPathToATrue(shortestPath);
                }
            }

            writeOutput(o, minsize, maxsize);
        }
        if (overlap) {
            writeOverlapOut(overlapOut);
        }
    }

    public void propagateGenes() {
        if (rootEntry != null) {
//            if (overlap) {
//                propagateGenesFromNodeOverlap(rootEntry);
//            } else {
            propagateGenesFromNode(rootEntry);
            //}
        }
    }

    private void propagateGenesFromNodeOverlap(GOEntry entry) {
        // If the entry has children, propagate genes from them first
        for (GOEntry child : entry.getChildren()) {
            propagateGenesFromNode(child);

            // Add all genes from the child to the current entry
            entry.addMappedGenes(child.getMappedGenes());
            entry.addOverlaps(child.getOverlaps());
//            if (child.getOverlaps() != null) {
//                for (String o : child.getOverlaps()) {
//                    idPairs.add(new IdPair(entry.getId(), o));
//                }
//            }
        }
        for (String o : entry.getOverlaps()) {
            idPairs.add(new IdPair(entry.getId(), o));
        }
    }

    private void propagateGenesFromNode(GOEntry entry) {
        // If the entry has children, propagate genes from them first
        for (GOEntry child : entry.getChildren()) {
            propagateGenesFromNode(child);
            // Add all genes from the child to the current entry
            entry.addMappedGenes(child.getMappedGenes());
            if (overlap) {
                for (String gene : child.getMappedGenes()) {
                    geneToGoMap.computeIfAbsent(gene, k -> new ArrayList<>()).add(entry.getId());
                }
            }
        }
    }

    public void filterGenesInGO() {
        for (GOEntry goEntry : dag.values()) {
            // Get the intersection of mapped genes and enrichGenes
            HashSet<String> filteredGenes = new HashSet<>(goEntry.getMappedGenes());
            filteredGenes.retainAll(enrichGenes);

            // Update the GOEntry with the filtered gene set size
            goEntry.setEnrichedGenes(filteredGenes);
        }
    }

    public void calculateNoverlap() {
        for (GOEntry goEntry : dag.values()) {
            int noverlapCount = 0;

            // Get the subset of enrichGenes mapped to this GOEntry
            HashSet<String> gofilteredGenes = new HashSet<>(goEntry.getMappedGenes());
            gofilteredGenes.retainAll(enrichGenes);

            // Count the number of significant genes
            for (String geneId : gofilteredGenes) {
                Gene gene = geneMap.get(geneId); // Retrieve the Gene object
                if (gene != null && gene.isSignif()) {
                    noverlapCount++;
                }
            }

            // Store the noverlap value in the GOEntry
            goEntry.setNoverlap(noverlapCount);
        }
    }

    public void calculateHgPvals() {
        int totalGenes = relevantGenes.size(); // N: Total measured genes
        int totalSignificantGenes = significantGenes; // K: Total significant genes

        for (GOEntry goEntry : dag.values()) {
            int mappedGenes = goEntry.getEnrichedGenes().size(); // n: Genes mapped to the GO entry
            int significantMappedGenes = goEntry.getNoverlap(); // k: Significant genes in the GO entry

            double pValue = computeLogHypergeometricPval(totalGenes, totalSignificantGenes, mappedGenes, significantMappedGenes);
            //}
            goEntry.setHgPval(pValue);
        }
    }

    public void applyBenjaminiHochbergCorrection(int minsize, int maxsize) {
        // Filter GO entries based on size constraints and create a mutable list
        List<GOEntry> filteredEntries = new ArrayList<>(
                dag.values().stream()
                        .filter(goEntry -> goEntry.getMappedGenes().size() >= minsize && goEntry.getMappedGenes().size() <= maxsize)
                        .toList()
        );

        // Sort by raw p-values in ascending order
        filteredEntries.sort(Comparator.comparingDouble(GOEntry::getHgPval));

        int m = filteredEntries.size(); // Total number of valid tests
        double[] correctedPvals = new double[m];

        // Apply Benjamini-Hochberg correction
        for (int i = 0; i < m; i++) {
            double rawPval = filteredEntries.get(i).getHgPval();
            correctedPvals[i] = rawPval * m / (i + 1); // BH formula
        }

        // Ensure corrected p-values are monotonically increasing
        for (int i = m - 2; i >= 0; i--) {
            correctedPvals[i] = Math.min(correctedPvals[i], correctedPvals[i + 1]);
        }

        // Assign corrected p-values back to the filtered GO entries
        for (int i = 0; i < m; i++) {
            filteredEntries.get(i).setHgFdr(Math.min(correctedPvals[i], 1.0)); // Ensure FDR ≤ 1.0
        }
    }

    public void calculateFejPvals(int minsize, int maxsize) {
        int totalGenes = relevantGenes.size(); // N: Total measured genes
        int totalSignificantGenes = significantGenes; // K: Total significant genes

        // Filter GO entries based on size constraints
        List<GOEntry> filteredEntries = dag.values().stream()
                .filter(goEntry -> goEntry.getSize() >= minsize && goEntry.getSize() <= maxsize)
                .toList();

        for (GOEntry goEntry : filteredEntries) {
            int mappedGenes = goEntry.getSize(); // n: Genes mapped to the GO entry
            int significantMappedGenes = goEntry.getNoverlap(); // a: Significant genes in the GO entry

            if (mappedGenes == 0 || significantMappedGenes == 0) {
                goEntry.setFejPval(1.0); // No enrichment possible
                continue;
            }

            double fejPval = computeLogHypergeometricPval(totalGenes - 1, totalSignificantGenes - 1, mappedGenes - 1, significantMappedGenes - 1);

            goEntry.setFejPval(fejPval);
        }
    }


    public void applyFejBenjaminiHochbergCorrection(int minsize, int maxsize) {
        // Filter GO entries based on size constraints and create a mutable list
        List<GOEntry> filteredEntries = new ArrayList<>(
                dag.values().stream()
                        .filter(goEntry -> goEntry.getSize() >= minsize && goEntry.getSize() <= maxsize)
                        .toList()
        );

        // Sort by Fisher's exact p-values in ascending order
        filteredEntries.sort(Comparator.comparingDouble(GOEntry::getFejPval));

        int m = filteredEntries.size(); // Total number of valid tests
        double[] correctedPvals = new double[m];

        // Apply Benjamini-Hochberg correction
        for (int i = 0; i < m; i++) {
            double rawPval = filteredEntries.get(i).getFejPval();
            correctedPvals[i] = rawPval * m / (i + 1); // BH formula
        }

        // Ensure corrected p-values are monotonically increasing
        for (int i = m - 2; i >= 0; i--) {
            correctedPvals[i] = Math.min(correctedPvals[i], correctedPvals[i + 1]);
        }

        // Assign corrected p-values back to the filtered GO entries
        for (int i = 0; i < m; i++) {
            filteredEntries.get(i).setFejFdr(Math.min(correctedPvals[i], 1.0)); // Ensure FDR ≤ 1.0
        }
    }

    //parsing
    private void parseEnsMappingFile(String mapping) {
        try (BufferedReader reader = new BufferedReader(new FileReader(mapping))) {
            String line = reader.readLine(); // Skip the header line
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\t");
                if (parts.length >= 3) {
                    String hgncSymbol = parts[1]; // HGNC symbol
                    if (hgncSymbol.isEmpty()) {
                        continue;
                    }
                    String[] goTerms = parts[2].split("\\|"); // Split GO terms

                    boolean isMapped = false;
                    for (String goTerm : goTerms) {
                        // Add the gene to the corresponding GO entry
                        if (dag.containsKey(goTerm)) {
                            dag.get(goTerm).addMappedGene(hgncSymbol);
                            isMapped = true;
                        }
                    }
                    if (isMapped) {
                        mappingGenes.add(hgncSymbol);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void overlapParseEnsMappingFile(String mapping) {
        try (BufferedReader reader = new BufferedReader(new FileReader(mapping))) {
            String line = reader.readLine(); // Skip the header line
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\t");
                if (parts.length >= 3) {
                    String hgncSymbol = parts[1]; // HGNC symbol
                    if (hgncSymbol.isEmpty()) {
                        continue;
                    }
                    String[] goTerms = parts[2].split("\\|"); // Split GO terms

                    boolean isMapped = false;

//                    for (int i = 0; i < goTerms.length; i++) {
//                        if (dag.containsKey(goTerms[i])) {
//                            for (int j = i + 1; j < goTerms.length; j++) {
//                                if (dag.containsKey(goTerms[j])) {
//                                    dag.get(goTerms[i]).addOverlap(goTerms[j]);
//                                    dag.get(goTerms[j]).addOverlap(goTerms[i]);
//                                    idPairs.add(new IdPair(goTerms[i], goTerms[j]));
//                                }
//                            }
//                        }
//                    }

                    for (String goTerm : goTerms) {
                        if (dag.containsKey(goTerm)) {
                            dag.get(goTerm).addMappedGene(hgncSymbol);
                            isMapped = true;
                        }
                    }
                    geneToGoMap.computeIfAbsent(hgncSymbol, k -> new ArrayList<>()).addAll(Arrays.asList(goTerms));


                    if (isMapped) {
                        mappingGenes.add(hgncSymbol);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void parseGoMappingFile(String mapping) {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(mapping))))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("!")) {
                    continue; // Skip comment lines
                }

                String[] parts = line.split("\t");
                if (parts.length >= 5) {
                    String geneName = parts[2];
                    String qualifier = parts[3];
                    String goTerm = parts[4];

                    // Skip the line if the qualifier field is not empty
                    if (!qualifier.isEmpty()) {
                        continue;
                    }

                    // Add the gene to the corresponding GO entry
                    if (dag.containsKey(goTerm)) {
                        dag.get(goTerm).addMappedGene(geneName);
                    }

                    // Keep track of genes mapped
                    mappingGenes.add(geneName);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void overlapParseGoMappingFile(String mapping) {

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(mapping))))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("!")) {
                    continue; // Skip comment lines
                }

                String[] parts = line.split("\t");
                if (parts.length >= 5) {
                    String geneName = parts[2];
                    String qualifier = parts[3];
                    String goTerm = parts[4];

                    // Skip the line if the qualifier field is not empty
                    if (!qualifier.isEmpty()) {
                        continue;
                    }

                    // Add the GO term to the gene's list in the map
                    geneToGoMap.computeIfAbsent(geneName, k -> new ArrayList<>()).add(goTerm);

                    // Add the gene to the corresponding GO entry
                    if (dag.containsKey(goTerm)) {
                        dag.get(goTerm).addMappedGene(geneName);
                    }

                    // Keep track of genes mapped
                    mappingGenes.add(geneName);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Process the gene-to-GO map to calculate overlaps
//        for (List<String> goTerms : geneToGoMap.values()) {
//            // Generate all unique pairs of GO terms for the current gene
//            for (int i = 0; i < goTerms.size(); i++) {
//                if (!dag.containsKey(goTerms.get(i))) continue;
//
//                for (int j = i + 1; j < goTerms.size(); j++) {
//                    if (!dag.containsKey(goTerms.get(j))) continue;
//
//                    // Add the overlap relationship to both entries
//                    dag.get(goTerms.get(i)).addOverlap(goTerms.get(j));
//                    dag.get(goTerms.get(j)).addOverlap(goTerms.get(i));
//
//                    // Add the pair to the global set
//                    idPairs.add(new IdPair(goTerms.get(i), goTerms.get(j)));
//                }
//            }
//        }
    }


    private void parseEnrichFile(String filePath) {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean readingGOEntries = true;

            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if (line.startsWith("#")) {
                    String goId = line.substring(1);
                    if (dag.containsKey(goId)) {
                        dag.get(goId).setTrue(true);
                    }
                } else if (line.startsWith("id")) {
                    readingGOEntries = false;
                } else if (!readingGOEntries && !line.isEmpty()) {
                    String[] parts = line.split("\t");
                    if (parts.length == 3) {
                        String id = parts[0];
                        double foldChange = Double.parseDouble(parts[1]);
                        boolean signif = Boolean.parseBoolean(parts[2]);
                        Gene gene = new Gene(id, foldChange, signif);
                        geneMap.put(id, gene);
                        enrichGenes.add(id);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeOutput(String outputFilePath, int minsize, int maxsize) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            // Write the header
            writer.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true\n");

            // Iterate through GO entries that satisfy the size constraints
            for (GOEntry goEntry : dag.values()) {
                int size = goEntry.getSize();
                if (size >= minsize && size <= maxsize) {
                    // Collect values
                    String term = goEntry.getId();
                    String name = goEntry.getName();
                    boolean isTrue = goEntry.isTrue();
                    int noverlap = goEntry.getNoverlap();
                    double hgPval = goEntry.getHgPval();
                    double hgFdr = goEntry.getHgFdr();
                    double fejPval = goEntry.getFejPval();
                    double fejFdr = goEntry.getFejFdr(); // Placeholder (set actual value if available)
                    double ksStat = goEntry.getKsStat(); // Placeholder (set actual value if available)
                    double ksPval = goEntry.getKsPval(); // Placeholder (set actual value if available)
                    double ksFdr = goEntry.getKsFdr(); // Placeholder (set actual value if available)
                    String shortestPathToATrue = goEntry.getShortestPathToATrue(); // Placeholder (set actual value if available)

//                     Write the row
                    writer.write(String.format(
                            "%s\t%s\t%d\t%s\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5f\t%.5e\t%.5e\t%s\n",
                            term, name, size, isTrue, noverlap, hgPval, hgFdr, fejPval, fejFdr, ksStat, ksPval, ksFdr, shortestPathToATrue
                    ));
//                    writer.write(term+"\t" +shortestPathToATrue+ "\n"
//                    );
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public double computeLogHypergeometricPval(int N, int K, int n, int k) {
        double logDenominator = newLogCombination(N, n); // log(denominator for all probabilities)
        double sum = 0.0;

        // Compute the sum of probabilities for X >= k
        for (int i = k; i <= Math.min(K, n); i++) {
            double logNumerator = newLogCombination(K, i) + newLogCombination(N - K, n - i); // log(numerator)
            double logP = logNumerator - logDenominator; // log(P(X=i))
            sum += Math.exp(logP); // Convert back to normal space and sum
        }

        return sum; // This is P(X >= k)
    }


    // Helper to calculate log(nCr)
    private double logCombination(int n, int r) {
        if (r > n || r < 0) return Double.NEGATIVE_INFINITY;
        double logResult = 0.0;
        for (int i = 1; i <= r; i++) {
            logResult += Math.log(n - i + 1) - Math.log(i);
        }
        return logResult;
    }

    private double newLogCombination(int n, int r) {
        if (r > n || r < 0) return Double.NEGATIVE_INFINITY; // Invalid combination
        double logResult = 0.0;
        for (int i = 1; i <= r; i++) {
            logResult += Math.log(n - i + 1) - Math.log(i);
        }
        return logResult;
    }

    public void calculateKsStats(int minsize, int maxsize) {
        KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();

        // Collect background fold changes (from relevantGenes) for valid GO entries
        HashSet<String> validBackgroundGenes = new HashSet<>();
        for (GOEntry goEntry : dag.values()) {
            if (goEntry.getSize() >= minsize && goEntry.getSize() <= maxsize) {
                validBackgroundGenes.addAll(goEntry.getEnrichedGenes());
            }
        }

        // Map valid background genes to their fold-change values
        double[] backgroundFoldChanges = relevantGenes.stream()
                //.filter(geneMap::containsKey) // Ensure the gene exists in geneMap
                .mapToDouble(geneId -> geneMap.get(geneId).getFoldChange())
                .toArray();

        for (GOEntry goEntry : dag.values()) {
            // Skip GO entries that don't meet size constraints
            if (goEntry.getSize() < minsize || goEntry.getSize() > maxsize) {
                //goEntry.setKsStat(0.0);
                //goEntry.setKsPval(1.0);
                continue;
            }

            // Get fold changes for genes in enrichedGenes (specific to this GO term)
            double[] goFoldChanges = goEntry.getEnrichedGenes().stream()
                    //.filter(geneMap::containsKey) // Ensure the gene exists in geneMap
                    .mapToDouble(geneId -> geneMap.get(geneId).getFoldChange())
                    .toArray();
            // Skip if there are no genes in the GO term
            if (goFoldChanges.length == 0) {
                goEntry.setKsStat(0.0);
                goEntry.setKsPval(1.0);
                continue;
            }

            // Perform the KS test
            double ksStat = ksTest.kolmogorovSmirnovStatistic(goFoldChanges, backgroundFoldChanges);
            double ksPval = ksTest.kolmogorovSmirnovTest(goFoldChanges, backgroundFoldChanges);

            // Store the results in the GO entry
            goEntry.setKsStat(ksStat);
            goEntry.setKsPval(ksPval);
        }
    }

    public void applyKsBenjaminiHochbergCorrection(int minsize, int maxsize) {
        // Filter GO entries based on size constraints and create a mutable list
        List<GOEntry> filteredEntries = new ArrayList<>(
                dag.values().stream()
                        .filter(goEntry -> goEntry.getSize() >= minsize && goEntry.getSize() <= maxsize)
                        .toList()
        );

        // Sort by KS p-values in ascending order
        filteredEntries.sort(Comparator.comparingDouble(GOEntry::getKsPval));

        int m = filteredEntries.size(); // Total number of valid tests
        double[] correctedPvals = new double[m];

        // Apply Benjamini-Hochberg correction
        for (int i = 0; i < m; i++) {
            double rawPval = filteredEntries.get(i).getKsPval();
            correctedPvals[i] = rawPval * m / (i + 1); // BH formula
        }

        // Ensure corrected p-values are monotonically increasing
        for (int i = m - 2; i >= 0; i--) {
            correctedPvals[i] = Math.min(correctedPvals[i], correctedPvals[i + 1]);
        }

        // Assign corrected p-values back to the filtered GO entries
        for (int i = 0; i < m; i++) {
            filteredEntries.get(i).setKsFdr(Math.min(correctedPvals[i], 1.0)); // Ensure FDR ≤ 1.0
        }
    }

    public String calculateShortestPathToTrue(GOEntry goEntry, List<GOEntry> trueEntries) {
        // Return empty if the entry is already true
        if (goEntry.isTrue()) {
            return "";
        }


        // Return empty if no true entries are present
        if (trueEntries.isEmpty()) {
            return "";
        }

        // Perform BFS to find the shortest path to a true entry
        Queue<List<GOEntry>> queue = new LinkedList<>();
        Set<GOEntry> visited = new HashSet<>();

        // Start BFS from the analyzed GO entry
        queue.add(Collections.singletonList(goEntry));
        visited.add(goEntry);

        List<GOEntry> shortestPath = null;
        GOEntry lca = null;

        while (!queue.isEmpty()) {
            List<GOEntry> path = queue.poll();
            GOEntry current = path.get(path.size() - 1);

            // If the current node is a true entry, we found the shortest path
            if (current.isTrue()) {
                shortestPath = path;
                lca = findLCA(path, current); // Find the LCA for the path
                break;
            }

            //Traverse neighbors (parents and children)
            List<GOEntry> neighbors = new ArrayList<>();
            neighbors.addAll(current.getParents());
            neighbors.addAll(current.getChildren());

            for (GOEntry neighbor : neighbors) {
                if (!visited.contains(neighbor)) {
                    visited.add(neighbor);

                    // Create a new path with the neighbor added
                    List<GOEntry> newPath = new ArrayList<>(path);
                    newPath.add(neighbor);

                    // Enqueue the new path
                    queue.add(newPath);
                }
            }
        }

        // If no path was found, return empty
        if (shortestPath == null) {
            return "";
        }

        // Format the path as a string with '|' separators and mark the LCA with '*'
        StringBuilder result = new StringBuilder();
        for (GOEntry entry : shortestPath) {
            result.append(entry.getName());
            if (entry.equals(lca)) {
                result.append(" * ");
            }
            result.append("|");
        }

        // Remove trailing '|'
        result.setLength(result.length() - 1);

        return result.toString();
    }

    /**
     * Find the Least Common Ancestor (LCA) for a given path.
     */
    private GOEntry findLCA(List<GOEntry> path, GOEntry trueEntry) {
        // Track ancestors of the true entry
        Set<GOEntry> trueAncestors = new HashSet<>();
        Queue<GOEntry> queue = new LinkedList<>();
        queue.add(trueEntry);

        while (!queue.isEmpty()) {
            GOEntry current = queue.poll();
            trueAncestors.add(current);
            queue.addAll(current.getParents());
        }

        // Traverse the path from the analyzed entry and find the first common ancestor
        for (GOEntry entry : path) {
            if (trueAncestors.contains(entry)) {
                return entry;
            }
        }

        // If no common ancestor is found (which shouldn't happen), return null
        return null;
    }

    public String newCalculateShortestPathToTrue(GOEntry goEntry, List<GOEntry> trueEntries) {
        // Return empty if the entry is already true
        if (goEntry.isTrue()) {
            return "";
        }

        // Return empty if no true entries are present
        if (trueEntries.isEmpty()) {
            return "";
        }

        // Perform BFS to find the shortest path to a true entry
        Queue<PathState> queue = new LinkedList<>();
        Set<GOEntry> visited = new HashSet<>();

        // Start BFS from the analyzed GO entry
        queue.add(new PathState(Collections.singletonList(goEntry), Direction.UP));
        visited.add(goEntry);

        List<GOEntry> shortestPath = null;
        GOEntry lca = null;

        while (!queue.isEmpty()) {
            PathState state = queue.poll();
            List<GOEntry> path = state.path;
            Direction direction = state.direction;
            GOEntry current = path.get(path.size() - 1);

            // If the current node is a true entry, we found the shortest path
            if (current.isTrue()) {
                shortestPath = path;
                lca = findLCA(path, current); // Find the LCA for the path
                break;
            }

            // Traverse neighbors (parents or children based on the current direction)
            List<GOEntry> neighbors = new ArrayList<>();
            if (direction == Direction.UP) {
                // Add parents (going up)
                for (GOEntry parent : current.getParents()) {
                    neighbors.add(parent);
                }
            }

            // Always add children (going down)
            for (GOEntry child : current.getChildren()) {
                neighbors.add(child);
            }

            for (GOEntry neighbor : neighbors) {
                if (!visited.contains(neighbor)) {
                    visited.add(neighbor);

                    // Determine the new direction
                    Direction newDirection = direction;
                    if (current.getParents().contains(neighbor)) {
                        newDirection = Direction.UP;
                    } else if (current.getChildren().contains(neighbor)) {
                        newDirection = Direction.DOWN;
                    }

                    // Enforce the directional constraint: Once going down, never go up again
                    if (direction == Direction.DOWN && newDirection == Direction.UP) {
                        continue; // Skip neighbors that violate the constraint
                    }

                    // Create a new path with the neighbor added
                    List<GOEntry> newPath = new ArrayList<>(path);
                    newPath.add(neighbor);

                    // Enqueue the new path
                    queue.add(new PathState(newPath, newDirection));
                }
            }
        }

        // If no path was found, return empty
        if (shortestPath == null) {
            return "";
        }

        // Format the path as a string with '|' separators and mark the LCA with '*'
        StringBuilder result = new StringBuilder();
        for (GOEntry entry : shortestPath) {
            result.append(entry.getName());
            if (entry.equals(lca)) {
                result.append(" * ");
            }
            result.append("|");
        }

        // Remove trailing '|'
        result.setLength(result.length() - 1);

        return result.toString();
    }

    enum Direction {
        UP, DOWN
    }

    class PathState {
        List<GOEntry> path;
        Direction direction;

        PathState(List<GOEntry> path, Direction direction) {
            this.path = path;
            this.direction = direction;
        }
    }

    private Map<String, List<GOEntry>> buildGeneToGoMap() {
        Map<String, List<GOEntry>> geneToGoMap = new HashMap<>();

        for (GOEntry goEntry : dag.values()) {
            // Only consider entries that satisfy the size constraints
            if (goEntry.getSize() >= minsize && goEntry.getSize() <= maxsize) {
                for (String gene : goEntry.getMappedGenes()) {
                    geneToGoMap.computeIfAbsent(gene, k -> new ArrayList<>()).add(goEntry);
                }
            }
        }
        return geneToGoMap;
    }

    public void writeOverlapOut(String overlapOutFile) {
        //Map<String, List<GOEntry>> geneToGoMap = buildGeneToGoMap();
//        System.out.println("start");
//        HashSet<String> overlapPairs = buildOverlapPairs();
//        System.out.println(overlapPairs.size());

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(overlapOutFile))) {
            // Write the header
            writer.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent\n");

            // Store visited pairs to avoid duplicates
            Set<String> processedPairs = new HashSet<>();
            int count = 0;

            // Iterate over each gene and its associated GO entries
            //for (String gene : overlapPairs) {
            for (IdPair pair : relPairs) {
                // Compare each pair of GO entries associated with the current gene
                GOEntry term1 = dag.get(pair.getFirst());
                GOEntry term2 = dag.get(pair.getSecond());
//                GOEntry term1 = goEntries.get(i);

                if (term1.getMappedGenes().size() < minsize || term1.getMappedGenes().size() > maxsize || term2.getMappedGenes().size() < minsize || term2.getMappedGenes().size() > maxsize) {
                    continue;
                }

                count++;
                if (count % 1000 == 0) {
                    //System.out.println(count);
                }
//                GOEntry term2 = goEntries.get(j);

                // Skip if this pair has already been processed
//                String pairKey = term1.getId() + "-" + term2.getId();
//                if (processedPairs.contains(pairKey)) {
//                    continue;
//                }
//                processedPairs.add(pairKey);

                // Calculate overlap metrics
                int numOverlapping = calculateNumOverlapping(term1, term2);
                if (numOverlapping == 0) {
                    continue; // Skip pairs with no overlapping genes
                }
                double maxOvPercent = calculateMaxOvPercent(term1, term2, numOverlapping);
                boolean isRelative = isRelative(term1, term2);
                //boolean isRelative = false;
                //int pathLength = calculatePathLength(term1, term2);
                int pathLength = 0;
                // Write the result
                writer.write(String.format("%s\t%s\t%s\t%d\t%d\t%.2f\n",
                        term1.getId(),
                        term2.getId(),
                        isRelative,
                        pathLength,
                        numOverlapping,
                        maxOvPercent
                ));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private int calculateNumOverlapping(GOEntry term1, GOEntry term2) {
//        Set<String> set1 = new HashSet<>(term1.getMappedGenes());
//        Set<String> set2 = new HashSet<>(term2.getMappedGenes());
//
//        set1.retainAll(set2); // Keep only the intersection
        return (int) term1.getMappedGenes().stream().filter(term2.getMappedGenes()::contains).count();
        //return set1.size();
    }

    private boolean isRelative(GOEntry term1, GOEntry term2) {
        return isAncestor(term1, term2) || isAncestor(term2, term1);
    }

    private boolean isAncestor(GOEntry ancestor, GOEntry descendant) {
        Set<GOEntry> visited = new HashSet<>();
        Queue<GOEntry> queue = new LinkedList<>();
        queue.add(descendant);

        while (!queue.isEmpty()) {
            GOEntry current = queue.poll();
            if (current.equals(ancestor)) {
                return true;
            }
            for (GOEntry parent : current.getParents()) {
                if (!visited.contains(parent)) {
                    visited.add(parent);
                    queue.add(parent);
                }
            }
        }

        return false;
    }


    private double calculateMaxOvPercent(GOEntry term1, GOEntry term2, int numOverlapping) {
        int size1 = term1.getMappedGenes().size();
        int size2 = term2.getMappedGenes().size();

        double percent1 = (numOverlapping * 100.0) / size1;
        double percent2 = (numOverlapping * 100.0) / size2;

        return Math.max(percent1, percent2);
    }

    private int calculatePathLength(GOEntry term1, GOEntry term2) {
        Queue<GOEntry> queue = new LinkedList<>();
        Map<GOEntry, Integer> distances = new HashMap<>();
        queue.add(term1);
        distances.put(term1, 0);

        while (!queue.isEmpty()) {
            GOEntry current = queue.poll();
            int currentDistance = distances.get(current);

            if (current.equals(term2)) {
                return currentDistance;
            }

            // Traverse neighbors
            for (GOEntry neighbor : current.getParents()) {
                if (!distances.containsKey(neighbor)) {
                    distances.put(neighbor, currentDistance + 1);
                    queue.add(neighbor);
                }
            }
            for (GOEntry neighbor : current.getChildren()) {
                if (!distances.containsKey(neighbor)) {
                    distances.put(neighbor, currentDistance + 1);
                    queue.add(neighbor);
                }
            }
        }

        // Return a large number if no path exists
        return Integer.MAX_VALUE;
    }
}