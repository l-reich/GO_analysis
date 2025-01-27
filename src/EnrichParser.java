import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class EnrichParser {
    private String filePath;
    private HashMap<String, GOEntry> goEntries;
    private HashMap<String, Gene> genes;

    public EnrichParser(String filePath, HashMap<String, GOEntry> goEntries) {
        this.filePath = filePath;
        this.goEntries = goEntries;
        this.genes = new HashMap<>();
    }

    public void parseFile() {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean readingGOEntries = true;

            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if (line.startsWith("#")) {
                    String goId = line.substring(1);
                    if (goEntries.containsKey(goId)) {
                        goEntries.get(goId).setTrue(true);
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
                        genes.put(id, gene);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public HashMap<String, Gene> getGenes() {
        return genes;
    }
}