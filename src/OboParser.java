import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class OboParser {
    private String obo;
    private String rootName;
    private HashMap<String, GOEntry> dag;
    private GOEntry rootEntry; // Store the root GOEntry

    public OboParser(String obo, String rootName) {
        this.obo = obo;
        this.rootName = rootName;
        this.dag = new HashMap<>();
        this.rootEntry = null;
    }

    public HashMap<String, GOEntry> getDag() {
        try (BufferedReader reader = new BufferedReader(new FileReader(obo))) {
            String line;
            GOEntry currentEntry = null;
            boolean isObsolete = false;
            boolean readMode = false;
            boolean toBeAdded = false;

            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if (line.startsWith("[Term]")) {
                    toBeAdded = true;
                    readMode = true;
                    currentEntry = null;
                    isObsolete = false;
                } else if (readMode && line.isEmpty()) {
                    readMode = false;
                    if (toBeAdded && currentEntry != null) {
                        dag.put(currentEntry.getId(), currentEntry);
                    }
                } else if (readMode && line.startsWith("id: ") && currentEntry == null) {
                    String id = line.substring(4);
                    if (dag.containsKey(id)) {
                        currentEntry = dag.get(id);
                        toBeAdded = false;
                    } else {
                        currentEntry = new GOEntry(id, null);
                    }
                } else if (readMode && line.startsWith("name: ") && currentEntry != null) {
                    String name = line.substring(6);
                    currentEntry.setName(name);

                    // Check if this entry is the root
                    if (name.equals(rootName)) {
                        rootEntry = currentEntry;
                    }
                } else if (readMode && line.startsWith("namespace: ") && currentEntry != null) {
                    if (!line.substring(11).equals(rootName)) {
                        currentEntry = null;
                        toBeAdded = false;
                    }
                } else if (readMode && line.startsWith("is_obsolete: true")) {
                    toBeAdded = false;
                    isObsolete = true;
                } else if (readMode && line.startsWith("is_a: ") && currentEntry != null && !isObsolete) {
                    String parentId = line.substring(6, 16);
                    GOEntry parent = dag.getOrDefault(parentId, new GOEntry(parentId, null));
                    dag.putIfAbsent(parentId, parent);
                    currentEntry.addParent(parent);
                    parent.addChild(currentEntry);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return dag;
    }

    // Method to get the root GOEntry
    public GOEntry getRootEntry() {
        return rootEntry;
    }
}
