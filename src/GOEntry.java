import java.util.ArrayList;
import java.util.List;

public class GOEntry {
    private String id;
    private String name;
    private List<GOEntry> parents;
    private List<GOEntry> children;

    public GOEntry(String id, String name) {
        this.id = id;
        this.name = name;
        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<GOEntry> getParents() {
        return parents;
    }

    public List<GOEntry> getChildren() {
        return children;
    }

    public void addParent(GOEntry parent) {
        this.parents.add(parent);
    }

    public void addChild(GOEntry child) {
        this.children.add(child);
    }
}