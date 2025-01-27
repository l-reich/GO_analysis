import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class GOEntry {
    private String id;
    private String name;
    private List<GOEntry> parents;
    private List<GOEntry> children;
    private boolean isTrue;
    private HashSet<String> mappedGenes;
    private HashSet<String> enrichedGenes;
    private int noverlap;
    private double hgPval;
    private double hgFdr;
    private double fejPval;
    private double fejFdr;
    private double ksStat;
    private double ksPval;
    private double ksFdr;
    private String shortestPathToATrue;
    private HashSet<String> overlaps;

    public GOEntry(String id, String name) {
        this.id = id;
        this.name = name;
        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
        this.isTrue = false;
        this.mappedGenes = new HashSet<>();
        this.noverlap = 0;
    }

    public void setNoverlap(int noverlap) {
        this.noverlap = noverlap;
    }

    public int getNoverlap() {
        return noverlap;
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

    public boolean isTrue() {
        return isTrue;
    }

    public void setTrue(boolean isTrue) {
        this.isTrue = isTrue;
    }

    public HashSet<String> getMappedGenes() {
        return mappedGenes;
    }

    //method to add mapped genes
    public void addMappedGene(String gene) {
        if (mappedGenes == null) {
            mappedGenes = new HashSet<>();
        }
        mappedGenes.add(gene);
    }

    public HashSet<String> getEnrichedGenes() {
        return enrichedGenes;
    }

    public void setEnrichedGenes(HashSet<String> enrichedGenes) {
        this.enrichedGenes = enrichedGenes;
    }

    public void addMappedGenes(Set<String> genes) {
        if (genes != null) {
            this.mappedGenes.addAll(genes);
        }
    }

    public double getHgPval() {
        return hgPval;
    }

    public void setHgPval(double hgPval) {
        this.hgPval = hgPval;
    }

    public double getHgFdr() {
        return hgFdr;
    }

    public void setHgFdr(double hgFdr) {
        this.hgFdr = hgFdr;
    }

    public double getFejPval() {
        return fejPval;
    }

    public void setFejPval(double fejPval) {
        this.fejPval = fejPval;
    }

    public double getFejFdr() {
        return fejFdr;
    }

    public void setFejFdr(double fejFdr) {
        this.fejFdr = fejFdr;
    }

    public int getSize() {
        return enrichedGenes.size();
    }


    public double getKsStat() {
        return ksStat;
    }

    public void setKsStat(double ksStat) {
        this.ksStat = ksStat;
    }

    public double getKsPval() {
        return ksPval;
    }

    public void setKsPval(double ksPval) {
        this.ksPval = ksPval;
    }

    public double getKsFdr() {
        return ksFdr;
    }

    public void setKsFdr(double ksFdr) {
        this.ksFdr = ksFdr;
    }

    public void setShortestPathToATrue(String shortestPath) {
        this.shortestPathToATrue = shortestPath;
    }

    public String getShortestPathToATrue() {
        return shortestPathToATrue;
    }

    public void setOverlaps(HashSet<String> overlaps) {
        this.overlaps = overlaps;
    }

    public HashSet<String> getOverlaps() {
        return overlaps;
    }

    public void addOverlap(String overlap) {
        if (overlaps == null) {
            overlaps = new HashSet<>();
        }
        overlaps.add(overlap);
    }

    public void addOverlaps(HashSet<String> overlaps) {
        if (overlaps != null) {
            this.overlaps.addAll(overlaps);
        }
    }
}