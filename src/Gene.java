public class Gene {
    private String id;
    private double foldChange;
    private boolean signif;

    public Gene(String id, double foldChange, boolean signif) {
        this.id = id;
        this.foldChange = foldChange;
        this.signif = signif;
    }

    public String getId() {
        return id;
    }

    public double getFoldChange() {
        return foldChange;
    }

    public boolean isSignif() {
        return signif;
    }
}