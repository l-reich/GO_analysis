import java.util.HashMap;

public class DAGHandler {
    private OboParser oboParser;
    private MappingParser mappingParser;
    private EnrichParser enrichParser;
    private String root;
    private String o;
    private int minsize;
    private int maxsize;
    private boolean overlap;
    private String overlapOut;
    private HashMap<String, GOEntry> dag;
    private HashMap<String, String> mapping;
    private HashMap<String, String> enrich;

    public DAGHandler(String obo, String root, String mapping, boolean mapGo, String enrich, String o, int minsize, int maxsize, boolean overlap, String overlapOut) {
        this.oboParser = new OboParser(obo, root);
//        this.mappingParser = new MappingParser(mapping, mapGo);
//        this.enrichParser = new EnrichParser(enrich);
        this.root = root;
        this.o = o;
        this.minsize = minsize;
        this.maxsize = maxsize;
        this.overlap = overlap;
        this.overlapOut = overlapOut;
        dag = oboParser.getDag();
//        mapping = mappingParser.getMapping();
//        enrich = enrichParser.getEnrich();
    }
}
