import org.apache.commons.cli.*;

public class Main {
    public static void main(String[] args) throws ParseException {
        System.out.println("Hello, World!");

        Options options = new Options();
        options.addOption("obo", true, "obo file");
        options.addOption("root", true, "namespace");
        options.addOption("mapping", true, "mapping file");
        options.addOption("mappingtype", true, "mapping type");
        options.addOption("enrich", true, "output path");
        options.addOption("o", true, "output path");
        options.addOption("minsize", true, "GTF file path");
        options.addOption("maxsize", true, "GTF file path");
        options.addOption(Option.builder("overlapout")
                .hasArg()
                .optionalArg(true)
                .desc("strandness of experiment")
                .build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);

        if (cmd.getOptions().length == 0) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Main", options);
            return;
        }

        boolean overlap = false;
        String overlapOut = null;
        if (cmd.hasOption("overlapout")) {
            overlap = true;
            overlapOut= cmd.getOptionValue("overlapout");
        }

        String obo = cmd.getOptionValue("obo");
        String root = cmd.getOptionValue("root");
        String mapping = cmd.getOptionValue("mapping");
        boolean mapGo = cmd.getOptionValue("mappingtype").equals("go");
        String enrich = cmd.getOptionValue("enrich");
        String o = cmd.getOptionValue("o");
        int minsize = Integer.parseInt(cmd.getOptionValue("minsize"));
        int maxsize = Integer.parseInt(cmd.getOptionValue("maxsize"));
    }
}