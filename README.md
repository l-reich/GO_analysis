# GO_analysis

GO_analysis is a command-line Java tool for performing Gene Ontology (GO) enrichment analysis and processing gene-to-GO term mappings. It supports statistical evaluation (hypergeometric, Fisher's exact, and Kolmogorov-Smirnov tests) and multiple testing correction, helping users identify enriched GO terms among gene sets.

## Features

- Parses GO term relationships from `.obo` files.
- Processes gene-to-GO mappings in various formats.
- Supports enrichment analysis using hypergeometric, Fisher's exact, and KS tests.
- Applies Benjamini-Hochberg correction for multiple testing.
- Handles custom gene lists and mapping types.
- Outputs results and enrichment statistics to files.

## Usage

1. **Build the project** (ensure you have Java and any required dependencies):
   ```bash
   javac -d bin src/*.java
   ```

2. **Run the program** with required options:
   ```bash
   java -cp bin Main \
     -obo path/to/go.obo \
     -root NAMESPACE \
     -mapping path/to/mapping_file \
     -mappingtype [go|ens] \
     -enrich path/to/enrichment_file \
     -o output_path \
     -minsize MIN_GO_SIZE \
     -maxsize MAX_GO_SIZE
   ```
   **Optional:**  
   Add `-overlapout path/to/overlap_output` to include overlap statistics.

### Arguments

- `-obo` : Path to GO `.obo` ontology file.
- `-root` : GO namespace/root term (e.g., BP/CC/MF).
- `-mapping` : File with gene-to-GO mappings.
- `-mappingtype` : Mapping type (`go` or `ens`).
- `-enrich` : File with input gene list (for enrichment).
- `-o` : Output file path for results.
- `-minsize` : Minimum GO term size to consider.
- `-maxsize` : Maximum GO term size to consider.
- `-overlapout` : (Optional) Output file for overlap results.

To display help:
```bash
java -cp bin Main
```

## Input File Formats

- **OBO file**: Standard GO ontology format.
- **Mapping file**: Tab-separated gene-to-GO relationships.
- **Enrichment file**: List of genes (optionally with fold change and significance).

## Output

The tool generates enrichment results including GO terms, associated genes, p-values, FDR, and optionally overlap statistics.
