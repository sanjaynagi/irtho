# Irtho

<div align="center">

[<img src="https://raw.githubusercontent.com/sanjaynagi/irtho/main/irtho.png" width="400"/>](https://raw.githubusercontent.com/sanjaynagi/irtho/main/irtho.png)   

</div>

## Overview

**IRTHO** is a python package to find orthologous loci between reference genomes, based on OrthoFinder.

## Features

- **Ortholog Identification**: Identify orthologous genes between species using OrthoFinder.
- **Determine synteny**: Determine synteny between orthologous genes in each species.
- **Sequence Alignment**: Align protein sequences to find orthologous positions of amino acid residues.

## Usage

```python
import irtho
import pandas as pd

# Initialize Orthologs object with OrthoFinder results
ortho = irtho.Orthologs(results_dir="path/to/orthofinder/results")

# Load input targets DataFrame
targets_df = pd.read_csv("targets.tsv", sep="\t")
```

## Key Features

### 1. Mapping Orthologs Between Species
```python
# Map genes from reference to target species
mapped_df = ortho.map_input_genes_to_orthologs(
    input_df=targets_df,
    reference_species="AgambiaePEST",
    target_species="AaegyptiLVP_AGWG"
)
```

### 2. Finding Orthologous Positions
```python
# Find orthologous genomic positions
final_targets = ortho.find_orthologous_targets(
    targets_df=mapped_df,
    reference_dir="path/to/reference/files",
    ref_genome="AgambiaePEST",
    target_genome="AaegyptiLVP_AGWG"
)
```

### 3. Random Exon Selection
For exploring new potential targets, Irtho can select random positions within exons:
```python
# Process targets with 'randomN' in codon column
random_targets = ortho.process_random_targets(
    targets_df=mapped_df,
    reference_dir="path/to/reference/files",
    target_species="AaegyptiLVP_AGWG"
)
```

## Input Format
The input targets DataFrame should contain the following columns:
- `gene`: Reference gene ID
- `transcript`: Reference transcript ID (optional)
- `codon`: Codon number or 'randomN' for random exon selection
- Additional columns as needed

## Example Workflow

```python
# Complete workflow example
import irtho

# Initialize
ortho = irtho.Orthologs("results/orthofinder/Results_Feb03_1")

# Load targets
targets_df = pd.read_csv("irtho-targets.tsv", sep="\t")

# Map orthologs
mapped_df = ortho.map_input_genes_to_orthologs(
    targets_df,
    reference_species="AgambiaePEST",
    target_species="AaegyptiLVP_AGWG"
)

# Process random targets (if any)
processed_df = ortho.process_random_targets(
    mapped_df,
    reference_dir="resources/reference",
    target_species="AaegyptiLVP_AGWG"
)

# Find orthologous positions
final_df = ortho.find_orthologous_targets(
    processed_df,
    reference_dir="resources/reference",
    ref_genome="AgambiaePEST",
    target_genome="AaegyptiLVP_AGWG"
)

# Export results
final_df.to_csv("orthologous_targets.tsv", sep="\t", index=False)
```