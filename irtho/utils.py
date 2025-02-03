import os
from Bio import SeqIO
import pandas as pd

def load_fasta(fasta_path, debug=False):
    """
    Load sequences from a FASTA file into a dictionary.
    
    Args:
        fasta_path: Path to the FASTA file.
        debug: Whether to print debug information.
    
    Returns:
        A dictionary of sequences keyed by protein ID.
    """
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    
    if debug:
        print(f"Loaded {len(sequences)} sequences from {fasta_path}")
        print(f"Sample IDs: {list(sequences.keys())[:3]}")
        print(f"Length range: {min(len(seq) for seq in sequences.values())} - {max(len(seq) for seq in sequences.values())}")
    
    return sequences

def create_gene_mapping(file_path):
    """
    Parse a GFF file and create a mapping from genes to transcripts to proteins.
    
    Args:
        file_path (str): Path to the GFF file.
    
    Returns:
        pd.DataFrame: A DataFrame with columns ['gene_id', 'transcript_id', 'protein_id'].
    """
    mappings = []
    with open(file_path, 'r') as file:
        gene_id = None
        transcript_id = None

        for line in file:
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip():
                continue

            # Split the line into columns
            columns = line.strip().split("\t")
            if len(columns) < 9:
                continue  # Skip malformed lines

            feature_type = columns[2]  # The third column indicates the feature type
            attributes = columns[8]  # The ninth column contains the attributes

            # Parse the attributes into a dictionary
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key.strip()] = value.strip()

            # Extract gene_id and transcript_id from mRNA or transcript rows
            if feature_type in {"mRNA", "transcript"}:
                transcript_id = attr_dict.get("ID")
                gene_id = attr_dict.get("Parent")

            # Extract protein_id from CDS rows
            elif feature_type == "CDS":
                protein_id = attr_dict.get("protein_source_id") or attr_dict.get("protein_id")
                if gene_id and transcript_id and protein_id:
                    mappings.append((gene_id, transcript_id, protein_id))

    # Create a DataFrame from the mappings
    df = pd.DataFrame(mappings, columns=["gene_id", "transcript_id", "protein_id"])
    df.loc[:, 'gene_id'] = df.gene_id.str.replace("gene-", "")
    df.loc[:, 'transcript_id'] = df.transcript_id.str.replace("rna-", "")
    return df.drop_duplicates().reset_index(drop=True)
    
def get_longest_transcripts(sequences, mapping_df, debug=False):
    """
    Find the longest transcript for each gene using the mapping DataFrame.

    Args:
        sequences: Dictionary of protein sequences keyed by protein ID.
        mapping_df: DataFrame with columns ['gene_id', 'transcript_id', 'protein_id'].
        debug: Whether to print debug information.

    Returns:
        A dictionary of longest transcripts keyed by gene ID.
    """
    # Filter the mapping DataFrame to include only rows where the protein_id exists in the sequences
    mapping_df = mapping_df[mapping_df["protein_id"].isin(sequences.keys())]

    if debug:
        print(f"Filtered mapping file to {len(mapping_df)} rows with valid protein IDs.")
        print(f"Sample mapping rows:\n{mapping_df.head()}")

    # Group transcripts by gene ID
    gene_transcripts = mapping_df.groupby("gene_id").apply(
        lambda group: group[["transcript_id", "protein_id"]].to_dict(orient="records")
    ).to_dict()

    if debug:
        print(f"Found {len(gene_transcripts)} unique genes in the mapping file.")
        print(f"Sample gene-transcript mappings: {list(gene_transcripts.items())[:3]}")

    # Keep the longest transcript for each gene
    longest_transcripts = {}
    for gene_id, transcript_records in gene_transcripts.items():
        # Map protein IDs to sequences and find the longest
        valid_transcripts = [
            (record["transcript_id"], record["protein_id"], sequences[record["protein_id"]])
            for record in transcript_records
            if record["protein_id"] in sequences
        ]

        if not valid_transcripts:
            if debug:
                print(f"Warning: No valid transcripts found for gene {gene_id}")
            continue

        # Find the longest transcript
        longest = max(valid_transcripts, key=lambda x: len(x[2]))
        longest_transcripts[gene_id] = (longest[0], longest[2])  # (transcript_id, sequence)

        if debug and len(valid_transcripts) > 1:
            lengths = [(t[0], len(t[2])) for t in valid_transcripts]
            print(f"\nGene {gene_id} transcripts:")
            for tid, length in lengths:
                print(f"  {tid}: {length} aa{'*' if tid == longest[0] else ''}")

    return longest_transcripts

def write_fasta(sequences, output_path, debug=False):
    """
    Write sequences to a FASTA file with line wrapping at 60 characters.

    Args:
        sequences: Dictionary of sequences keyed by gene ID.
        output_path: Path to the output FASTA file.
        debug: Whether to print debug information.
    """
    with open(output_path, 'w') as f:
        for gene_id, (transcript_id, sequence) in sequences.items():
            f.write(f">{transcript_id}\n")
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")
    
    if debug:
        print(f"Wrote {len(sequences)} sequences to {output_path}")
        print(f"Total residues written: {sum(len(seq[1]) for seq in sequences.values())}")

def write_longest_isoforms(genome, debug=False, path_to_references="../resources/reference"):
    """
    Combines the workflow to extract and write the longest isoforms for a given genome.

    Args:
        genome (str): The genome identifier (e.g., "AgambiaePEST").
        debug (bool): Whether to print debug information.
        path_to_references (str): Path to the reference directory containing GFF and FASTA files.

    Returns:
        None
    """
    # Define file paths
    fasta_path = os.path.join(path_to_references, f"{genome}_AnnotatedProteins.fasta")
    gff_path = os.path.join(path_to_references, f"{genome}.gff")
    output_fasta_path = f"../results/proteome/{genome}.fasta"

    if debug:
        print(f"Processing genome: {genome}")
        print(f"FASTA file: {fasta_path}")
        print(f"GFF file: {gff_path}")
        print(f"Output FASTA file: {output_fasta_path}")

    # Step 1: Load protein sequences from the FASTA file
    sequences = load_fasta(fasta_path, debug=debug)

    # Step 2: Create a gene-to-transcript mapping from the GFF file
    df_mapping = create_gene_mapping(gff_path)

    # Step 3: Identify the longest isoforms for each gene
    longest_transcripts = get_longest_transcripts(sequences, df_mapping, debug=debug)

    # Step 4: Write the longest isoforms to a new FASTA file
    write_fasta(longest_transcripts, output_fasta_path, debug=debug)

    if debug:
        print(f"Longest isoforms written to: {output_fasta_path}")
