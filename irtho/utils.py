import os
from Bio import SeqIO
import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_tqdm():
    """
    Return the appropriate tqdm version depending on the environment.
    
    In a Jupyter Notebook (IPython with ZMQInteractiveShell) this returns 
    tqdm.notebook.tqdm; otherwise, it returns the standard tqdm.tqdm.
    
    Returns:
        The tqdm function.
    """
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            # likely running in a Jupyter notebook or qtconsole
            from tqdm.notebook import tqdm
        else:
            # likely running in a terminal or other interface
            from tqdm import tqdm
    except (ImportError, NameError):
        # Fallback to the default CLI version if IPython is not available.
        from tqdm import tqdm

    return tqdm

def split_one_to_many_orthologs(targets_df, target_genome):
    # split one to many ortholog targets 
    """
    Split a column with one-to-many ortholog targets into individual rows for each possible ortholog.
    
    Args:
        targets_df (pd.DataFrame): DataFrame containing a column of target genes.
        target_genome (str): The name of the column to split.
    
    Returns:
        pd.DataFrame: A DataFrame with the same columns as the original,
            but with the one to many orthologs split into individual rows.
    """
    split_targets_df = targets_df.copy()
    split_targets_df[target_genome] = split_targets_df[target_genome].str.split(',')
    split_targets_df = split_targets_df.explode(target_genome).reset_index(drop=True)
    split_targets_df[target_genome] = split_targets_df[target_genome].str.replace(" ", "")
    
    return split_targets_df


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

import pandas as pd

def parse_attributes(attr_string):
    """Parse GFF attribute string into dictionary."""
    attributes = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    return attributes

def create_gene_mapping(gff_file):
    """
    Create a pandas DataFrame mapping genes to their transcripts and proteins.
    Handles both VectorBase and RefSeq GFF formats.
    
    Args:
        gff_file (str): Path to the GFF file
        
    Returns:
        pandas.DataFrame: DataFrame with columns for gene_id, transcript_id, and protein_id
    """
    gene_transcript_protein_map = []
    seen = set()  # To track unique combinations
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            
            source = fields[1]
            attributes = parse_attributes(fields[8])
            
            # Handle VectorBase format
            if source == 'VEuPathDB':
                gene_id = attributes.get('gene_id')
                transcript_id = attributes.get('Parent')
                protein_id = attributes.get('protein_source_id')
                
            # Handle RefSeq format - only process if source is RefSeq
            elif source == 'RefSeq':
                gene_id = attributes.get('gene')  # LOC ID
                transcript_id = attributes.get('Parent')
                if transcript_id and transcript_id.startswith('rna-'):
                    transcript_id = transcript_id[4:]  # Remove 'rna-' prefix
                protein_id = attributes.get('protein_id')
                if protein_id and protein_id.startswith('cds-'):
                    protein_id = protein_id[4:]  # Remove 'cds-' prefix
                    
            else:
                continue  # Skip all other sources
                
            # Skip if any required field is missing
            if not all([gene_id, transcript_id, protein_id]):
                continue
                
            # Create tuple of the mapping
            mapping = (gene_id, transcript_id, protein_id)
            
            # Only add if we haven't seen this combination before
            if mapping not in seen:
                gene_transcript_protein_map.append({
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                    'protein_id': protein_id
                })
                seen.add(mapping)
    
    # Convert to DataFrame
    df = pd.DataFrame(gene_transcript_protein_map)
    
    # Sort by gene_id to maintain consistent order
    df = df.sort_values('gene_id').reset_index(drop=True)
    
    return df
    
def get_longest_transcripts(sequences, mapping_df, debug=False):
    """
    Find the longest transcript for each gene using the mapping DataFrame and return
    sequences keyed by gene ID.

    Args:
        sequences: Dictionary of protein sequences keyed by protein ID.
        mapping_df: DataFrame with columns ['gene_id', 'transcript_id', 'protein_id'].
        debug: Whether to print debug information.

    Returns:
        A dictionary of longest protein sequences keyed by gene ID.
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
            (record["protein_id"], sequences[record["protein_id"]])
            for record in transcript_records
            if record["protein_id"] in sequences
        ]

        if not valid_transcripts:
            if debug:
                print(f"Warning: No valid transcripts found for gene {gene_id}")
            continue

        # Find the longest transcript
        longest = max(valid_transcripts, key=lambda x: len(x[1]))
        longest_transcripts[gene_id] = longest[1]  # Store the sequence keyed by gene_id

        if debug and len(valid_transcripts) > 1:
            lengths = [(t[0], len(t[1])) for t in valid_transcripts]
            print(f"\nGene {gene_id} transcripts:")
            for pid, length in lengths:
                print(f"  {pid}: {length} aa{'*' if pid == longest[0] else ''}")

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
        for gene_id, sequence in sequences.items():
            f.write(f">{gene_id}\n")  # Use gene_id as the FASTA header
            for i in range(0, len(sequence), 60):  # Wrap lines at 60 characters
                f.write(f"{sequence[i:i+60]}\n")
    
    if debug:
        print(f"Wrote {len(sequences)} sequences to {output_path}")
        print(f"Total residues written: {sum(len(seq) for seq in sequences.values())}")


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
    longest_transcripts = get_longest_transcripts(sequences, df_mapping, debug=False)

    # Step 4: Write the longest isoforms to a new FASTA file
    write_fasta(longest_transcripts, output_fasta_path, debug=debug)

    if debug:
        print(f"Longest isoforms written to: {output_fasta_path}")



from Bio import SeqIO
import pandas as pd

def extract_protein_sequences(ref_transcript_id, target_gene_id, 
                          ref_gene_map, target_gene_map,
                          ref_sequences, target_sequences,
                          debug=False):
    """
    Extract amino acid sequences for a reference transcript and the longest protein 
    sequence from its orthologous target gene.
    
    Args:
        ref_transcript_id (str): Transcript ID from reference genome
        target_gene_id (str): Gene ID from target genome (found through orthology)
        ref_gene_map (pd.DataFrame): Gene mapping DataFrame for reference genome
        target_gene_map (pd.DataFrame): Gene mapping DataFrame for target genome
        ref_sequences (dict): Dictionary of reference protein sequences
        target_sequences (dict): Dictionary of target protein sequences
        debug (bool): Whether to print debug information
        
    Returns:
        tuple: (ref_sequence, target_sequence, target_protein_id)
            - ref_sequence: The amino acid sequence for the reference transcript
            - target_sequence: The longest protein sequence for the target gene
            - target_protein_id: The ID of the longest target protein sequence
    """
    # Get reference protein ID from transcript ID
    ref_protein_id = ref_gene_map[
        ref_gene_map['transcript_id'] == ref_transcript_id
    ]['protein_id'].iloc[0]
    
    # Get all target protein IDs for the target gene
    target_protein_ids = target_gene_map[
        target_gene_map['gene_id'] == target_gene_id
    ]['protein_id'].tolist()
    
    if debug:
        print(f"Reference protein ID: {ref_protein_id}")
        print(f"Target protein IDs: {target_protein_ids}")
    
    # Extract reference sequence
    ref_sequence = ref_sequences.get(ref_protein_id)
    if ref_sequence is None:
        raise KeyError(f"Reference protein {ref_protein_id} not found in sequences")
        
    # Find the longest target protein sequence
    longest_target_seq = ""
    longest_target_id = None
    for protein_id in target_protein_ids:
        if protein_id in target_sequences:
            seq = target_sequences[protein_id]
            if len(seq) > len(longest_target_seq):
                longest_target_seq = seq
                longest_target_id = protein_id
    
    if not longest_target_seq:
        raise KeyError(f"No target proteins found in sequences for gene {target_gene_id}")
        
    if debug:
        print(f"Selected longest target protein: {longest_target_id} "
              f"(length: {len(longest_target_seq)})")
        
    return ref_sequence, longest_target_seq, longest_target_id


def load_genome_data(genome_dir, genome_name):
    """Load protein sequences and create gene mapping for a genome."""
    logger.info(f"Loading genome data for {genome_name}")
    
    protein_path = os.path.join(genome_dir, f"{genome_name}_AnnotatedProteins.fasta")
    logger.debug(f"Reading protein sequences from {protein_path}")
    
    proteins = SeqIO.to_dict(SeqIO.parse(protein_path, "fasta"))
    protein_seqs = {k: str(v.seq) for k, v in proteins.items()}
    logger.debug(f"Loaded {len(protein_seqs)} protein sequences")
    
    gff_path = os.path.join(genome_dir, f"{genome_name}.gff")
    logger.debug(f"Reading GFF from {gff_path}")
    gene_map = create_gene_mapping(str(gff_path))
    logger.debug(f"Created gene mapping with {len(gene_map)} entries")
    
    return protein_seqs, gene_map