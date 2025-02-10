import os
import pandas as pd

from .utils import get_tqdm

tqdm = get_tqdm()

def evaluate_synteny(ref_gene, target_gene, ref_genes_df, target_genes_df, orthologs_df, target_species, window_size=5, n_required_orthologs=3):
    """
    Evaluate if two genes are in syntenic regions by checking orthology of neighboring genes.
    
    Args:
        ref_gene (str): Reference gene ID
        target_gene (str): Target gene ID
        ref_gff (str): Path to reference GFF file
        target_gff (str): Path to target GFF file
        orthologs_df (pd.DataFrame): DataFrame of orthologs from get_orthologs()
        target_species (str): Name of target species
        window_size (int): Number of genes to check on each side
        n_required_orthologs: n required orthologs
        
    Returns:
        tuple: (is_syntenic, shared_orthologs)
            - is_syntenic: Boolean indicating if region is syntenic
            - shared_orthologs: Number of shared orthologous pairs in window
    """
    
    # Get chromosome and position for reference gene
    ref_row = ref_genes_df[ref_genes_df['ID'] == ref_gene]
    if len(ref_row) == 0:
        return False, 0
    ref_chrom = ref_row['seqid'].iloc[0]
    
    # Get chromosome and position for target gene
    target_row = target_genes_df[target_genes_df['ID'] == target_gene]
    if len(target_row) == 0:
        return False, 0
    target_chrom = target_row['seqid'].iloc[0]
    
    # Get neighboring genes
    ref_neighbors = ref_genes_df[
        (ref_genes_df['seqid'] == ref_chrom) &
        (ref_genes_df['start'] >= ref_row['start'].iloc[0] - window_size * 10000) &
        (ref_genes_df['start'] <= ref_row['start'].iloc[0] + window_size * 10000)
    ]['ID'].tolist()
    
    target_neighbors = target_genes_df[
        (target_genes_df['seqid'] == target_chrom) &
        (target_genes_df['start'] >= target_row['start'].iloc[0] - window_size * 10000) &
        (target_genes_df['start'] <= target_row['start'].iloc[0] + window_size * 10000)
    ]['ID'].tolist()
    
    # Create sets of orthologous pairs
    ortholog_pairs = set()
    for _, row in orthologs_df.iterrows():
        ref_genes = set(row['gene'].split(', '))
        target_genes = set(row[target_species].split(', '))
        for r in ref_genes:
            for t in target_genes:
                ortholog_pairs.add((r, t))
    
    # Count shared orthologous pairs in window
    shared_orthologs = 0
    for ref_n in ref_neighbors:
        for target_n in target_neighbors:
            if (ref_n, target_n) in ortholog_pairs:
                shared_orthologs += 1
    
    # Define synteny threshold (you can adjust this)
    is_syntenic = shared_orthologs >= n_required_orthologs  # At least 3 shared orthologous pairs in window
    
    return is_syntenic, shared_orthologs

def add_synteny_information(targets_df, reference_dir, reference_species, target_species, ortho_obj):
    """
    Add synteny information to the targets DataFrame.
    
    Args:
        targets_df (pd.DataFrame): Input targets DataFrame
        reference_dir (str): Directory containing reference files
        reference_species (str): Name of reference species
        target_species (str): Name of target species
        ortho_obj (Orthologs): Orthologs object for getting ortholog relationships
        
    Returns:
        pd.DataFrame: DataFrame with additional 'syntenic' and 'shared_orthologs' columns
    """
    # Create a copy of the input DataFrame
    result_df = targets_df.copy()
    
    # Add columns
    result_df['syntenic'] = False
    result_df['shared_orthologs'] = 0
    
    # Get paths to GFF files
    ref_gff = os.path.join(reference_dir, f"{reference_species}.gff")
    target_gff = os.path.join(reference_dir, f"{target_species}.gff")

    # Load and process GFFs
    import allel
    ref_gff_df = allel.gff3_to_dataframe(ref_gff, attributes=['ID']).sort_values(['seqid', 'start'])
    target_gff_df = allel.gff3_to_dataframe(target_gff, attributes=['ID']).sort_values(['seqid', 'start'])

    # Filter GFFs based on source
    if 'RefSeq' in ref_gff_df['source'].unique():
        ref_genes_df = ref_gff_df.query("type == 'gene' and source == 'RefSeq'")
    else:
        ref_genes_df = ref_gff_df.query("type == 'protein_coding_gene'")
        
    if 'RefSeq' in target_gff_df['source'].unique():
        target_genes_df = target_gff_df.query("type == 'gene' and source == 'RefSeq'")
    else:
        target_genes_df = target_gff_df.query("type == 'protein_coding_gene'")
    
    # Get orthologs DataFrame
    orthologs_df = ortho_obj.get_orthologs(reference_species, target_species)
    
    # Evaluate synteny for each row
    for idx, row in tqdm(result_df.iterrows()):
        if pd.notna(row['gene']) and pd.notna(row[target_species]):
            is_syntenic, shared_count = evaluate_synteny(
                ref_gene=row['gene'],
                target_gene=row[target_species],
                ref_genes_df=ref_genes_df,
                target_genes_df=target_genes_df,
                orthologs_df=orthologs_df,
                target_species=target_species
            )
            result_df.at[idx, 'syntenic'] = is_syntenic
            result_df.at[idx, 'shared_orthologs'] = shared_count
    
    return result_df
