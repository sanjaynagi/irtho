import warnings
import logging
import pandas as pd
from Bio import Align

from .utils import logger

def align_protein_sequences(ref_seq, target_seq, match_score=2, mismatch_score=-1, 
                          open_gap_score=-10, extend_gap_score=-0.5):
    """Perform pairwise alignment of protein sequences."""
    logger.debug("Performing sequence alignment")
    logger.debug(f"Reference sequence length: {len(ref_seq)}")
    logger.debug(f"Target sequence length: {len(target_seq)}")
    
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    
    alignments = aligner.align(ref_seq, target_seq)
    best_alignment = alignments[0]
    logger.debug(f"Alignment score: {best_alignment.score}")
    
    return best_alignment


def map_codon_position(ref_codon_num, alignment_str, window_size=10):
    """
    Maps a reference codon number to the corresponding codon number in the target sequence
    using a pairwise protein alignment. If target position is a gap, finds closest non-gap
    position within specified window and raises warning.
    
    Args:
        ref_codon_num (int): Codon number in reference sequence (1-based)
        alignment_str (str): String representation of the pairwise alignment
        window_size (int): Maximum distance to look for non-gap positions
        
    Returns:
        int: Corresponding codon number in target sequence (1-based)
        or None if no suitable position found within window
    """
    # Extract aligned sequences
    lines = alignment_str.split('\n')[::2]  # Skip match lines
    
    # Process the chunked alignment lines
    query_chunks = []
    target_chunks = []
    
    for line in lines:
        if len(line.split(None, 2)) == 3:
            if line.startswith('query'):
                query_chunks.append(line.split(None, 2)[2])
            elif line.startswith('target'):
                target_chunks.append(line.split(None, 2)[2])
            
    # Join the chunks
    query_seq = ''.join(query_chunks)   # Reference sequence
    target_seq = ''.join(target_chunks) # Target sequence
    
    if len(query_seq) != len(target_seq):
        raise ValueError("Aligned sequences must be the same length")
    
    # Convert codon number to amino acid position (1-based)
    aa_pos = ref_codon_num
    
    # Count non-gap positions up to our reference position
    ref_gaps = 0
    target_gaps = 0
    aligned_pos = 0
    
    # Find corresponding position in alignment
    for i in range(len(query_seq)):
        if query_seq[i] != '-':
            ref_gaps += 1
        if target_seq[i] != '-':
            target_gaps += 1
            
        if ref_gaps == aa_pos:
            aligned_pos = i
            break
            
    # If we found a valid position
    if ref_gaps == aa_pos:
        # Count non-gap positions in target up to aligned position
        target_pos = target_gaps
        
        # Check if target position is a gap
        if target_seq[aligned_pos] == '-':
            # Look for closest non-gap position within window
            left = max(0, aligned_pos - window_size)
            right = min(len(target_seq), aligned_pos + window_size + 1)
            
            # Count target positions up to left boundary
            left_target_count = sum(1 for x in target_seq[:left] if x != '-')
            
            closest_pos = None
            min_distance = window_size + 1
            best_target_pos = None
            
            # Check positions in both directions
            for i in range(left, right):
                if target_seq[i] != '-':
                    distance = abs(i - aligned_pos)
                    if distance < min_distance:
                        min_distance = distance
                        closest_pos = i
                        # Count non-gaps up to this position for target position
                        best_target_pos = left_target_count + sum(1 for x in target_seq[left:i+1] if x != '-')
            
            if closest_pos is not None:
                warnings.warn(
                    f"Target position was gap, using closest non-gap position {best_target_pos} "
                    f"({min_distance} positions away)",
                    UserWarning
                )
                return best_target_pos
                
            return None
            
        return target_pos
    else:
        raise ValueError(f"Reference codon number {ref_codon_num} is out of range")
    

def get_genomic_position_from_codon(codon_num, gff_path, protein_id, gene_map):
    """
    Convert a codon number to genomic coordinates using GFF file and protein ID.
    
    Args:
        codon_num (int): The codon number (1-based)
        gff_path (str): Path to the GFF file
        protein_id (str): The protein ID to look up
        gene_map (pd.DataFrame): Gene mapping DataFrame containing protein_id to transcript_id mappings
        
    Returns:
        tuple: (chromosome, start, end, strand)
            - chromosome: The chromosome name
            - start: Start position of the codon in genomic coordinates
            - end: End position of the codon in genomic coordinates
            - strand: The strand ('+' or '-')
    """
    # Get transcript ID for this protein
    transcript_id = gene_map[gene_map['protein_id'] == protein_id]['transcript_id'].iloc[0]
    
    # Parse GFF to get CDS entries for this transcript
    cds_regions = []
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            if fields[2] == 'CDS':
                attributes = dict(x.split('=') for x in fields[8].split(';'))
                # Handle both Parent and protein_id attributes
                matches = False
                if attributes.get('Parent') == transcript_id:
                    matches = True
                if attributes.get('protein_id') == protein_id:
                    matches = True
                if matches:
                    cds_regions.append({
                        'chrom': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'phase': int(fields[7]) if fields[7] != '.' else 0
                    })
    
    if not cds_regions:
        raise ValueError(f"No CDS regions found for protein {protein_id} (transcript {transcript_id})")
        
    # Sort CDS regions by genomic coordinates
    strand = cds_regions[0]['strand']
    chrom = cds_regions[0]['chrom']
    
    if strand == '+':
        cds_regions.sort(key=lambda x: x['start'])
    else:
        cds_regions.sort(key=lambda x: x['start'], reverse=True)
    
    # Convert codon number to nucleotide position (0-based)
    target_nt = (codon_num - 1) * 3
    
    # Account for phase in first CDS
    first_phase = cds_regions[0]['phase']
    if first_phase > 0:
        target_nt += first_phase
    
    # Find which CDS region contains our position
    current_pos = 0
    for cds in cds_regions:
        cds_length = cds['end'] - cds['start'] + 1
        
        if current_pos + cds_length > target_nt:
            # Found the CDS containing our position
            offset = target_nt - current_pos
            
            if strand == '+':
                codon_start = cds['start'] + offset
                codon_end = codon_start + 2  # +2 because codon is 3 bases
            else:
                codon_end = cds['end'] - offset
                codon_start = codon_end - 2  # -2 because codon is 3 bases
            
            return (chrom, codon_start, codon_end, strand)
            
        current_pos += cds_length
    
    raise ValueError(f"Codon number {codon_num} is outside CDS regions")


