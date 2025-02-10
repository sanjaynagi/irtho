import os
import pandas as pd

from  .utils  import (
    get_tqdm, 
    extract_protein_sequences,
    load_genome_data,
    get_random_exon_positions
)

from .align import logger, align_protein_sequences, get_genomic_position_from_codon, map_codon_position

tqdm = get_tqdm()

class Orthologs:
    """
    A class to parse and store relevant information from OrthoFinder output.
    """

    def __init__(self, results_dir, debug=False):
        """
        Initialize the Orthologs class.

        Args:
            results_dir (str): Path to the OrthoFinder results directory.
            debug (bool): Whether to print debug information.
        """
        self.results_dir = results_dir
        self.debug = debug
        self.hierarchical_orthogroups = None
        self.one_to_one_orthologs = None

        if debug:
            print(f"Initializing OrthoFinderResults for directory: {results_dir}")

        # Load relevant files
        file_path = os.path.join(
            self.results_dir, "Comparative_Genomics_Statistics", "OrthologuesStats_one-to-one.tsv"
        )
        self.stats = pd.read_csv(file_path, sep="\t", index_col=0)
        self.orthogroups = self._load_hierarchical_orthogroups()

    def _load_hierarchical_orthogroups(self):
        """
        Load the hierarchical orthogroups file (N0.tsv).
        """
        file_path = os.path.join(
            self.results_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv"
        )
        if os.path.exists(file_path):
            if self.debug:
                print(f"Loading hierarchical orthogroups from: {file_path}")
            self.hierarchical_orthogroups = pd.read_csv(file_path, sep="\t")
        else:
            if self.debug:
                print(f"File not found: {file_path}")
            self.hierarchical_orthogroups = None

    def get_orthologs(self, species_a, species_b):
        """
        Get orthologs (one-to-one, one-to-many, and many-to-many) between two species.

        Args:
            species_a (str): Name of the first species.
            species_b (str): Name of the second species.

        Returns:
            pd.DataFrame: A DataFrame containing orthologs between the two species.
                        Includes a column 'ortholog_type' to specify 'one-to-one',
                        'one-to-many', or 'many-to-many'.
        """
        # Define the path to the pairwise orthologs file
        orthologs_file = os.path.join(
            self.results_dir,
            "Orthologues",
            f"Orthologues_{species_a}",
            f"{species_a}__v__{species_b}.tsv",
        )

        # Check if the file exists
        if not os.path.exists(orthologs_file):
            raise FileNotFoundError(
                f"Orthologs file not found: {orthologs_file}. Ensure the species names are correct."
            )

        if self.debug:
            print(f"Loading orthologs from: {orthologs_file}")

        # Load the orthologs file
        orthologs_df = pd.read_csv(orthologs_file, sep="\t")

        # Ensure the required columns are present
        required_columns = [species_a, species_b, "Orthogroup"]
        if not all(col in orthologs_df.columns for col in required_columns):
            raise ValueError(
                f"Orthologs file does not contain the required columns: {required_columns}"
            )

        if self.debug:
            print(f"Orthologs file loaded with {len(orthologs_df)} rows.")

        # Identify one-to-one orthologs
        one_to_one = orthologs_df[
            (orthologs_df[species_a].str.contains(",") == False)  # No commas in species_a
            & (orthologs_df[species_b].str.contains(",") == False)  # No commas in species_b
        ].copy()
        one_to_one["ortholog_type"] = "one-to-one"

        if self.debug:
            print(f"Found {len(one_to_one)} one-to-one orthologs.")

        # Identify one-to-many orthologs
        one_to_many = orthologs_df[
            (orthologs_df[species_a].str.contains(",") == False)  # Single gene in species_a
            & (orthologs_df[species_b].str.contains(","))  # Multiple genes in species_b
        ].copy()
        one_to_many["ortholog_type"] = "one-to-many"

        if self.debug:
            print(f"Found {len(one_to_many)} one-to-many orthologs.")


        # Identify many-to-one orthologs
        many_to_one = orthologs_df[
            (orthologs_df[species_a].str.contains(","))  # Multiple genes in species_a
            & (orthologs_df[species_b].str.contains(",") == False)  # Single gene in species_b
        ].copy()
        many_to_one["ortholog_type"] = "many-to-one"

        if self.debug:
            print(f"Found {len(many_to_one)} many-to-one orthologs.")

        # Identify many-to-many orthologs
        many_to_many = orthologs_df[
            (orthologs_df[species_a].str.contains(","))  # Multiple genes in species_a
            & (orthologs_df[species_b].str.contains(","))  # Multiple genes in species_b
        ].copy()
        many_to_many["ortholog_type"] = "many-to-many"

        if self.debug:
            print(f"Found {len(many_to_many)} many-to-many orthologs.")

        # Combine all ortholog types into a single DataFrame
        all_orthologs = pd.concat([one_to_one, one_to_many, many_to_one, many_to_many], ignore_index=True)

        if self.debug:
            print(f"Total orthologs retrieved: {len(all_orthologs)}")

        return all_orthologs.rename(columns={species_a: "gene"})


    def get_hierarchical_orthogroup(self, orthogroup_id):
        """
        Get the genes in a specific hierarchical orthogroup.

        Args:
            orthogroup_id (str): The ID of the hierarchical orthogroup.

        Returns:
            pd.DataFrame: A DataFrame containing the genes in the specified orthogroup.
        """
        if self.hierarchical_orthogroups is None:
            raise ValueError("Hierarchical orthogroups file not loaded.")

        if self.debug:
            print(f"Retrieving hierarchical orthogroup: {orthogroup_id}")

        # Filter for the specified orthogroup
        orthogroup = self.hierarchical_orthogroups[
            self.hierarchical_orthogroups["HOG"] == orthogroup_id
        ]

        return orthogroup

    def list_species(self):
        """
        List all species in the OrthoFinder results.

        Returns:
            list: A list of species names.
        """
        if self.stats is None:
            raise ValueError("One-to-one orthologs statistics file not loaded.")

        return list(self.stats.columns)
    
    def map_input_genes_to_orthologs(
        self, 
        input_df,
        reference_species,
        target_species,
        debug=False,
    ):
        """
        Map genes from an input DataFrame to their orthologs in a target species.

        Args:
            input_df (pd.DataFrame): Input DataFrame with a 'gene' column containing gene IDs.
            reference_species (str): Name of the reference species.
            target_species (str): Name of the target species.
            gene_mapping_df (pd.DataFrame): DataFrame with columns ['gene_id', 'transcript_id', 'protein_id'].
            debug (bool): Whether to print debug information.

        Returns:
            pd.DataFrame: A DataFrame with the original input and additional columns:
                        'ortholog_gene' and 'ortholog_type'.
        """
        # Step 1: Validate species in OrthoFinder results
        available_species = self.list_species()
        if reference_species not in available_species or target_species not in available_species:
            raise ValueError(
                f"One or both species ({reference_species}, {target_species}) are not in the OrthoFinder results."
            )

        if debug:
            print(f"Reference species: {reference_species}")
            print(f"Target species: {target_species}")
            print(f"Available species: {available_species}")

        # Step 2: Retrieve orthologs (one-to-one, one-to-many, and many-to-many)
        orthologs = self.get_orthologs(reference_species, target_species)

        if debug:
            print(f"Retrieved {len(orthologs)} orthologs (one-to-one, one-to-many, and many-to-many).")

        # Step 3: Merge the orthologs back into the input DataFrame
        final_df = input_df.merge(orthologs, on="gene", how="left")

        if debug:
            print(f"Final result contains {len(final_df)} rows.")
            print(final_df.head())

        return final_df


    def find_orthologous_target(self, ref_transcript_id, target_ortholog_id, codon_num,
                                ref_proteins, target_proteins, ref_gene_map, target_gene_map, 
                                target_gff_path, alignment_params=None):
        """Locate the orthologous position in the target genome."""
        logger.info(f"Processing ortholog mapping for reference transcript {ref_transcript_id}")
        logger.info(f"Target ortholog: {target_ortholog_id}, Codon: {codon_num}")
        
        result = {
            'target_codon': None,
            'target_contig': None,
            'target_start': None,
            'target_end': None,
            'target_strand': None
        }
        
        # Check for required inputs
        if pd.isna(codon_num) or pd.isna(target_ortholog_id) or "random" in str(codon_num):
            logger.warning("Missing codon or ortholog info - skipping")
            return result
            
        # Get protein sequences
        logger.debug("Extracting protein sequences")
        ref_seq, target_seq, target_id = extract_protein_sequences(
            ref_transcript_id=ref_transcript_id,
            target_gene_id=target_ortholog_id,
            ref_gene_map=ref_gene_map,
            target_gene_map=target_gene_map,
            ref_sequences=ref_proteins,
            target_sequences=target_proteins
        )
        logger.debug(f"Target protein ID: {target_id}")
        
        # Perform alignment
        alignment_params = alignment_params or {}
        logger.debug(f"Alignment parameters: {alignment_params}")
        alignment = align_protein_sequences(target_seq, ref_seq, **alignment_params)
        
        # Map codon position
        target_codon = map_codon_position(int(codon_num), str(alignment))
        logger.debug(f"Mapped reference codon {codon_num} to target codon {target_codon}")
        
        if target_codon is not None:
            # Get genomic coordinates
            logger.debug("Getting genomic coordinates")
            contig, start, end, strand = get_genomic_position_from_codon(
                target_codon,
                target_gff_path,
                target_id,
                target_gene_map
            )
            
            result.update({
                'target_codon': target_codon,
                'target_contig': contig,
                'target_start': start,
                'target_end': end,
                'target_strand': strand
            })
            logger.info(f"Found target location: {contig}:{start}-{end} ({strand})")
        else:
            logger.warning("Could not map codon position")
        
        return result


    def find_orthologous_targets(self, targets_df, reference_dir, ref_genome, target_genome, 
                            transcript_col='transcript', codon_col='codon'):
        """Map resistance sites from reference to target genome."""
        logger.info(f"Starting resistance site mapping from {ref_genome} to {target_genome}")
        logger.info(f"Processing {len(targets_df)} sites")
        
        # Load genome data
        logger.info("Loading genome data")
        ref_proteins, ref_gene_map = load_genome_data(reference_dir, ref_genome)
        target_proteins, target_gene_map = load_genome_data(reference_dir, target_genome)
        
        # Get target GFF path
        target_gff = os.path.join(reference_dir, f"{target_genome}.gff")
        logger.debug(f"Target GFF path: {target_gff}")
        
        # Add target columns
        for col in ['target_codon', 'target_contig', 'target_start', 'target_end', 'target_strand']:
            targets_df[col] = None
        
        # Process each row
        mapped_count = 0
        for idx, row in tqdm(targets_df.iterrows()):
            logger.debug(f"Processing row {idx}")
            result = self.find_orthologous_target(
                ref_transcript_id=row[transcript_col],
                target_ortholog_id=row[target_genome],
                codon_num=row[codon_col],
                ref_proteins=ref_proteins,
                target_proteins=target_proteins,
                ref_gene_map=ref_gene_map,
                target_gene_map=target_gene_map,
                target_gff_path=str(target_gff)
            )
            
            # Update DataFrame
            for key, value in result.items():
                targets_df.at[idx, key] = value
            
            if result['target_codon'] is not None:
                mapped_count += 1
        
        logger.info(f"Mapping complete. Successfully mapped {mapped_count}/{len(targets_df)} sites")
        return targets_df
    
    def process_random_targets(self, targets_df, reference_dir, target_species):
        """
        Process targets with 'randomN' in the codon column and valid orthologs.
        
        Args:
            targets_df (pd.DataFrame): Input targets DataFrame
            reference_dir (str): Directory containing reference files
            target_species (str): Name of target species
            
        Returns:
            pd.DataFrame: Updated targets DataFrame with random positions
        """
        # Create a copy of the input DataFrame
        result_df = targets_df.copy()
        
        # Keep track of rows to remove and add
        rows_to_remove = []
        new_rows = []
        
        # Iterate through rows where codon contains 'random' AND has a valid ortholog
        for idx, row in result_df.iterrows():
            # Check both conditions:
            # 1. codon contains 'random'
            # 2. target species column has a valid gene ID (not NaN or empty)
            if (isinstance(row['codon'], str) and 'random' in str(row['codon']).lower() and 
                pd.notna(row[target_species]) and str(row[target_species]).strip()):
                
                try:
                    # Extract number of random positions needed
                    n = int(str(row['codon']).lower().replace('random', ''))
                except ValueError:
                    print(f"Warning: Could not parse number from codon value {row['codon']}")
                    continue
                
                # Get GFF path for target species
                target_gff = os.path.join(reference_dir, f"{target_species}.gff")
                
                # Get random positions for the target gene
                random_positions = get_random_exon_positions(
                    target_gff, 
                    row[target_species],  # Use the target gene ID
                    n
                )
                
                # If we found positions, mark row for removal and create new rows
                if random_positions:
                    # Mark original row for removal
                    rows_to_remove.append(idx)
                    
                    # Create new rows for each random position
                    for i, (contig, start, end, strand) in enumerate(random_positions):
                        new_row = row.copy()
                        new_row['target_contig'] = contig
                        new_row['target_start'] = start
                        new_row['target_end'] = end
                        new_row['target_strand'] = strand
                        new_row['codon'] = f"random{i+1}"
                        new_rows.append(new_row)
                else:
                    print(f"Warning: No exons found for gene {row[target_species]}")
        
        # Remove all marked rows at once
        result_df = result_df.drop(rows_to_remove)
        
        # Add all new rows at once
        if new_rows:
            result_df = pd.concat([result_df, pd.DataFrame(new_rows)], ignore_index=True)
        
        return result_df.sort_values(['contig', 'start'])





    def __repr__(self):
        """
        Provide a pretty string representation of the Orthologs object.
        """
        repr_str = "Orthologs Summary\n"
        repr_str += f"Results Directory: {self.results_dir}\n"
        repr_str += "-" * 40 + "\n"

        if self.stats is not None:
            repr_str += f"Orthologs Loaded: {self.stats.shape[0]} genes x {self.stats.shape[1]} species\n"
        else:
            repr_str += "Orthologs: Not Loaded\n"

        repr_str += "-" * 40 + "\n"
        repr_str += "Available Methods:\n"
        repr_str += " - list_species(): List all species in the results\n"
        repr_str += " - get_orthologs(species_a, species_b): Retrieve orthologs between two species\n"
        repr_str += " - get_hierarchical_orthogroup(orthogroup_id): Retrieve genes in a specific hierarchical orthogroup\n"

        return repr_str