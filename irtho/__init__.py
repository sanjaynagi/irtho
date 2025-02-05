
# src/irtho/__init__.py
from .orthofinder import run_orthofinder
from .utils import logger, get_tqdm, load_fasta, write_fasta, write_longest_isoforms, create_gene_mapping, get_longest_transcripts, split_one_to_many_orthologs
from .orthologs import Orthologs
import logging 

tqdm = get_tqdm()

def enable_debug():
    """Enable debug logging"""
    logging.basicConfig(level=logging.DEBUG, 
                       format='%(asctime)s - %(levelname)s - %(message)s')
    
def disable_debug():
    """Disable all logging"""
    logging.getLogger(__name__).setLevel(logging.CRITICAL + 1)

__version__ = "0.1.0"