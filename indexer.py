from Bio import SeqIO
from collections import defaultdict

def create_index(transcriptome_file, k=31):
    """
    Create an index of k-mers from the transcriptome.
    
    Parameters:
    transcriptome_file (str): Path to the transcriptome file in FASTA format.
    k (int): Length of k-mers to index.
    
    Returns:
    defaultdict: A dictionary with k-mers as keys and sets of transcript IDs as values.
    """
    index = defaultdict(set)
    for record in SeqIO.parse(transcriptome_file, "fasta"):
        seq = str(record.seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            index[kmer].add(record.id)  # k-mer as key and isoform list as value
    return index
