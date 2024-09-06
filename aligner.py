from tqdm import tqdm
from collections import defaultdict
from utils import reverse_complement
from Bio import SeqIO
import gzip

def open_reads_file(reads_file):
    """
    Automatically handle gzipped or plain FASTA files.

    Parameters:
    reads_file (str): Path to the gzipped or plain FASTA file.

    Returns:
    generator: SeqIO parser object.
    """
    if reads_file.endswith('.gz'):
        return gzip.open(reads_file, "rt", encoding='utf-8')
    else:
        return open(reads_file, "rt", encoding='utf-8')

def pseudoalignment(reads_file, index, k=31):
    """
    Perform pseudoalignment of reads to the k-mer index.

    Parameters:
    reads_file (str): Path to the gzipped or plain FASTA file with the reads.
    index (dict): K-mer index to match reads against.
    k (int): Length of k-mers to match.

    Returns:
    tuple: Counts of equivalence classes and the equivalence classes themselves.
    """
    equivalence_class_counts = defaultdict(int)
    equivalence_classes = defaultdict(set)

    with open_reads_file(reads_file) as f:
        records = list(SeqIO.parse(f, "fasta"))
        total_records = len(records)
        
        if total_records == 0:
            print("No records found in the reads file.")
            return equivalence_class_counts, equivalence_classes

        for record in tqdm(records, total=total_records, desc="Pseudoalignment"):
            seq = str(record.seq)
            seq_rc = reverse_complement(seq)
            current_kmer_eq_class = None
            current_kmer_rc_eq_class = None

            # Perform k-mer alignment and reverse complement alignment
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_rc = seq_rc[i:i+k]

                if kmer in index:
                    if current_kmer_eq_class is None:
                        current_kmer_eq_class = index[kmer]
                    else:
                        current_kmer_eq_class &= index[kmer]  # Efficient intersection

                if kmer_rc in index:
                    if current_kmer_rc_eq_class is None:
                        current_kmer_rc_eq_class = index[kmer_rc]
                    else:
                        current_kmer_rc_eq_class &= index[kmer_rc]

            # Assign forward and reverse equivalence class keys
            forward_ec_key = tuple(sorted(current_kmer_eq_class)) if current_kmer_eq_class else "NA"
            reverse_ec_key = tuple(sorted(current_kmer_rc_eq_class)) if current_kmer_rc_eq_class else "NA"

            # Update counts and equivalence classes
            if current_kmer_eq_class:
                equivalence_class_counts[forward_ec_key] += 1
                equivalence_classes[forward_ec_key] = current_kmer_eq_class

            if current_kmer_rc_eq_class:
                equivalence_class_counts[reverse_ec_key] += 1
                equivalence_classes[reverse_ec_key] = current_kmer_rc_eq_class

    return equivalence_class_counts, equivalence_classes
