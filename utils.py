from Bio.Seq import Seq

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    
    Parameters:
    seq (str): DNA sequence to reverse complement.

    Returns:
    str: Reverse complement of the input sequence.
    """
    return str(Seq(seq).reverse_complement())
