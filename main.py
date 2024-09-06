from indexer import create_index
from aligner import pseudoalignment
from plotter import plot_equivalence_class_stats
import pandas as pd

def main():
    transcriptome_file = "data/chr11_transcriptome.fasta"
    reads_file = "data/reads.fasta.gz"
    k = 100

    # Create k-mer index
    index = create_index(transcriptome_file, k)

    # Perform pseudoalignment
    eq_class_counts, eq_classes = pseudoalignment(reads_file, index, k)

    # Save results to a TSV file
    data = []
    for eq_class, count in eq_class_counts.items():
        num_items = len(eq_classes[eq_class])
        data.append([count, num_items, eq_class])

    df = pd.DataFrame(data, columns=['counts', 'number of items in equivalence class', 'isoforms in equivalence class']).sort_values(by='counts', ascending=False)
    df.to_csv('outputs/equivalence_classes.tsv', index=False, sep='\t')

    # Print top 15 equivalence classes
    print("Top 15 Equivalence Class Counts:")
    top_15 = df.head(15)
    print(top_15.to_string(index=False))

    # Plot the equivalence class statistics
    plot_equivalence_class_stats(eq_class_counts, eq_classes)

if __name__ == "__main__":
    main()
