import matplotlib.pyplot as plt
from collections import defaultdict

def plot_equivalence_class_stats(equivalence_class_counts, equivalence_classes):
    """
    Plot the equivalence class size distribution.

    Parameters:
    equivalence_class_counts (dict): Dictionary with equivalence class sizes as keys and their counts as values.
    equivalence_classes (dict): Dictionary with equivalence classes.
    """
    class_size_to_counts = defaultdict(int)

    for ec_key, count in equivalence_class_counts.items():
        class_size = len(equivalence_classes[ec_key])
        class_size_to_counts[class_size] += count

    sizes = list(class_size_to_counts.keys())
    frequencies = list(class_size_to_counts.values())

    plt.bar(sizes, frequencies, log=True)
    plt.xlabel('Equivalence Class Size (Number of Isoforms)')
    plt.ylabel('Frequency (Number of Reads)')
    plt.title('Equivalence Class Size Distribution')
    plt.savefig('outputs/equivalence_class_size_distribution.png')
    plt.show()
