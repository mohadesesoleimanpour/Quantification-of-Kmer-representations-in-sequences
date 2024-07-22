import itertools
import pandas as pd
from Bio import SeqIO
import os

"""
Created on Fri Jul 19 15:30:24 2024

@author: MSoleimanpour
"""

'''Title: Calculating and Comparing Observed and Expected Probabilities of DNA Motifs Using Markov Chains

This Python script calculates the observed and expected probabilities of DNA motifs in a given sequence using first-order 
and second-order Markov chains. The sequence is defined, and probabilities for single, di-, tri-, and tetranucleotides are calculated.
The script uses a function to count occurrences of k-mers and compute their probabilities.
Expected probabilities for di-, tri- and tetranucleotides are determined based on Markov chain principles: 
first-order (considering adjacent nucleotide pairs) and second-order (considering triplet contexts). 
Then compares observed and expected probabilities by computing their ratios, organizing the data into a DataFrame, 
and saving it to an Excel file, facilitating the identification of motifs with significant deviations from expected patterns.'''

# Create output directory if it doesn't exist
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Function to clean sequences
def clean_sequence(sequence):
    sequence = sequence.upper().replace('U', 'T')
    return sequence

# Step 1: Read sequences from FASTA file and concatenate them
sequences = []
with open('hiv-db.fasta', 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences.append(clean_sequence(str(record.seq)))

# Concatenate all sequences into a single string
sequence = ''.join(sequences)

# Function to calculate probabilities and counts with correct initialization
def calculate_probabilities_and_counts(sequence, k):
    counts = {''.join(kmer): 0 for kmer in itertools.product('ATCG', repeat=k)}

    total_kmers = 0
    removed_kmers = 0
    all_kmers = 0

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        all_kmers += 1
        if set(kmer) <= set('ATCG'):  # Check if kmer contains only A, T, C, G
            counts[kmer] += 1
            total_kmers += 1
        else:
            # print(kmer)

            removed_kmers += 1


        # print("Kmer", kmer)
        # print("counts", counts)

    # print("************************************")

    probabilities = {k: v / total_kmers for k, v in counts.items()} if total_kmers > 0 else {}

    # print("-------------------------")

    # Troubleshooting

    # print(f"For {k}-mers:")
    # print(f"Total number of {k}-mers (including invalid): {all_kmers}")
    # print(f"Number of valid {k}-mers (A, T, C, G only): {total_kmers}")
    # print(f"Number of removed {k}-mers: {removed_kmers}")
    # print(f"Number of unique {k}-mers observed: {len([v for v in counts.values() if v > 0])}")

    # Troubleshooting

    # # Check for "CGCG" regardless of k-mer size
    # cgcg_positions = [i for i in range(len(sequence) - 3) if sequence[i:i + 4] == "CGCG"]
    # if cgcg_positions:
    #     for pos in cgcg_positions:
    #         print(f"Found 'CGCG' at position {pos}, value: '{sequence[pos:pos + 20]}'")
    # else:
    #     print("'CGCG' not found in the sequence")
    # print("-------------------------")

    return probabilities, counts

# Calculate probabilities and counts for 1-mers to 4-mers
for k in range(1, 5):
    prob, count = calculate_probabilities_and_counts(sequence, k)

# Calculate probabilities and counts for 1-mers to 4-mers
prob_nucleotides, count_nucleotides = calculate_probabilities_and_counts(sequence, 1)
prob_dinucleotides, count_dinucleotides = calculate_probabilities_and_counts(sequence, 2)
prob_trinucleotides, count_trinucleotides = calculate_probabilities_and_counts(sequence, 3)
prob_tetranucleotides, count_tetranucleotides = calculate_probabilities_and_counts(sequence, 4)

# Expected probability of dinucleotides (independent probability)
expected_dinucleotides = {
    ''.join(kmer): prob_nucleotides.get(kmer[0], 0) * prob_nucleotides.get(kmer[1], 0)
    for kmer in itertools.product('ATCG', repeat=2)
}

# Expected probability of trinucleotides using first-order Markov chains
expected_trinucleotides_first_order = {}
for kmer in itertools.product('ATCG', repeat=3):
    kmer = ''.join(kmer)
    prob1 = prob_dinucleotides.get(kmer[:2], 0)
    prob2 = prob_dinucleotides.get(kmer[1:], 0)
    prob_middle = prob_nucleotides.get(kmer[1], 1)
    expected_trinucleotides_first_order[kmer] = (prob1 * prob2) / prob_middle if prob_middle != 0 else 0

# Expected probability of tetranucleotides using first-order Markov chains
expected_tetranucleotides_first_order = {}
for kmer in itertools.product('ATCG', repeat=4):
    kmer = ''.join(kmer)
    prob1 = prob_dinucleotides.get(kmer[:2], 0)
    prob2 = prob_dinucleotides.get(kmer[1:3], 0)
    prob3 = prob_dinucleotides.get(kmer[2:], 0)
    prob_middle1 = prob_nucleotides.get(kmer[1], 1)
    prob_middle2 = prob_nucleotides.get(kmer[2], 1)
    expected_tetranucleotides_first_order[kmer] = (prob1 * prob2 * prob3) / (prob_middle1 * prob_middle2) if (prob_middle1 * prob_middle2) != 0 else 0

# Expected probability of tetranucleotides using second-order Markov chains
expected_tetranucleotides_second_order = {}
for kmer in itertools.product('ATCG', repeat=4):
    kmer = ''.join(kmer)
    prob1 = prob_trinucleotides.get(kmer[:3], 0)
    prob2 = prob_trinucleotides.get(kmer[1:], 0)
    prob_middle = prob_dinucleotides.get(kmer[1:3], 1)
    expected_tetranucleotides_second_order[kmer] = (prob1 * prob2) / prob_middle if prob_middle != 0 else 0

# Calculate ratios and prepare DataFrame
data = []

# Add 1-mers
for kmer in itertools.product('ATCG', repeat=1):
    kmer = ''.join(kmer)
    observed_count = count_nucleotides.get(kmer, 0)
    observed_prob = prob_nucleotides.get(kmer, 0)
    expected_prob = 0.25
    ratio = observed_prob / expected_prob if expected_prob else 0
    data.append([kmer, observed_count, observed_prob, expected_prob, ratio, expected_prob, ratio])



# Add 2-mers
for kmer in itertools.product('ATCG', repeat=2):
    kmer = ''.join(kmer)
    observed_count = count_dinucleotides.get(kmer, 0)
    observed_prob = prob_dinucleotides.get(kmer, 0)
    expected_prob = expected_dinucleotides.get(kmer, 0)
    ratio = observed_prob / expected_prob if expected_prob else 0

    data.append([kmer, observed_count, observed_prob, expected_prob, ratio, expected_prob, ratio])



# Add 3-mers
for kmer in itertools.product('ATCG', repeat=3):
    kmer = ''.join(kmer)
    observed_count = count_trinucleotides.get(kmer, 0)
    observed_prob = prob_trinucleotides.get(kmer, 0)
    expected_prob = expected_trinucleotides_first_order.get(kmer, 0)
    ratio = observed_prob / expected_prob if expected_prob else 0

    data.append([kmer, observed_count, observed_prob, expected_prob, ratio, expected_prob, ratio])

# Add 4-mers
for kmer in itertools.product('ATCG', repeat=4):
    kmer = ''.join(kmer)
    observed_count = count_tetranucleotides.get(kmer, 0)
    observed_prob = prob_tetranucleotides.get(kmer, 0)
    expected_prob_first = expected_tetranucleotides_first_order.get(kmer, 0)
    expected_prob_second = expected_tetranucleotides_second_order.get(kmer, 0)
    ratio_first = observed_prob / expected_prob_first if expected_prob_first else 0
    ratio_second = observed_prob / expected_prob_second if expected_prob_second else 0
    data.append([kmer, observed_count, observed_prob, expected_prob_first, ratio_first, expected_prob_second, ratio_second])

import pandas as pd

columns = ["Motif", "Observed Count", "Observed Probability", "Expected Probability (First Order)", "Ratio (First Order)", "Expected Probability (Second Order)", "Ratio (Second Order)"]
df = pd.DataFrame(data, columns=columns)

# Split DataFrame into three different subsets
df_first_order = df[["Motif", "Ratio (First Order)"]]
df_second_order = df[["Motif", "Ratio (Second Order)"]]
df_all_columns = df

# Create output directory if it doesn't exist
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Save each DataFrame to a different Excel file
output_path_first_order = os.path.join(output_dir, "Kmer_Ratio_First_Order.xlsx")
output_path_second_order = os.path.join(output_dir, "Kmer_Ratio_Second_Order.xlsx")
output_path_all_columns = os.path.join(output_dir, "Motif_Probabilities_Information.xlsx")

df_first_order.to_excel(output_path_first_order, index=False)
df_second_order.to_excel(output_path_second_order, index=False)
df_all_columns.to_excel(output_path_all_columns, index=False)

print(f"Results saved to {output_dir}")