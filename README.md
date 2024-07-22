# Quantification-of-Kmer-representations-in-sequences

This Python script calculates the observed and expected probabilities of DNA motifs in a given sequence using first-order 
and second-order Markov chains. The sequence is defined, and probabilities for single, di-, tri-, and tetranucleotides are calculated.
The script uses a function to count occurrences of k-mers and compute their probabilities.
Expected probabilities for di-, tri- and tetranucleotides are determined based on Markov chain principles: 
first-order (considering adjacent nucleotide pairs) and second-order (considering triplet contexts). 
Then compares observed and expected probabilities by computing their ratios, organizing the data into a DataFrame, 
and saving it to an Excel file, facilitating the identification of motifs with significant deviations from expected patterns.
