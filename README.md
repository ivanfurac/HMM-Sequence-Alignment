# HMM Sequence Alignment
This project was done as a part of the Bioinformatics 2 course at the Faculty of Electrical Engineering and Computing. The goal was to implement a pairwise nucleotide sequence alignment model based on hidden Markov model theory. 

## Project Overview
### Nucleotide Sequence Alignment
Nucleotide sequences (DNA and RNA) consist of a large number of nucleotides. There are four nucleotides that can appear in a DNA sequence, whose symbols are A, C, T, and G. In RNA sequences, nucleotide T is replaced with nucleotide U. Nucleotide sequences in bioinformatics can therefore be viewed as large strings, for example:

ATGAAAGTGAAGGAGATCATGAAGAATTATCGGCACTTATGGAA...

To compare two (or more) protein sequences (e.g., to determine regions of similarity, predict functional and structural similarities, or determine an evolutionary relationship), they need to be aligned. Aligning two or more sequences means inserting gap characters (gap character is "-") at appropriate positions in the sequences to achieve as many matches as possible. There are many computational algorithms based on dynamic programming that have been applied to the sequence alignment problem. These algorithms produce optimal alignments but suffer from large time and space complexities when dealing with longer sequences. The algorithm implemented in this project uses probabilities and hidden Markov model theory to produce alignments.

### Hidden Markov Models
Hidden Markov models include states and observations that are emitted in those states. The parameters of such a model are matrix π, which represents probabilities for starting from each state; matrix A, which represents probabilities for state transition; and matrix E, which represents probabilities that a certain symbol will be emitted in a certain state. In the case of a pairwise sequence alignment, possible states are match and insertion in the first (second sequence contains gap) or second sequence (first sequence contains gap), while symbols that can be emitted are nucleotide pairs (columns of an alignment). The picture below shows a hidden Markov model for pairwise nucleotide sequence alignment:

![image](https://github.com/ivanfurac/HMM-Sequence-Alignment/assets/73389887/778e1e01-20fd-4bdf-bdeb-b415d659eabd)

### Project Structure
There are four steps in training an HMM model for sequence alignment and checking its performance:
* **data preprocessing**: The HIV sequences dataset which was used in this project can be downloaded [here](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html).





## Repository and Usage

## Authors and Acknowledgment
Show your appreciation to those who have contributed to the project.
