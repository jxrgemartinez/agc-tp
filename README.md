# OTU Calculation

You can find the complete description of the lab [here](https://docs.google.com/document/d/1qWNqPZ9Ecd-yZ5Hpl6n2zd7ZGtHPjf3yaW1ulKRdWnk/edit?usp=sharing).

## Introduction

The objective of this lab will be to calculate the OTUs obtained from a "mock" sequencing. We have only amplified bacteria (not fungi). 8 species are thus expected.

You will need to develop a program that performs full-length dereplication, chimeric sequence detection, and clustering based on a greedy algorithm ("Abundance Greedy Clustering").

## Installing Dependencies

You will use the nwalign3, pytest, and pylint libraries from Python:
```
pip3 install --user nwalign3 pytest pylint pytest-cov
```

## Usage

You will need to develop a Python 3 program that performs full-length dereplication, chimeric sequence detection, and clustering based on a greedy algorithm ("Abundance Greedy Clustering"). It will take the following arguments:

 -i, -amplicon_file file containing sequences in FASTA format
 -s, -minseqlen Minimum sequence length (optional - default value 400)
 -m, -mincount Minimum sequence count (optional - default value 10)
 -c, -chunk_size Size of sequence partitions (optional - default value 100)
 -k, -kmer_size Length of "kmers" (optional - default value 8)
 -o, -output_file output file with OTUs in FASTA format

 ## Tests

You will test your functions using the command pytest --cov=agc to be executed in the agc-tp/ folder. Due to this constraint, the function names will not be free. It will therefore be imperative to respect the names of the "imposed" functions, as well as their characteristics and parameters.
You will also check the syntactic quality of your program by executing the command: pylint agc.py

## Contact

If you have any questions, you can contact me by email: amine.ghozlane[at]pasteur.fr