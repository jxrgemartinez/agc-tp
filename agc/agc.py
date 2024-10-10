#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line

        if len(seq) >= minseqlen:
            yield seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequence_counts = {}

    # Count occurrences of each sequence
    for sequence in read_fasta(amplicon_file, minseqlen):
        sequence_counts[sequence] = sequence_counts.get(sequence, 0) + 1

    filtered_sequences = [
        (seq, count) for seq, count in sequence_counts.items() if count >= mincount
    ]

    filtered_sequences.sort(key=lambda x: x[1], reverse=True)
   
    for sequence, count in filtered_sequences:
        yield [sequence, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq1, seq2 = alignment_list
    identical = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
    return (identical / len(seq1)) * 100

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otus = []
    
    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        if not otus:
            otus.append([sequence, count])
            continue
        
        is_otu = True
        for otu_seq, _ in otus:
            alignment = nw.global_align(sequence, otu_seq, gap_open=-1, gap_extend=-1, 
                                        matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)
            if identity > 97:
                is_otu = False
                break
        
        if is_otu:
            otus.append([sequence, count])
    
    return otus


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w') as f:
        for i, (sequence, count) in enumerate(OTU_list, start=1):
            f.write(f">OTU_{i} occurrence:{count}\n")
            
            wrapped_sequence = textwrap.fill(sequence, width=80)
            f.write(f"{wrapped_sequence}\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    # Perform abundance greedy clustering
    print(f"Starting abundance greedy clustering on {args.amplicon_file}")
    otu_list = abundance_greedy_clustering(
        args.amplicon_file, 
        args.minseqlen, 
        args.mincount,
        chunk_size=0,
        kmer_size=0
    )
    
    # Write OTUs to the output file
    print(f"Writing {len(otu_list)} OTUs to {args.output_file}")
    write_OTU(otu_list, args.output_file)
    
    print("OTU clustering completed successfully.")

    # Print some statistics
    total_sequences = sum(count for _, count in otu_list)
    print(f"Total number of OTUs: {len(otu_list)}")
    print(f"Total number of sequences represented: {total_sequences}")
    
    if len(otu_list) > 0:
        max_count = max(count for _, count in otu_list)
        min_count = min(count for _, count in otu_list)
        avg_count = total_sequences / len(otu_list)
        print(f"Largest OTU size: {max_count}")
        print(f"Smallest OTU size: {min_count}")
        print(f"Average OTU size: {avg_count:.2f}")

if __name__ == '__main__':
    main()
