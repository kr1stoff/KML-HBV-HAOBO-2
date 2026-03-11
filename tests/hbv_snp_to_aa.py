#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq

# HBV ORF. GENBANK 注释 S 基因就是 PreS1. HBVdb 中用 PreS1 不要用 S
ORFS = {
    "P": (2307, 1623),
    "S": (2848, 835),
    "C": (1901, 2452),
    "X": (1374, 1838),
}

GENOME_LEN = 3215


def circular_range(start, end):
    """Return genome positions in circular genome"""
    if start <= end:
        return list(range(start, end + 1))
    else:
        return list(range(start, GENOME_LEN + 1)) + list(range(1, end + 1))


def get_orf_seq(genome, start, end):
    pos = circular_range(start, end)
    seq = "".join(genome[p - 1] for p in pos)
    return seq, pos


def annotate_snp(genome, pos, alt):
    genome_mut = list(genome)
    genome_mut[pos - 1] = alt

    results = []

    for gene, (start, end) in ORFS.items():

        orf_seq, pos_map = get_orf_seq(genome, start, end)
        if pos not in pos_map:
            continue

        idx = pos_map.index(pos)

        codon_start = (idx // 3) * 3
        codon_ref = orf_seq[codon_start:codon_start + 3]

        orf_mut = list(orf_seq)
        orf_mut[idx] = alt
        codon_alt = "".join(orf_mut[codon_start:codon_start + 3])

        aa_ref = str(Seq(codon_ref).translate())
        aa_alt = str(Seq(codon_alt).translate())

        aa_pos = codon_start // 3 + 1

        if aa_ref != aa_alt:
            results.append(f"{gene}:{aa_ref}{aa_pos}{aa_alt}")

    return results


def main():

    ref = SeqIO.read(sys.argv[1], "fasta")
    genome = str(ref.seq).upper()

    vcf = open(sys.argv[2])

    print("POS\tREF\tALT\tAA_change")

    for line in vcf:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        pos = int(cols[1])
        ref = cols[3]
        alt = cols[4]

        if len(ref) != 1 or len(alt) != 1:
            continue

        aa = annotate_snp(genome, pos, alt)

        if aa:
            print(f"{pos}\t{ref}\t{alt}\t{';'.join(aa)}")


if __name__ == "__main__":
    main()
