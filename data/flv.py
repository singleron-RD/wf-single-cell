import gzip

def reverse_complement(dna: str) -> str:
    """Returns the reverse complement of a DNA sequence, allowing 'N' bases."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(dna))

# Read barcodes
with open("/SGRNJ06/randd/USER/zhouyiqi/work/repo/CeleScope/celescope/data/chemistry/flv_rna/bc.txt", 'r') as f:
    bcs = f.read().splitlines()
    bcs = [bc.strip() for bc in bcs]

# Write joint reverse-complemented barcodes to gzip
with gzip.open('884k-flv.txt.gz', 'wt') as f:
    for bc1 in bcs:
        rc1 = reverse_complement(bc1)
        for bc2 in bcs:
            rc2 = reverse_complement(bc2)
            for bc3 in bcs:
                rc3 = reverse_complement(bc3)
                f.write(f"{rc1}{rc2}{rc3}\n")
