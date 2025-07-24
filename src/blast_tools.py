from Bio.Blast.Applications import NcbiblastnCommandline as bn
import pandas as pd
import io
import numpy as np
import os


def blast_probes_to_db(probe_seqs, db_path, e_thresh=1e-13, min_identity=75.0):
    """
    BLASTs probe sequences against a database and filters good matches.

    Parameters:
        probe_seqs (dict): {idx: sequence string}
        db_path (str): Path to FASTA database.
        e_thresh (float): Max e-value for keeping matches.
        min_identity (float): Minimum % identity required.

    Returns:
        accepted (list of int): Indices of accepted probes.
        rejected (list of int): Indices of flagged low-quality probes.
        blast_good (DataFrame), blast_bad (DataFrame): BLAST hit info
    """
    # Write temp FASTA file
    tmp_fa = "temp_probes.fa"
    with open(tmp_fa, "w") as f:
        for i, seq in probe_seqs.items():
            f.write(f">{i}\n{seq}\n")

    # Run BLAST
    cline = bn(query=tmp_fa, subject=db_path, outfmt=6, task='blastn-short')  
    stdout, stderr = cline()

    # Parse output
    dt = [
        (np.unicode_, 8), (np.unicode_, 40), (np.float32), (np.int32), (np.int32),
        (np.int32), (np.int32), (np.int32), (np.int32), (np.int32),
        (np.float32), (np.float32)
    ]
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    arr = np.genfromtxt(io.StringIO(stdout), delimiter='\t', dtype=dt)

    df = pd.DataFrame(arr, columns=cols)

    good_hits = df[(df["pident"] >= min_identity) & (df["evalue"] <= e_thresh)]
    bad_hits = df[(df["pident"] >= 60) & (df["evalue"] > 1e-12)]

    accepted_ids = list(pd.unique(good_hits["qseqid"].astype(int)))
    rejected_ids = list(pd.unique(bad_hits["qseqid"].astype(int)))

    if os.path.exists(tmp_fa):
        os.remove(tmp_fa)

    return accepted_ids, rejected_ids, good_hits, bad_hits


def find_longest_orf(dna_seq):
    """
    Identifies the longest ORF in all three reading frames of the sense strand.
    Only considers ATG as start and TAA, TAG, TGA as stop codons.

    Parameters:
        dna_seq (str): Sense-strand cDNA sequence.

    Returns:
        tuple: (orf_start, orf_end, orf_seq) in base-pair coordinates.
    """
    stop_codons = {"TAA", "TAG", "TGA"}
    max_orf = (0, 0, "")  # (start, end, sequence)

    for frame in range(3):
        i = frame
        while i < len(dna_seq) - 2:
            codon = dna_seq[i:i+3]
            if codon == "ATG":
                j = i + 3
                while j < len(dna_seq) - 2:
                    stop = dna_seq[j:j+3]
                    if stop in stop_codons:
                        length = j + 3 - i
                        if length > (max_orf[1] - max_orf[0]):
                            max_orf = (i, j + 3, dna_seq[i:j+3])
                        break
                    j += 3
                i = j  # Skip to after the stop codon
            else:
                i += 3

    return max_orf
