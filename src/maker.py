import numpy as np
import pandas as pd
from Bio.Seq import Seq
from datetime import date
from utils import maximize_probe_set
from blast_tools import blast_probes_to_db, find_longest_orf
from Bio.SeqUtils import MeltingTemp as mt


def get_amplifier(ampl):
    """Return amplifier initiator sequences and spacers."""
    amplifiers = {
        "B1":  ("aa", "ta", "GAGGAGGGCAGCAAACGG", "GAAGAGTCTTCCTTTACG"),  # double-checked
        "B2":  ("aa", "aa", "CCTCGTAAATCCTCATCA", "ATCATCCAGTAAACCGCC"),  # double-checked
        "B3":  ("tt", "tt", "GTCCCTGCCTCTATATCT", "CCACTCAACTTTAACCCG"),  # double-checked
        "B4":  ("aa", "at", "CCTCAACCTACCTCCAAC", "TCTCACCATATTCGCTTC"),  # double-checked
        "B5":  ("aa", "aa", "CTCACTCCCAATCTCTAT", "CTACCCTACAAATCCAAT"),  # double-checked
        "B7":  ("ww", "ww", "CTTCAACCTCCACCTACC", "TCCAATCCCTACCCTCAC"),
        "B9":  ("ww", "ww", "CACGTATCTACTCCACTC", "TCAGCACACTCCCAACCC"),
        "B10": ("ww", "ww", "CCTCAAGATACTCCTCTA", "CCTACTCGACTACCCTAG"),
        "B11": ("ww", "ww", "CGCTTAGATATCACTCCT", "ACGTCGACCACACTCATC"),
        "B13": ("ww", "ww", "AGGTAACGCCTTCCTGCT", "TTATGCTCAACATACAAC"),
        "B14": ("ww", "ww", "AATGTCAATAGCGAGCGA", "CCCTATATTTCTGCACAG"),
        "B15": ("ww", "ww", "CAGATTAACACACCACAA", "GGTATCTCGAACACTCTC"),
        "B17": ("ww", "ww", "CGATTGTTTGTTGTGGAC", "GCATGCTAATCGGATGAG"),
    }
    if ampl not in amplifiers:
        raise ValueError(f"Amplifier '{ampl}' not recognized.")
    return amplifiers[ampl]


def filter_homopolymers(seq, polyAT, polyCG):
    """Remove regions with homopolymer runs (A/T/C/G) beyond allowed limits."""
    hpA, hpT = "A"*(polyAT+1), "T"*(polyAT+1)
    hpC, hpG = "C"*(polyCG+1), "G"*(polyCG+1)

    cdna_len = len(seq)
    start = np.arange(0, cdna_len - 52)
    end = start + 52
    valid_regions = []

    for s, e in zip(start, end):
        subseq = seq[s:e]
        if not any(hp in subseq for hp in [hpA, hpT, hpC, hpG]):
            valid_regions.append((s, e))

    return valid_regions


def generate_probes(positions, pause, cdna_len, min_spacing=2):
    """
    Selects non-overlapping probes with a minimum spacing requirement.

    Parameters:
        positions (list of (start, end)): Candidate probe regions (52bp each).
        pause (int): 5' pause from start of cDNA (used to trim max region).
        cdna_len (int): Total length of cDNA.
        min_spacing (int): Minimum number of bases between probes (default = 2).

    Returns:
        list of (start, end): Filtered, non-overlapping probe positions.
    """
    max_pos = cdna_len - pause
    filtered = [p for p in positions if p[1] < max_pos]

    if not filtered:
        return []

    selected = [filtered[0]]
    last_end = filtered[0][1]

    for s, e in filtered[1:]:
        if s > last_end + min_spacing:
            selected.append((s, e))
            last_end = e

    return selected


def downsample_probes(probe_coords, maxprobe, num_requested, scores=None):
    """
    Uniformly downsample the list of probe coordinates to a desired number.

    Parameters:
        probe_coords (list of tuple): List of (start, end) probe coordinates.
        maxprobe (str): 'y' to enable limiting, 'n' to return full list.
        num_requested (int): Desired number of probes. Defaults to 33 if 0.
        scores (list of float, optional): Unused here, kept for API compatibility.

    Returns:
        list of tuple: Downsampled list of (start, end) coordinates.
    """
    if maxprobe != 'y' or not probe_coords:
        return probe_coords

    total_probes = len(probe_coords)
    num_to_keep = num_requested if num_requested > 0 else 33

    if num_to_keep >= total_probes:
        print(f"There were fewer than {num_to_keep} pairs. No downsampling applied.")
        return probe_coords

    # Select evenly spaced probe indices
    indices = np.linspace(0, total_probes - 1, num=num_to_keep, dtype=int)
    return [probe_coords[i] for i in indices]


def format_output(probe_data, fullseq, cdna_len, amplifier_id, pause, gene_name, upinit, dninit, upspc, dnspc):
    """
    Formats the final probe output for CSV export and downstream visualization.

    Parameters:
        probe_data (dict): {index: (start, sequence, end)}
        fullseq (str): Antisense (target) sequence
        cdna_len (int): Length of antisense cDNA
        amplifier_id (str): Amplifier tag (e.g. B1)
        pause (int): 5' pause offset
        gene_name (str): User-defined gene label
        upinit, dninit (str): Amplifier initiator sequences
        upspc, dnspc (str): Spacer sequences

    Returns:
        tuple:
            - List of IDT-compatible CSV lines
            - Antisense sequence (as string)
            - Reconstructed sense sequence (string)
            - Dict of probe halves for visualization/data
    """
    graphic = ['n'] * cdna_len
    submission_lines = {}
    pool_name = f"{amplifier_id}_{gene_name}_{len(probe_data)}_Dla{pause}"

    # Build per-oligo lines
    for i, (idx, (start, seq, end)) in enumerate(probe_data.items()):
        probe_5p = seq[:25]
        probe_3p = seq[-25:]

        # Save 3′ oligo (oligo A)
        oligo_A = upinit + upspc + probe_3p
        submission_lines[f"{i+1}_A"] = f"{pool_name},{oligo_A}"

        # Save 5′ oligo (oligo B)
        oligo_B = probe_5p + dnspc + dninit
        submission_lines[f"{i+1}_B"] = f"{pool_name},{oligo_B}"

        # Track positions for visualization
        graphic[start:end] = seq

    # Assemble final sequences
    g = ''.join(graphic)
    antisense = Seq(g)
    sense = antisense.reverse_complement()

    return list(submission_lines.values()), str(antisense), str(sense), probe_data


def select_probes_wisp(candidates, spacing=2, max_count=33):
    """
    Solves the Weighted Interval Scheduling Problem (WISP) for non-overlapping probe windows.

    Parameters:
        candidates (list): list of tuples (start, end, score)
        spacing (int): min spacing in bp between probes
        max_count (int): maximum number of probes to return

    Returns:
        List of (start, end) tuples representing optimal probe set.
    """
    import bisect

    # Sort candidates by end position
    candidates = sorted(candidates, key=lambda x: x[1])
    starts = [c[0] for c in candidates]
    ends = [c[1] for c in candidates]
    scores = [c[2] for c in candidates]
    n = len(candidates)

    # Compute p(j): index of last compatible probe before j
    p = [0] * n
    for j in range(n):
        i = bisect.bisect_right(ends, starts[j] - spacing) - 1
        p[j] = i if i >= 0 else -1

    # DP table: M[j] = max score using first j probes
    M = [0] * n
    decision = [False] * n

    for j in range(n):
        incl = scores[j] + (M[p[j]] if p[j] != -1 else 0)
        excl = M[j-1] if j > 0 else 0
        if incl > excl:
            M[j] = incl
            decision[j] = True
        else:
            M[j] = excl

    # Reconstruct solution
    selected = []

    def trace(j):
        while j >= 0:
            if decision[j]:
                selected.append((candidates[j][0], candidates[j][1]))
                j = p[j]
            else:
                j -= 1

    trace(n - 1)
    selected.reverse()

    return selected[:max_count]


def smart_probe_selector(
    fullseq,
    valid_regions,
    max_count,
    orf_range,
    spacing=2,
    tm_target=50.0,
    gc_target=50.0,
    tm_bounds=(40.0, 70.0),
    gc_bounds=(30.0, 70.0)
):
    """
    Rank and select non-overlapping probe regions based on thermodynamic and positional criteria,
    falling back to relaxed bounds or greedy if needed.

    Returns:
        List of (start, end) tuples selected as probe positions
    """

    seq_len = len(fullseq)
    orf_start_rc = seq_len - orf_range[1]
    orf_end_rc = seq_len - orf_range[0]

    def score_region(start, end):
        seq_5p = fullseq[start: start + 25]
        seq_3p = fullseq[start + 27: end]

        tm_5p = mt.Tm_NN(seq_5p, nn_table=mt.DNA_NN4, Na=50, dnac1=50, dnac2=50)
        tm_3p = mt.Tm_NN(seq_3p, nn_table=mt.DNA_NN4, Na=50, dnac1=50, dnac2=50)
        gc_5p = 100 * (seq_5p.count("G") + seq_5p.count("C")) / 25
        gc_3p = 100 * (seq_3p.count("G") + seq_3p.count("C")) / 25

        # Hard cutoffs
        if not (tm_bounds[0] <= tm_5p <= tm_bounds[1] and tm_bounds[0] <= tm_3p <= tm_bounds[1]):
            return -float('inf')
        if not (gc_bounds[0] <= gc_5p <= gc_bounds[1] and gc_bounds[0] <= gc_3p <= gc_bounds[1]):
            return -float('inf')

        tm_score = 1 - (abs(tm_5p - tm_target) + abs(tm_3p - tm_target)) / 40
        gc_score = 1 - (abs(gc_5p - gc_target) + abs(gc_3p - gc_target)) / 100
        orf_score = (1 if orf_start_rc <= start <= orf_end_rc else 0) + \
                    (1 if orf_start_rc <= end <= orf_end_rc else 0)
        
        return tm_score + gc_score + orf_score 

    def prepare_candidates():
        return [(s, e, score_region(s, e)) for s, e in valid_regions if score_region(s, e) > -float('inf')]

    def try_with_bounds(tm_b, gc_b):
        nonlocal tm_bounds, gc_bounds
        tm_bounds = tm_b
        gc_bounds = gc_b
        scored = prepare_candidates()
        return select_probes_wisp(scored, spacing=spacing, max_count=max_count), scored

    # Attempt 1
    selected, all_scored = try_with_bounds(tm_bounds, gc_bounds)
    if len(selected) >= int(0.8 * max_count):
        print("probe sets found with user-specified limits")
        pass
    else:
        # Attempt 2
        relaxed_tm = (tm_bounds[0] - 2, tm_bounds[1] + 2)
        relaxed_gc = (gc_bounds[0] - 5, gc_bounds[1] + 5)
        selected, all_scored = try_with_bounds(relaxed_tm, relaxed_gc)
        if len(selected) < int(0.8 * max_count):
            # Attempt 3
            relaxed_tm = (tm_bounds[0] - 5, tm_bounds[1] + 5)
            relaxed_gc = (gc_bounds[0] - 10, gc_bounds[1] + 10)
            selected, all_scored = try_with_bounds(relaxed_tm, relaxed_gc)

    # Fallback to greedy if too few
    if len(selected) < int(0.8 * max_count):
        print("Falling back to greedy selection due to low probe yield.")
        candidates = prepare_candidates()
        # candidates.sort(key=lambda x: x[2], reverse=True)
        selected = []
        last_end = -spacing
        for start, end, score in candidates:
            if start >= last_end + spacing:
                selected.append((start, end))
                last_end = end
            if len(selected) >= max_count:
                break
        score_map = {(s, e): sc for s, e, sc in candidates}

    else:
        score_map = {(s, e): sc for s, e, sc in all_scored}

    # Optional downsampling
    if len(selected) > max_count:
        score_map = [score_map[(s, e)] for s, e in selected]
        selected = downsample_probes(selected, maxprobe='y', num_requested=max_count, scores=score_map)

    return selected, score_map


def maker(
    name,
    fullseq,
    amplifier,
    pause,
    polyAT,
    polyCG,
    BlastProbes,
    db,
    dropout,
    show,
    report,
    maxprobe,
    numbr,
    min_spacing,
    tm_bounds,
    gc_bounds,
    tm_target,
    gc_target
):
    """
    Main HCR probe design function. Handles all steps from input to formatted output.

    Returns:
        submission_lines (list): IDT-ready CSV lines
        probe_data (dict): indexed probe metadata
    """

    # Convert to antisense (target)
    fullseq = str(Seq(fullseq).reverse_complement())
    cdna_len = len(fullseq)

    # Get amplifier sequences
    upspc, dnspc, upinit, dninit = get_amplifier(amplifier.upper())

    # Step 1: Remove poly-n runs
    valid_regions = []
    hpA, hpT = "A"*(polyAT+1), "T"*(polyAT+1)
    hpC, hpG = "C"*(polyCG+1), "G"*(polyCG+1)
    for s in range(0, cdna_len - 52):
        subseq = fullseq[s:s+52]
        if not any(hp in subseq for hp in [hpA, hpT, hpC, hpG]):
            valid_regions.append((s, s+52))

    # Step 2: Find longest ORF (sense strand)
    sense_seq = str(Seq(fullseq).reverse_complement())
    orf_start, orf_end, _ = find_longest_orf(sense_seq)

    # Step 3: Smart selection of probe regions with scoring
    selected_coords, score_map = smart_probe_selector(
        fullseq=fullseq,
        valid_regions=valid_regions,
        max_count=999,  # get many and prune later
        orf_range=(orf_start, orf_end),
        spacing=min_spacing,
        tm_bounds=tm_bounds,
        gc_bounds=gc_bounds,
        tm_target=tm_target,
        gc_target=gc_target
    )

    # Step 4: Downsample to user-specified max
    selected_coords = downsample_probes(
        selected_coords,
        maxprobe=maxprobe,
        num_requested=numbr,
        scores=score_map
    )

    # Step 5: Build probe data
    probe_data = {}
    for i, (start, end) in enumerate(selected_coords):
        seq = fullseq[start:start+25] + "nn" + fullseq[start+27:end]
        probe_data[i] = (start, seq, end)

    # Step 6: Optional BLAST filtering
    if BlastProbes == 'y':
        probe_seqs = {i: v[1] for i, v in probe_data.items()}
        accepted, rejected, good_df, bad_df = blast_probes_to_db(probe_seqs, db)

        if dropout == 'y':
            probe_data = {i: probe_data[i] for i in accepted}
    else:
        good_df = pd.DataFrame()
        bad_df = pd.DataFrame()

    # Step 7: Output formatting
    submission_lines, antisense_str, sense_str, probe_data = format_output(
        probe_data, fullseq, cdna_len, amplifier, pause, name, upinit, dninit, upspc, dnspc
    )

    # Step 8: Optional report
    if report == 'y':
        print("\nRun Summary:")
        print(f"Date: {date.today()}")
        print(f"Pause: {pause}")
        print(f"PolyA/T max: {polyAT}, PolyC/G max: {polyCG}")
        print(f"Amplifier: {amplifier}")
        print(f"BLAST: {BlastProbes}, Dropout: {dropout}, Max probes: {maxprobe} ({numbr})")
        print(f"Min spacing: {min_spacing} bp")
        print(f"Tm target: {tm_target} (range: {tm_bounds})")
        print(f"GC target: {gc_target} (range: {gc_bounds})")
        print(f"Total probe pairs: {len(probe_data)}")

    return submission_lines, probe_data, good_df, bad_df
