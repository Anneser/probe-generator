import pandas as pd
import plotly.express as px
from blast_tools import find_longest_orf
import plotly.graph_objects as go
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
from Bio import SeqIO


def write_annotated_genbank(sense_seq, probe_data, orf_range, output_path):
    """
    Writes a GenBank file containing the input sequence annotated with:
    - Longest ORF
    - Each probe oligo (25 nt) as a misc_feature
    

    Parameters:
        sense_seq (str): Sense-strand cDNA sequence
        probe_data (dict): {index: (start_rc, sequence_with_nn, end_rc)} on antisense
        orf_range (tuple): (start, end) of longest ORF on sense strand
        output_path (str): Output .gb file path
    """
    seq_len = len(sense_seq)

    record = SeqRecord(
        Seq(sense_seq),
        id="HCR_ProbeSet",
        name="HCR_Probes",
        description="Annotated HCR probe binding regions + ORF",
        annotations={"molecule_type": "DNA", "date": str(datetime.today().date())}
    )

    # ORF annotation
    orf_start, orf_end = int(orf_range[0]), int(orf_range[1])
    orf_feature = SeqFeature(
        FeatureLocation(orf_start, orf_end),
        type="CDS",
        qualifiers={"label": "Longest ORF", "ApEinfo_color": "#66ff66"}
    )
    record.features.append(orf_feature)

    # Annotate each probe set
    for i, (start_rc, seq, end_rc) in probe_data.items():
        # Convert to sense-coordinates
        start_sense = seq_len - end_rc            # start of 5′ oligo
        mid_sense = seq_len - (start_rc + 27)     # start of 3′ oligo
        end_sense = seq_len - start_rc            # end of 3′ oligo

        # Oligo B (5′)
        feature_b = SeqFeature(
            FeatureLocation(int(start_sense), int(start_sense + 25)),
            type="misc_feature",
            qualifiers={
                "label": f"Probe {i+1}B",
                "ApEinfo_color": "#6495ED"  # Cornflower blue
            }
        )
        record.features.append(feature_b)

        # Spacer (2 nt)
        # feature_spacer = SeqFeature(
        #   FeatureLocation(int(start_sense + 25), int(start_sense + 27)),
        #    type="misc_feature",
        #    qualifiers={
        #        "label": f"Spacer {i+1}",
        #        "ApEinfo_color": "#CCCCCC"
        #    }
        # )
        # record.features.append(feature_spacer)

        # Oligo A (3′)
        feature_a = SeqFeature(
            FeatureLocation(int(mid_sense), int(end_sense)),
            type="misc_feature",
            qualifiers={
                "label": f"Probe {i+1}A",
                "ApEinfo_color": "#2E8B57"  # Sea green
            }
        )
        record.features.append(feature_a)

    # Write GenBank
    with open(output_path, "w") as handle:
        SeqIO.write(record, handle, "genbank")


def plot_probe_coverage(probe_data, cdna_len):
    """
    Generates a Plotly bar plot showing probe positions across the full cDNA.

    Parameters:
        probe_data (dict): {index: (start, sequence, end)}
        cdna_len (int): Length of the full reverse-complemented sequence

    Returns:
        plotly Figure
    """
    rows = []
    for i, (start, seq, end) in probe_data.items():
        rows.append({
            "Probe #": f"Probe {i+1}",
            "Start": start,
            "End": end,
            "Midpoint": (start + end) // 2,
            "Sequence": seq
        })

    df = pd.DataFrame(rows)

    fig = px.timeline(
        df,
        x_start="Start",
        x_end="End",
        y="Probe #",
        hover_data=["Sequence"],
        title="HCR Probe Coverage Across cDNA",
        height=500
    )
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(xaxis_title="cDNA Position (bp)", yaxis_title="Probes")
    return fig


def plot_probe_coverage_with_orf(probe_data, cdna_len, sense_seq):
    """
    Generates a bar plot with probes overlaid on top of the longest ORF.

    Parameters:
        probe_data (dict): {index: (start, sequence, end)}
        cdna_len (int): Length of the full sense cDNA sequence
        sense_seq (str): The original sense-strand cDNA sequence

    Returns:
        plotly Figure
    """
    orf_start, orf_end, orf_seq = find_longest_orf(sense_seq)

    fig = go.Figure()

    # ORF background (single bar)
    fig.add_trace(go.Bar(
        x=[orf_end - orf_start],
        y=["Longest ORF"],
        base=[orf_start],
        orientation='h',
        marker_color="rgba(0, 200, 0, 0.3)",
        hovertext=f"Longest ORF: {orf_start}–{orf_end} ({len(orf_seq)} bp)",
        name="Longest ORF"
    ))

    # Probes
    for i, (start, seq, end) in probe_data.items():
        fig.add_trace(go.Bar(
            x=[end - start],
            y=[f"Probe {i+1}"],
            base=start,
            orientation='h',
            name=f"Probe {i+1}",
            hovertext=seq,
            marker=dict(color="steelblue")
        ))

    fig.update_layout(
        title="🧬 HCR Probe Positions and Longest ORF",
        xaxis_title="cDNA Position (bp)",
        yaxis_title="",
        height=600,
        barmode='overlay',
        showlegend=False
    )

    return fig


def plot_probe_coverage_with_oligos(probe_data, cdna_len, sense_seq):
    """
    Plots each individual oligo (A and B) as separate color-coded bars on the cDNA axis.

    Parameters:
        probe_data (dict): {index: (start, combined_seq_with_nn, end)}
        cdna_len (int): Length of cDNA sequence
        sense_seq (str): Original sense-strand cDNA sequence

    Returns:
        plotly.graph_objects.Figure
    """
    fig = go.Figure()
    orf_start, orf_end, orf_seq = find_longest_orf(sense_seq)

    # ORF track
    fig.add_trace(go.Bar(
        x=[orf_end - orf_start],
        y=["ORF"],
        base=orf_start,
        orientation='h',
        marker_color="rgba(0, 200, 0, 0.2)",
        hovertext=f"Longest ORF: {orf_start}–{orf_end} ({len(orf_seq)} bp)",
        name="Longest ORF"
    ))

    for i, (start, seq, end) in probe_data.items():
        probe_5p = seq[:25]
        probe_3p = seq[-25:]

        # Oligo B: 5' half (probe_B)
        fig.add_trace(go.Bar(
            x=[25],
            y=[f"Set {i+1}B"],
            base=start,
            orientation='h',
            marker_color="cornflowerblue",
            name=f"Probe Set {i+1}B",
            hovertext=probe_5p
        ))

        # Oligo A: 3' half (probe_A)
        fig.add_trace(go.Bar(
            x=[25],
            y=[f"Set {i+1}A"],
            base=start + 27,
            orientation='h',
            marker_color="seagreen",
            name=f"Probe Set {i+1}A",
            hovertext=probe_3p
        ))

    fig.update_layout(
        title="🧬 Probe Oligo Coverage Across cDNA",
        xaxis_title="cDNA Position (bp)",
        yaxis_title="",
        barmode="overlay",
        height=700,
        showlegend=False,
    )

    return fig


def compute_probe_features(probe_data, orf_range, upinit, dninit, upspc, dnspc):
    """
    Computes GC content and Tm for each oligo in each probe set,
    using only the hybridizing region (25 nt), not initiator or spacer.

    Parameters:
        probe_data (dict): {index: (start, sequence_with_nn, end)}
        orf_range (tuple): (orf_start, orf_end)
        upinit, dninit (str): Initiator sequences
        upspc, dnspc (str): Spacer sequences

    Returns:
        pd.DataFrame: Table with one row per oligo
    """

    rows = []

    for i, (start, seq, end) in probe_data.items():
        probe_5p = seq[:25]
        probe_3p = seq[-25:]

        start_b = start
        end_b = start + 25
        start_a = start + 27
        end_a = end

        # Oligo B (5′ half)
        hybrid_seq_b = probe_5p
        gc_b = 100 * (hybrid_seq_b.count("G") + hybrid_seq_b.count("C")) / len(hybrid_seq_b)
        # tm_b = mt.Tm_Wallace(hybrid_seq_b)
        tm_b = mt.Tm_NN(hybrid_seq_b, nn_table=mt.DNA_NN4, Na=50, dnac1=50, dnac2=50)
        in_orf_b = start_b >= orf_range[0] and end_b <= orf_range[1]
        full_oligo_b = probe_5p + dnspc + dninit

        rows.append({
            "Probe ID": f"Probe {i+1}B",
            "Start": start_b,
            "End": end_b,
            "Length": 25,
            "GC (%)": round(gc_b, 1),
            "Tm (°C)": round(tm_b, 1),
            "Sequence": full_oligo_b,
            "In ORF": "✅" if in_orf_b else ""
        })

        # Oligo A (3′ half)
        hybrid_seq_a = probe_3p
        gc_a = 100 * (hybrid_seq_a.count("G") + hybrid_seq_a.count("C")) / len(hybrid_seq_a)
        # tm_a = mt.Tm_Wallace(hybrid_seq_a)
        tm_a = mt.Tm_NN(hybrid_seq_a, nn_table=mt.DNA_NN4, Na=50, dnac1=50, dnac2=50)
        in_orf_a = start_a >= orf_range[0] and end_a <= orf_range[1]
        full_oligo_a = upinit + upspc + probe_3p

        rows.append({
            "Probe ID": f"Probe {i+1}A",
            "Start": start_a,
            "End": end_a,
            "Length": 25,
            "GC (%)": round(gc_a, 1),
            "Tm (°C)": round(tm_a, 1),
            "Sequence": full_oligo_a,
            "In ORF": "✅" if in_orf_a else ""
        })

    return pd.DataFrame(rows)


def maximize_probe_set(valid_regions, spacing=2, cdna_len=None):
    """
    Given a list of homopolymer-free regions, select the longest set of
    non-overlapping 52-nt regions with at least `spacing` bp between them.

    Parameters:
        valid_regions (list of tuples): [(start, end), ...] coordinates of candidate probes
        spacing (int): minimum bp distance between adjacent probes
        cdna_len (int or None): optional, used for debugging or padding

    Returns:
        best_path (list of tuples): selected [(start, end)] pairs maximizing count
    """

    all_paths = []

    for x in range(len(valid_regions)):
        path = [valid_regions[x]]
        last_end = valid_regions[x][1]
        for y in range(x + 1, len(valid_regions)):
            if valid_regions[y][0] >= last_end + spacing:
                path.append(valid_regions[y])
                last_end = valid_regions[y][1]
        all_paths.append(path)

    # Select longest path
    longest = max(all_paths, key=len)
    return longest
