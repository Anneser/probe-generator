import streamlit as st
from src.maker import maker
from src.utils import plot_probe_coverage, plot_probe_coverage_with_orf, plot_probe_coverage_with_oligos
from src.utils import compute_probe_features
from src.blast_tools import find_longest_orf
from src.maker import get_amplifier
import os
from src.utils import write_annotated_genbank

st.set_page_config(page_title="HCR Probe Generator", layout="wide")

st.title("🔬 HCR Probe Generator")
st.markdown("Generate HCR v3.0 probe sets from your cDNA sequence with optional BLAST filtering.")

# Live-updating toggles (not inside the form)
BlastProbes = st.checkbox("BLAST probes against a database", value=False)
maxprobe = st.checkbox("Limit number of probes", value=True)

with st.form("probe_form"):
    st.header("🧬 Input Parameters")

    name = st.text_input("Gene name (e.g. eGFP)", "eGFP")

    fullseq = st.text_area(
        "Sense strand cDNA sequence",
        placeholder="Paste your sense sequence here (A/T/C/G only, no spaces or headers)",
        height=150
    ).upper().replace("\n", "").replace(" ", "")

    amplifier = st.selectbox("Amplifier type", ["B1", "B2", "B3", "B4", "B5", "B7", "B9", "B10", "B11", "B13", "B14", "B15", "B17"])

    pause = st.number_input("5' Pause before hybridization (in bases)", value=30, min_value=0)
    polyAT = st.number_input("Max polyA/T length", value=4, min_value=1)
    polyCG = st.number_input("Max polyC/G length", value=4, min_value=1)

    st.subheader("⚙️ Optional Settings")

    db = ""
    show = False
    dropout = False
    if BlastProbes:
        db = st.text_input("Path to local FASTA file", placeholder="e.g. /path/to/transcriptome.fasta")
        show = st.checkbox("Show detailed BLAST hits", value=False)
        dropout = st.checkbox("Drop low-quality BLAST hits", value=True)

    report = st.checkbox("Show parameter summary", value=True)

    numbr = 33
    if maxprobe:
        numbr = st.number_input("Number of probe pairs to generate", value=33, min_value=1)

    min_spacing = st.number_input("Minimum spacing between probe starts (in bp)", value=2, min_value=0)

    with st.expander("⚙️ Advanced Probe Scoring Options", expanded=False):
        tm_bounds = st.slider("Accepted Tm Range (°C)", 30, 80, (40, 60), key="tm_slider")
        gc_bounds = st.slider("Accepted GC Content Range (%)", 20, 80, (40, 60), key="gc_slider")
        tm_target = st.slider("Ideal Tm (°C)", 30, 80, 50, key="tm_target_slider")
        gc_target = st.slider("Ideal GC Content (%)", 20, 80, 50, key="gc_target_slider")

    submitted = st.form_submit_button("🚀 Run Probe Generator")

# Output section
if submitted:
    if not fullseq:
        st.error("You must enter a cDNA sequence.")
    else:
        st.success("Running probe generator...")
        submission_lines, probe_data, good_df, bad_df = maker(
            name=name,
            fullseq=fullseq,
            amplifier=amplifier,
            pause=pause,
            polyAT=polyAT,
            polyCG=polyCG,
            BlastProbes='y' if BlastProbes else 'n',
            db=db,
            dropout='y' if dropout else 'n',
            show='y' if show else 'n',
            report='y' if report else 'n',
            maxprobe='y' if maxprobe else 'n',
            numbr=numbr,
            min_spacing=min_spacing,
            tm_bounds=tm_bounds,
            gc_bounds=gc_bounds,
            tm_target=tm_target,
            gc_target=gc_target
        )
        if show and not good_df.empty:
            st.subheader("✅ Good BLAST Matches")
            st.dataframe(good_df)

        if show and not bad_df.empty:
            st.subheader("⚠️ Potential Off-Target Matches")
            st.dataframe(bad_df)

        # Convert submission_lines to CSV string
        csv_output = "Pool name,Sequence\n" + "\n".join(submission_lines)
        st.download_button(
            label="📥 Download IDT Submission CSV",
            data=csv_output,
            file_name=f"{name}_probes.csv",
            mime="text/csv"
        )

        orf_start, orf_end, _ = find_longest_orf(fullseq)
        gb_filename = f"{name}_annotated.gb"
        gb_path = os.path.join("outputs", gb_filename)

        os.makedirs("outputs", exist_ok=True)
        write_annotated_genbank(fullseq, probe_data, (orf_start, orf_end), gb_path)

        with open(gb_path, "rb") as f:
            st.download_button(
                label="📥 Download Annotated GenBank File",
                data=f.read(),
                file_name=gb_filename,
                mime="application/genbank"
            )

        # Visualize probe coverage
        fig = plot_probe_coverage_with_oligos(probe_data, len(fullseq), fullseq)
        with st.expander("🧬 View Probe Coverage Plot"):
            st.plotly_chart(fig, use_container_width=True)

        # Show interactive data table
        upspc, dnspc, upinit, dninit = get_amplifier(amplifier)
        probe_df = compute_probe_features(
            probe_data,
            (orf_start, orf_end),
            upinit,
            dninit,
            upspc,
            dnspc
        )

        st.subheader("📋 Probe Summary Table")
        st.dataframe(probe_df, use_container_width=True)
