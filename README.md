# HCR Probe Generator

The **HCR Probe Generator** is a full-featured, browser-based tool to design HCR v3.0 probe sets based on user-supplied cDNA sequences. It enables researchers and students to generate probes with optimal hybridization properties and biological relevance, with optional off-target filtering using BLAST.

---

## 🚀 Features

- Accepts any sense-strand cDNA sequence
- Select from multiple amplifier initiators (B1, B2, ..., B17)
- Filters out homopolymer-rich regions
- Computes Tm, GC%, and CDS coverage per probe
- Prioritizes probes within the coding sequence (CDS)
- Avoids overlapping probe pairs
- Scores and ranks candidate probes by biological and thermodynamic criteria
- Optional BLAST filtering against transcriptomic FASTA
- Generates:
  - IDT-compatible CSV for oPool ordering
  - Annotated GenBank file
  - Coverage plot
  - Probe summary table with all computed metrics

---

## 🧬 Input Parameters

| Field            | Description                                            |
| ---------------- | ------------------------------------------------------ |
| Gene name        | Any identifier for your gene                           |
| cDNA sequence    | Sense-strand DNA, no headers/spaces                    |
| Amplifier type   | Choose one of B1–B17                                   |
| Pause            | Offset from 5' before first probe is allowed           |
| PolyN thresholds | Max allowable length of poly-A/T or poly-C/G stretches |
| Tm / GC filters  | Set bounds and ideal target for probe evaluation       |
| Min spacing      | Minimum distance between probe pair starts             |

---

## 🧪 BLAST Filtering (Optional)

- Input any FASTA-formatted transcriptome or gene set
- BLASTn-short is used internally
- You can choose to:
  - See good/bad hits
  - Drop flagged probes with weak or ambiguous matches

---

## 📦 Outputs

- **IDT oPool CSV**: Split oligos with shared pool name
- **GenBank file**: Probe + CDS annotations for viewing in ApE
- **Interactive Plot**: Visual layout of all probe sets
- **Summary Table**: Tm, GC%, position, CDS status for each probe

---

## 🛠️ How to Run Locally

```bash
# Clone the repo
$ git clone https://github.com/Anneser/probe-generator
$ cd hcr-probe-generator

# Create virtual environment (recommended)
$ python -m venv venv
$ source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install dependencies
$ pip install -r requirements.txt

# Launch the app
$ streamlit run streamlit_app.py
```

---

## 📄 License & Citation

- Based on methodology from:
  - Choi et al., 2018 (HCR v3.0)
  - Wang et al., 2020 (Amplifier sets)

Please cite the appropriate literature and/or this repository when using the tool.

---

## 🤝 Contributing

Pull requests and feedback are welcome! Bug reports, UI suggestions, and feature requests can be filed via GitHub Issues.

---

## 👥 Credits

- Developed by: Lukas Anneser, Friedrich Lab
- Original codebase adapted from: Ryan W. Null, Ozpolat Lab



