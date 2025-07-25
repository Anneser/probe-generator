Introduction
============

Welcome. If you want to visualize gene expression patterns using hybridization chain reaction (HCR), look no further.  
HCR is a powerful technique developed by Harry Choi during his time at Caltech:

- `Programmable in situ amplification for multiplexed imaging of mRNA expression <https://www.nature.com/articles/nbt.1692>`_ — HCR v1.0  
- `Next-Generation in Situ Hybridization Chain Reaction: Higher Gain, Lower Cost, Greater Durability <https://pubs.acs.org/doi/full/10.1021/nn405717p>`_ — HCR v2.0  
- `Third-generation in situ hybridization chain reaction: multiplexed, quantitative, sensitive, versatile, robust <https://journals.biologists.com/dev/article/145/12/dev165753/48466/Third-generation-in-situ-hybridization-chain>`_ — HCR v3.0

By now, his company is pushing the continuous development of the technique and the easiest way to get acquainted with the process is to get in touch with them and get one of the introduction kits. Up until generation v3.0, it was very common for researchers to 
create the probe sets binding to their genes of interest themselves. A great game-changer facilitating this approach was the `insitu-probe-generator <https://github.com/rwnull/insitu_probe_generator>`_ by the `Özpolat lab <https://bduyguozpolat.org/>`_.

Concept 
-------
For the longest time, visualization of gene expression in-situ was performed with either `chromogenic <https://en.wikipedia.org/wiki/Chromogenic_in_situ_hybridization>`_ or `fluorescent in-situ hybridization <https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization>`_. 
In a nutshell, you would clone a sequence characteristic for a gene of interest. These sequences could, for example, be tagged by an easily recognized molecule (biotin -streptavidin is a classic). First, the probes are added to a fixed tissue sample and hybridize to the sequence of
interest (if that sequence is expressed). In downstream steps, the tagged molecule is used to start a reaction leading to the amplification of signal, which can then be read out under a microscope. 

HCR invites you to tile your gene of interest with small stretches of probe sequences (50 nt). The probes contain a split inititator sequence both at the 5' and 3' end. This inititator sequence is necessary to break up the fluorescent hairpins h1 and h2 that subsequently form 
longer sequences at the bound mRNA of interest, leading to the accumulation of fluorescent probes and thus the build-up of signal indicating the expression of the gene of interest. The only things you need to make it work, are thus:

- a probe set specific to your gene of interest with the probes containing a specific initiator sequence
- the hairpin set recognizing this initiator sequence
- a set of buffers (see papers above)

All of this can be ordered via `molecular instruments <https://www.molecularinstruments.com/>`_. However, you can also design the probe set yourself using the tools provided here via IDT oPools. 

Tool Design 
-----------
As specified above, there already exists a nice tool to automatically create these probe sets, written by the Özpolat lab. However, this tool requires you to use the terminal or an IDE. Additionally, I always thought it would be nice to get some futher functionalities, such as

- upload of a fastafile with several genes or even just a list of refseq numbers instead of the sequence of only one gene of interest.
- some control over the sequence coverage and the probe set properties

For now (v0.1), I have implemented the basic probe set logic by the original probe generator and improved it a little bit, more updates to come.

- The sequence is tiled into all possible 52 nt windows, excluding windows containing user-specified homopolymer stretches (e.g., if you want to exclude GGGGGG, you can do that).
- Each window is now split into two 25 nt parts with a 2 nt spacer in the middle and basic properties such as GC content, Tm and whether the window is within the CDS or in the UTR are computed to create a combined score.
- A probe set is selected to maximize the number of probes (up to a user-specified maximum) and the combined score of the probes, while overlapping windows are excluded.
- The resulting probe set can be blasted against a user-provided fasta file to exclude potential off-targets
- The probe set can be downloaded as a .csv file that can be directly used for submission and ordering at IDT
- Additionally, the provided sequence can be downloaded as annotated .gb file. 