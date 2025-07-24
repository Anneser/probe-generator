def start():

    # These are required inputs
    name        = str(input("What is the gene name? (ex. eGFP) "))
    fullseq     = str(input("Enter the sense sequence of your cDNA without spaces or returns. "))
    fullseq = fullseq.upper()
    amplifier   = str(input("What is the amplifier to be used with this probe set? B1,B2,B3,B4,B5,B7,B9,B10,B11,B13,B14,B15,or B17  ").upper())
    if amplifier not in ['B1','B2','B3','B4','B5','B7','B9','B10','B11','B13','B14','B15','B17']:
        print("That choice was not recognized. Try again. ")
        amplifier   = str(input("What is the amplifier to be used with this probe set? B1,B2,B3,B4,B5,B7,B9,B10,B11,B13,B14,B15,or B17  ").upper())

    pause       = int(input("How many bases from 5' end of the Sense RNA before starting to hybridize? ex. 100 "))
    polyAT      = int(input("What is the max acceptable length for polyA or polyT homopolymers? "))
    polyCG      = int(input("What is the max acceptable length for polyC or polyG homopolymers? "))


    if (input("Do you want to choose program options?  Y -or- N ")).upper() == 'Y':

        # Optional inputs
        choose1      = input("Do you want to select between potential longest probe sets? (Y or N) ").lower()
        BlastProbes  = input("Would you like to BLAST potential probes against a FASTA file? (Y or N) ").lower()

        if BlastProbes == 'y':
            db       = input("Enter path to FASTA file to BLAST against (e.g. C:/users/user/file.fasta): ")
            show     = input("Show detailed BLAST output? (Y or N) ").lower()
            dropout  = input("Eliminate probes that appear in low-quality BLAST? (Y or N) ").lower()
        else:
            db = "n"
            show = "n"
            dropout = "n"

        report      = input("Display chosen parameters in output? (Y or N) ").lower()
        maxprobe    = input("Do you want to limit the number of probes made? (Y or N) ").lower()
        if maxprobe == 'y':
            numbr = int(input("Enter desired number of probes (default 33): ") or "33")
        else:
            numbr = 0

        min_spacing = int(input("Minimum spacing between probes (in bp, default = 2): ") or "2")

    else:
        choose1 = 'n'
        BlastProbes = 'n'
        db = "n"
        dropout = 'n'
        show = 'n'
        report = 'n'
        maxprobe = 'n'
        numbr = 0
        min_spacing = 2  # default

    return (name, fullseq, amplifier, pause, choose1, polyAT, polyCG,
        BlastProbes, db, dropout, show, report, maxprobe, numbr, min_spacing)



